
#include <stdio.h>
#include "Macromolecule.h"
//#include "ResIni.h"
//#include "Atomo.h"




///Global variable to know if the aminoacid table is created
bool aminoacids_is_init=false;
///Global Variable. Collection of descriptions of all the possible aminoacids
t_aa AA[N_AMINO];
///Collection of Carbon-large VDW factors for all atoms at different distances
float VDW_C[N_ATOM_MAX][MAX_VDW];
///Collection of Carbon-large VDW factors for all atoms at different distances
float VDW_H[N_ATOM_MAX][MAX_VDW];
///Collection of electrostatic factors for all atoms at different distances
float ELE[MAX_ELE];


// inline ??
int resnum_from_resname( char * resname )
{

  int k;

//  for ( k = 0; k < 30; k++ )
  for ( k = 0; k < N_AMINO; k++ ) // Mon modified (1/12/2009) --> rURA = N_AMINO-1
  {
    if ( strcmp( resname, AA[k].aa_name3 ) == 0 )
      return k;
  }

//  if ( k >= 29 )
  if ( k >= N_AMINO-1 ) // Mon modified (1/12/2009) --> rURA = 39
    fprintf( stderr, "  Warning \"%s\" residue not identified set to GLY\n", resname );

  return GLY;

}

// Gets the residue name (in 3/4 characteres format). Watch out! You should allocate 5 chars (last one is '\0')
void resname_from_resnum(int num_aa, char *resname)
{

  switch (num_aa) {
   case 0:
      strcpy( resname, "ALA" );
      break;
   case 1:
      strcpy( resname, "CYS" );
      break;
   case 2:
      strcpy( resname, "ASP" );
      break;
   case 3:
      strcpy( resname, "GLU" );
      break;
   case 4:
      strcpy( resname, "PHE" );
      break;
   case 5:
      strcpy( resname, "GLY" );
      break;
   case 6:
      strcpy( resname, "HIS" );
      break;
   case 7:
      strcpy( resname, "ILE" );
      break;
   case 8:
      strcpy( resname, "LYS" );
      break;
   case 9:
      strcpy( resname, "LEU" );
      break;
   case 10:
      strcpy( resname, "MET" );
      break;
   case 11:
      strcpy( resname, "ASN" );
      break;
   case 12:
      strcpy( resname, "PRO" );
      break;
   case 13:
      strcpy( resname, "GLN" );
       break;
   case 14:
      strcpy( resname, "ARG" );
      break;
   case 15:
      strcpy( resname, "SER" );
      break;
   case 16:
      strcpy( resname, "THR" );
      break;
   case 17:
      strcpy( resname, "VAL" );
      break;
   case 18:
      strcpy( resname, "TRP" );
      break;
   case 19:
      strcpy( resname, "TYR" );
      break;
   case 20:
      strcpy( resname, "ASH" );
      break;
   case 21:
      strcpy( resname, "CYX" );
      break;
   case 22:
      strcpy( resname, "CYM" );
      break;
   case 23:
      strcpy( resname, "GLH" );
      break;
   case 24:
      strcpy( resname, "HIP" );
      break;
   case 25:
      strcpy( resname, "HID" );
      break;
   case 26:
      strcpy( resname, "HIE" );
      break;
   case 27:
      strcpy( resname, "LYN" );
      break;
   case 28:
      strcpy( resname, "TYM" );
      break;
   case 29:
      strcpy( resname, "MSE" );
      break;
   case 30:
      strcpy( resname, "NtE" );
      break;
   case 31:
      strcpy( resname, "CtE" );
      break;
   case 32:
      strcpy( resname, "DGUA" );
      break;
   case 33:
      strcpy( resname, "DADE" );
      break;
   case 34:
      strcpy( resname, "DCYT" );
      break;
   case 35:
      strcpy( resname, "DTHY" );
      break;
   case 36:
      strcpy( resname, "GUA" );
      break;
   case 37:
      strcpy( resname, "ADE" );
      break;
   case 38:
      strcpy( resname, "CYT" );
      break;
   case 39:
      strcpy( resname, "URA" );
      break;
  }
}



// index for fullatom_type_complete
// fullatom type number -> human readable abbreviation
//  1   // CNH2   CD
//  1   // CNH2   CG
//  2   // COO    CD
//  2   // COO    CG
//  3   // CH1    CB
//  3   // CH1    CG
//  4   // CH2    CB
//  4   // CH2    CD
//  4   // CH2    CE
//  4   // CH2    CG
//  4   // CH2    CG1
//  5   // CH3    CB
//  5   // CH3    CD1
//  5   // CH3    CD2
//  5   // CH3    CE
//  5   // CH3    CG1
//  5   // CH3    CG2
//  6   // aroC   CD1
//  6   // aroC   CD2
//  6   // aroC   CE1
//  6   // aroC   CE2
//  6   // aroC   CE3
//  6   // aroC   CG
//  6   // aroC   CH2
//  6   // aroC   CZ
//  6   // aroC   CZ2
//  6   // aroC   CZ3
//  7   // Ntrp   NE1
//  8   // Nhis   ND1
//  8   // Nhis   NE2
//  9   // NH2O   ND2
//  9   // NH2O   NE2
// 10   // Nlys   NZ
// 11   // Narg   NE
// 11   // Narg   NH1
// 11   // Narg   NH2
// 12   // Npro   N
// 13   // OH     OG
// 13   // OH     OG1
// 13   // OH     OH
// 14   // ONH2   OD1
// 14   // ONH2   OE1
// 15   // OOC    OD1
// 15   // OOC    OD2
// 15   // OOC    OE1
// 15   // OOC    OE2
// 16   // S      SD
// 16   // S      SG
// 17   // Nbb    N
// 18   // CAbb   CA
// 19   // CObb   CS
// 20   // OCbb   O
// 21   // Phos   P
// 22   // Hpol  1HD2
// 22   // Hpol  1HE2
// 22   // Hpol  1HH1
// 22   // Hpol  1HH2
// 22   // Hpol  1HZ
// 22   // Hpol  2HD2
// 22   // Hpol  2HE2
// 22   // Hpol  2HH1
// 22   // Hpol  2HH2
// 22   // Hpol  2HZ
// 22   // Hpol  3HZ
// 22   // Hpol   HE
// 22   // Hpol   HE1
// 22   // Hpol   HG
// 22   // Hpol   HG1
// 22   // Hpol   HH
// 23   // Hapo  1HB
// 23   // Hapo  1HD1
// 23   // Hapo  1HD2
// 23   // Hapo  1HE
// 23   // Hapo  1HG1
// 23   // Hapo  1HG2
// 23   // Hapo  2HA
// 23   // Hapo  2HB
// 23   // Hapo  2HD
// 23   // Hapo  2HD1
// 23   // Hapo  2HD2
// 23   // Hapo  2HE
// 23   // Hapo  2HG
// 23   // Hapo  2HG1
// 23   // Hapo  2HG2
// 23   // Hapo  3HA
// 23   // Hapo  3HB
// 23   // Hapo  3HD
// 23   // Hapo  3HD1
// 23   // Hapo  3HD2
// 23   // Hapo  3HE
// 23   // Hapo  3HG
// 23   // Hapo  3HG1
// 23   // Hapo  3HG2
// 23   // Hapo   HA
// 23   // Hapo   HB
// 23   // Hapo   HD2
// 23   // Hapo   HE2
// 23   // Hapo   HG
// 24   // Haro   HD1
// 24   // Haro   HD2
// 24   // Haro   HE1
// 24   // Haro   HE2
// 24   // Haro   HE3
// 24   // Haro   HH2
// 24   // Haro   HZ
// 24   // Haro   HZ2
// 24   // Haro   HZ3
// 25   // HNbb   H
// 26   // H2O


/**
* Initializes Alanine
*
* @param opt: Rosseta o ICM
*/

void init_NtE(Convention opt)
{

  //AMINOACIDO GLUTAMATO  ("N-terminal")
  int k, p = NtE;

  strcpy( AA[p].aa_name3, "NtE" );
  AA[p].aa_name1 = 'N';
  AA[p].mass = 128.11;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = true;
  AA[p].aa_is_nonpolar = false;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = true;

  AA[p].natoms = 17;
  AA[p].nheavyatoms = 9;
  AA[p].nchi = 3;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );

  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;

  AA[p].atom = ( t_aa_atom_p * ) malloc( AA[p].natoms * sizeof( t_aa_atom_p ) );


  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " CG " );
  strcpy( AA[p].atom[6].atom_name, " CD " );
  strcpy( AA[p].atom[7].atom_name, " OE1" );
  strcpy( AA[p].atom[8].atom_name, " OE2" );
  strcpy( AA[p].atom[9].atom_name, " H  " );
  strcpy( AA[p].atom[10].atom_name, " HA ");
  strcpy( AA[p].atom[11].atom_name, "2HB " );
  strcpy( AA[p].atom[12].atom_name, "3HB " );
  strcpy( AA[p].atom[13].atom_name, "2HG " );
  strcpy( AA[p].atom[14].atom_name, "3HG " );
  strcpy( AA[p].atom[15].atom_name,"1H  " );
  strcpy( AA[p].atom[16].atom_name,"2H  " );


  switch (opt) {
  	case EEF1:
  	  AA[p].atom[0].charge = -1.35; //    N
	  AA[p].atom[1].charge =  0.00; //    CA
	  AA[p].atom[2].charge =  0.55; //    C
	  AA[p].atom[3].charge = -0.55; //    O
	  AA[p].atom[4].charge =  0.00; //    CB
	  AA[p].atom[5].charge = -0.15; //    CG
	  AA[p].atom[6].charge =  1.35; //    CD
	  AA[p].atom[7].charge = -0.60; //    OE1
	  AA[p].atom[8].charge = -0.60; //    OE2
	  AA[p].atom[9].charge =  0.45; //    H
	  AA[p].atom[10].charge = 0.00; //    HA
	  AA[p].atom[11].charge = 0.00; //    2HB
	  AA[p].atom[12].charge = 0.00; //    3HB
	  AA[p].atom[13].charge = 0.00; //    2HG
	  AA[p].atom[14].charge = 0.00; //    3HG
	  AA[p].atom[15].charge = 0.45; //    1H
	  AA[p].atom[16].charge = 0.45; //    2H
  	  break;
  	default:
	  AA[p].atom[0].charge = -0.30; //    N
	  AA[p].atom[1].charge = 0.21; //    CA
	  AA[p].atom[2].charge = 0.51; //    C
	  AA[p].atom[3].charge = -0.51; //    O
	  AA[p].atom[4].charge = -0.18; //    CB
	  AA[p].atom[5].charge = -0.28; //    CG
	  AA[p].atom[6].charge = 0.62; //    CD
	  AA[p].atom[7].charge = -0.76; //    OE1
	  AA[p].atom[8].charge = -0.76; //    OE2
	  AA[p].atom[9].charge = 0.33; //    H
	  AA[p].atom[10].charge = 0.10; //    HA
	  AA[p].atom[11].charge = 0.09; //    2HB
	  AA[p].atom[12].charge = 0.09; //    3HB
	  AA[p].atom[13].charge = 0.09; //    2HG
	  AA[p].atom[14].charge = 0.09; //    3HG
	  AA[p].atom[15].charge = 0.33; //  1H
	  AA[p].atom[16].charge = 0.33; //  2H
	  break;
  }

  switch (opt) {
  	case EEF1:
	  AA[p].atom[0].fullatom_type =  14; // NH3    N
	  AA[p].atom[1].fullatom_type =   5; // CH1E   CA
	  AA[p].atom[2].fullatom_type =   4; // C      C
	  AA[p].atom[3].fullatom_type =  16; // O      O
	  AA[p].atom[4].fullatom_type =   6; // CH2E   CB
	  AA[p].atom[5].fullatom_type =   6; // CH2E   CG
	  AA[p].atom[6].fullatom_type =   4; // C      CD
	  AA[p].atom[7].fullatom_type =  17; // OC     OE1
	  AA[p].atom[8].fullatom_type =  17; // OC     OE2
	  AA[p].atom[9].fullatom_type =   1; // HC     1H
	  AA[p].atom[10].fullatom_type = -1; //     HA
	  AA[p].atom[11].fullatom_type = -1; //     2HB
	  AA[p].atom[12].fullatom_type = -1; //     3HB
	  AA[p].atom[13].fullatom_type = -1; //     2HG
	  AA[p].atom[14].fullatom_type = -1; // HA     3HG
	  AA[p].atom[15].fullatom_type =  1; // HC     2H
	  AA[p].atom[16].fullatom_type =  1; // HC     3H
	  break;
  	case Sybil:
  		  AA[p].atom[0].fullatom_type = 6; // Nbb    N
  		  AA[p].atom[1].fullatom_type = 2; // CAbb   CA
  		  AA[p].atom[2].fullatom_type = 1; // CObb   C
  		  AA[p].atom[3].fullatom_type = 8; // OCbb   O
  		  AA[p].atom[4].fullatom_type = 1; // CH2    CB
  		  AA[p].atom[5].fullatom_type = 1; // CH2    CG
  		  AA[p].atom[6].fullatom_type = 1; // COO    CD
  		  AA[p].atom[7].fullatom_type = 10; // OOC    OE1
  		  AA[p].atom[8].fullatom_type = 10; // OOC    OE2
  		  AA[p].atom[9].fullatom_type = 11; // HNbb   H
  		  AA[p].atom[10].fullatom_type = 14; // Hapo   HA
  		  AA[p].atom[11].fullatom_type = 14; // Hapo  2HB
  		  AA[p].atom[12].fullatom_type = 14; // Hapo  3HB
  		  AA[p].atom[13].fullatom_type = 14; // Hapo  2HG
  		  AA[p].atom[14].fullatom_type = 14; // Hapo  3HG
  		  AA[p].atom[15].fullatom_type = 13; // Hpol
  		  AA[p].atom[16].fullatom_type =13; // Hpol
  		  break;
  	default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 4; // CH2    CG
	  AA[p].atom[6].fullatom_type = 2; // COO    CD
	  AA[p].atom[7].fullatom_type = 15; // OOC    OE1
	  AA[p].atom[8].fullatom_type = 15; // OOC    OE2
	  AA[p].atom[9].fullatom_type = 25; // HNbb   H
	  AA[p].atom[10].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[11].fullatom_type = 23; // Hapo  2HB
	  AA[p].atom[12].fullatom_type = 23; // Hapo  3HB
	  AA[p].atom[13].fullatom_type = 23; // Hapo  2HG
	  AA[p].atom[14].fullatom_type = 23; // Hapo  3HG
	  AA[p].atom[15].fullatom_type = 22; // Hpol
	  AA[p].atom[16].fullatom_type = 22; // Hpol
	  break;
  }

  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 9; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 10; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--CG
  AA[p].atom[4].bonded_neighbor[2] = 11; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 12; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 4; // CG
  AA[p].atom[5].bonded_neighbor[0] = 4; // CG--CB
  AA[p].atom[5].bonded_neighbor[1] = 6; // CG--CD
  AA[p].atom[5].bonded_neighbor[2] = 13; // CG--2HG
  AA[p].atom[5].bonded_neighbor[3] = 14; // CG--3HG
  AA[p].atom[6].nbonded_neighbors = 3; // CD
  AA[p].atom[6].bonded_neighbor[0] = 5; // CD--CG
  AA[p].atom[6].bonded_neighbor[1] = 7; // CD--OE1
  AA[p].atom[6].bonded_neighbor[2] = 8; // CD--OE2
  AA[p].atom[7].nbonded_neighbors = 1; // OE1
  AA[p].atom[7].bonded_neighbor[0] = 6; // OE1--CD
  AA[p].atom[8].nbonded_neighbors = 1; // OE2
  AA[p].atom[8].bonded_neighbor[0] = 6; // OE2--CD
  AA[p].atom[9].nbonded_neighbors = 1; // H
  AA[p].atom[9].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[10].nbonded_neighbors = 1; // HA
  AA[p].atom[10].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[11].nbonded_neighbors = 1; //2HB
  AA[p].atom[11].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[12].nbonded_neighbors = 1; //3HB
  AA[p].atom[12].bonded_neighbor[0] = 4; //3HB--CB
  AA[p].atom[13].nbonded_neighbors = 1; //2HG
  AA[p].atom[13].bonded_neighbor[0] = 5; //2HG--CG
  AA[p].atom[14].nbonded_neighbors = 1; //3HG
  AA[p].atom[14].bonded_neighbor[0] = 5; //3HG--CG
  AA[p].atom[15].nbonded_neighbors = 1; // 1H
  AA[p].atom[15].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[16].nbonded_neighbors = 1; // 2H
  AA[p].atom[16].bonded_neighbor[0] = 0; // H-N


  for ( k = 0; k < 10; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 10; k < 15; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template for building  HA
  AA[p].atom[10].ta[0] = 1; //  CA
  AA[p].atom[10].ta[1] = 0; //  N
  AA[p].atom[10].ta[2] = 2; //  C

  //bk   template for building 2HB
  AA[p].atom[11].ta[0] = 4; //  CB
  AA[p].atom[11].ta[1] = 1; //  CA
  AA[p].atom[11].ta[2] = 5; //  CG

  //bk   template for building 3HB
  AA[p].atom[12].ta[0] = 4; //  CB
  AA[p].atom[12].ta[1] = 1; //  CA
  AA[p].atom[12].ta[2] = 5; //  CG

  //bk   template for building 2HG
  AA[p].atom[13].ta[0] = 5; //  CG
  AA[p].atom[13].ta[1] = 4; //  CB
  AA[p].atom[13].ta[2] = 6; //  CD

  //bk   template for building 3HG
  AA[p].atom[14].ta[0] = 5; //  CG
  AA[p].atom[14].ta[1] = 4; //  CB
  AA[p].atom[14].ta[2] = 6; //  CD


  //bk   chi angles required to build atoms GLU
  //bk   chi angles needed for building  CG
  AA[p].chi_required[0] [5] = true;

  //bk   chi angles needed for building  CD
  AA[p].chi_required[0] [6] = true;
  AA[p].chi_required[1] [6] = true;

  //bk   chi angles needed for building  OE1
  AA[p].chi_required[0] [7] = true;
  AA[p].chi_required[1] [7] = true;
  AA[p].chi_required[2] [7] = true;

  //bk   chi angles needed for building  OE2
  AA[p].chi_required[0] [8] = true;
  AA[p].chi_required[1] [8] = true;
  AA[p].chi_required[2] [8] = true;

  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [11] = true;

  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [12] = true;

  //bk   chi angles needed for building 2HG
  AA[p].chi_required[0] [13] = true;
  AA[p].chi_required[1] [13] = true;

  //bk   chi angles needed for building 3HG
  AA[p].chi_required[0] [14] = true;
  AA[p].chi_required[1] [14] = true;




  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   CG
  //bk   four atoms that define chi angle  2
  AA[p].chi_atoms[1] [0] = 1; //   CA
  AA[p].chi_atoms[1] [1] = 4; //   CB
  AA[p].chi_atoms[1] [2] = 5; //   CG
  AA[p].chi_atoms[1] [3] = 6; //   CD
  //bk   four atoms that define chi angle  3
  AA[p].chi_atoms[2] [0] = 4; //   CB
  AA[p].chi_atoms[2] [1] = 5; //   CG
  AA[p].chi_atoms[2] [2] = 6; //   CD
  AA[p].chi_atoms[2] [3] = 7; //   OE1


// Coordinates for "free" GLU. (From: /sbg/MON/ADP/fragments/GLU.db)

  AA[p].atom[0].icoor[0] = -0.000; //   N
  AA[p].atom[0].icoor[1] =  7.300; //   N
  AA[p].atom[0].icoor[2] =  12.666; //   N
  AA[p].atom[1].icoor[0] = -0.000; //   CA
  AA[p].atom[1].icoor[1] = 7.207; //   CA
  AA[p].atom[1].icoor[2] = 14.137; //   CA
  AA[p].atom[2].icoor[0] = -0.000; //   C
  AA[p].atom[2].icoor[1] = 8.613; //   C
  AA[p].atom[2].icoor[2] = 14.748; //   C
  AA[p].atom[3].icoor[0] = -0.958; //   O
  AA[p].atom[3].icoor[1] = 9.340; //   O
  AA[p].atom[3].icoor[2] = 14.399; //   O
  AA[p].atom[4].icoor[0] = 1.155; //   CB
  AA[p].atom[4].icoor[1] = 6.319; //   CB
  AA[p].atom[4].icoor[2] = 14.624; //   CB
  AA[p].atom[5].icoor[0] = 1.043; //   CG
  AA[p].atom[5].icoor[1] = 5.985; //   CG
  AA[p].atom[5].icoor[2] = 16.117; //   CG
  AA[p].atom[6].icoor[0] = 2.236; //   CD
  AA[p].atom[6].icoor[1] = 5.158; //   CD
  AA[p].atom[6].icoor[2] = 16.573; //   CD
  AA[p].atom[7].icoor[0] = 3.324; //   OE1
  AA[p].atom[7].icoor[1] = 5.762; //   OE1
  AA[p].atom[7].icoor[2] = 16.680; //   OE1
  AA[p].atom[8].icoor[0] = 2.038; //   OE2
  AA[p].atom[8].icoor[1] = 3.944; //   OE2
  AA[p].atom[8].icoor[2] = 16.792; //   OE2
  AA[p].atom[9].icoor[0] = -0.744; //   H
  AA[p].atom[9].icoor[1] = 7.932; //   H
  AA[p].atom[9].icoor[2] = 12.399; //   H
  AA[p].atom[10].icoor[0] = -0.935; //   HA
  AA[p].atom[10].icoor[1] = 6.727; //   HA
  AA[p].atom[10].icoor[2] = 14.427; //   HA
  AA[p].atom[11].icoor[0] = 1.134; //  2HB
  AA[p].atom[11].icoor[1] = 5.378; //  2HB
  AA[p].atom[11].icoor[2] = 14.073; //  2HB
  AA[p].atom[12].icoor[0] = 2.115; //  3HB
  AA[p].atom[12].icoor[1] = 6.800; //  3HB
  AA[p].atom[12].icoor[2] = 14.429; //  3HB
  AA[p].atom[13].icoor[0] = 1.023; //  2HG
  AA[p].atom[13].icoor[1] = 6.890; //  2HG
  AA[p].atom[13].icoor[2] = 16.722; //  2HG
  AA[p].atom[14].icoor[0] = 0.130; //  3HG
  AA[p].atom[14].icoor[1] = 5.418; //  3HG
  AA[p].atom[14].icoor[2] = 16.301; //  3HG
  AA[p].atom[15].icoor[0] = 0.885; //  1H
  AA[p].atom[15].icoor[1] = 7.670; //  1H
  AA[p].atom[15].icoor[2] = 12.349; //  1H
  AA[p].atom[16].icoor[0] = -0.155; //  2H
  AA[p].atom[16].icoor[1] = 6.388; //  2H
  AA[p].atom[16].icoor[2] = 12.261; //  2H
}

//FIN DE AMINOACIDO GLUTAMATO N-Terminal

void init_CtE(Convention opt)
{

  //AMINOACIDO GLUTAMATO  ("C-terminal")
  int k, p = CtE;


  strcpy( AA[p].aa_name3, "CtE" );
  AA[p].aa_name1 = 'C';
  AA[p].mass = 128.11;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = true;
  AA[p].aa_is_nonpolar = false;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = true;

  AA[p].natoms = 16;
  AA[p].nheavyatoms = 10;
  AA[p].nchi = 3;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );

  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;

  AA[p].atom = ( t_aa_atom_p * ) malloc( AA[p].natoms * sizeof( t_aa_atom_p ) );

  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " CG " );
  strcpy( AA[p].atom[6].atom_name, " CD " );
  strcpy( AA[p].atom[7].atom_name, " OE1" );
  strcpy( AA[p].atom[8].atom_name, " OE2" );
  strcpy( AA[p].atom[9].atom_name, " H  " );
  strcpy( AA[p].atom[10].atom_name, " HA " );
  strcpy( AA[p].atom[11].atom_name, "2HB " );
  strcpy( AA[p].atom[12].atom_name, "3HB " );
  strcpy( AA[p].atom[13].atom_name, "2HG " );
  strcpy( AA[p].atom[14].atom_name, "3HG " );
  strcpy( AA[p].atom[15].atom_name," OXT" );
  switch (opt) {
  	case EEF1:
  	  AA[p].atom[0].charge =  -0.35; //    N
	  AA[p].atom[1].charge =   0.10; //    CA
	  AA[p].atom[2].charge =   1.20; //    C
	  AA[p].atom[3].charge =  -0.60; //    O
	  AA[p].atom[4].charge =   0.00; //    CB
	  AA[p].atom[5].charge =  -0.15; //    CG
	  AA[p].atom[6].charge =   1.35; //    CD
	  AA[p].atom[7].charge =  -0.60; //    OE1
	  AA[p].atom[8].charge =  -0.60; //    OE2
	  AA[p].atom[9].charge =   0.25; //    H
	  AA[p].atom[10].charge =  0.00; //    HA
	  AA[p].atom[11].charge =  0.00; //    2HB
	  AA[p].atom[12].charge =  0.00; //    3HB
	  AA[p].atom[13].charge =  0.00; //    2HG
	  AA[p].atom[14].charge =  0.00; //    3HG
	  AA[p].atom[15].charge = -0.60; //    0XT
  	  break;
  	default:
	  AA[p].atom[0].charge = -0.47; //    N
	  AA[p].atom[1].charge = 0.07; //    CA
	  AA[p].atom[2].charge = 0.34; //    C
	  AA[p].atom[3].charge = -0.67; //    O
	  AA[p].atom[4].charge = -0.18; //    CB
	  AA[p].atom[5].charge = -0.28; //    CG
	  AA[p].atom[6].charge = 0.62; //    CD
	  AA[p].atom[7].charge = -0.76; //    OE1
	  AA[p].atom[8].charge = -0.76; //    OE2
	  AA[p].atom[9].charge = 0.31; //    H
	  AA[p].atom[10].charge = 0.09; //    HA
	  AA[p].atom[11].charge = 0.09; //    2HB
	  AA[p].atom[12].charge = 0.09; //    3HB
	  AA[p].atom[13].charge = 0.09; //    2HG
	  AA[p].atom[14].charge = 0.09; //    3HG
	  AA[p].atom[15].charge = -0.67; //   0XT
	  break;
  }

  switch (opt) {
  	case EEF1:
	  AA[p].atom[0].fullatom_type =   12; // NH1    N
	  AA[p].atom[1].fullatom_type =    5; // CH1E   CA
	  AA[p].atom[2].fullatom_type =    4; // C      C
	  AA[p].atom[3].fullatom_type =   17; // OC     O
	  AA[p].atom[4].fullatom_type =    6; // CH2E   CB
	  AA[p].atom[5].fullatom_type =    6; // CH2E   CG
	  AA[p].atom[6].fullatom_type =    4; // C      CD
	  AA[p].atom[7].fullatom_type =   17; // OC     OE1
	  AA[p].atom[8].fullatom_type =   17; // OC     OE2
	  AA[p].atom[9].fullatom_type =    0; // H      1H
	  AA[p].atom[10].fullatom_type =  -1; //      HA
	  AA[p].atom[11].fullatom_type =  -1; //      2HB
	  AA[p].atom[12].fullatom_type =  -1; //      3HB
	  AA[p].atom[13].fullatom_type =  -1; //      2HG
	  AA[p].atom[14].fullatom_type =  -1; //      3HG
	  AA[p].atom[15].fullatom_type =  17; // OC	    OXT
	  break;
  	case Sybil:
  	  AA[p].atom[0].fullatom_type = 6; // Nbb    N
  	  AA[p].atom[1].fullatom_type = 2; // CAbb   CA
  	  AA[p].atom[2].fullatom_type = 1; // CObb   C
  	  AA[p].atom[3].fullatom_type = 8; // OCbb   O
  	  AA[p].atom[4].fullatom_type = 1; // CH2    CB
  	  AA[p].atom[5].fullatom_type = 1; // CH2    CG
  	  AA[p].atom[6].fullatom_type = 1; // COO    CD
  	  AA[p].atom[7].fullatom_type = 10; // OOC    OE1
  	  AA[p].atom[8].fullatom_type = 10; // OOC    OE2
  	  AA[p].atom[9].fullatom_type = 11; // HNbb   H
  	  AA[p].atom[10].fullatom_type = 14; // Hapo   HA
  	  AA[p].atom[11].fullatom_type = 14; // Hapo  2HB
  	  AA[p].atom[12].fullatom_type = 14; // Hapo  3HB
  	  AA[p].atom[13].fullatom_type = 14; // Hapo  2HG
  	  AA[p].atom[14].fullatom_type = 14; // Hapo  3HG
  	  AA[p].atom[15].fullatom_type = 14; // OXT
  	  break;
  	default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 4; // CH2    CG
	  AA[p].atom[6].fullatom_type = 2; // COO    CD
	  AA[p].atom[7].fullatom_type = 15; // OOC    OE1
	  AA[p].atom[8].fullatom_type = 15; // OOC    OE2
	  AA[p].atom[9].fullatom_type = 25; // HNbb   H
	  AA[p].atom[10].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[11].fullatom_type = 23; // Hapo  2HB
	  AA[p].atom[12].fullatom_type = 23; // Hapo  3HB
	  AA[p].atom[13].fullatom_type = 23; // Hapo  2HG
	  AA[p].atom[14].fullatom_type = 23; // Hapo  3HG
	  AA[p].atom[15].fullatom_type = 20; // OXT
	  break;
  }


  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 9; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 10; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--CG
  AA[p].atom[4].bonded_neighbor[2] = 11; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 12; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 4; // CG
  AA[p].atom[5].bonded_neighbor[0] = 4; // CG--CB
  AA[p].atom[5].bonded_neighbor[1] = 6; // CG--CD
  AA[p].atom[5].bonded_neighbor[2] = 13; // CG--2HG
  AA[p].atom[5].bonded_neighbor[3] = 14; // CG--3HG
  AA[p].atom[6].nbonded_neighbors = 3; // CD
  AA[p].atom[6].bonded_neighbor[0] = 5; // CD--CG
  AA[p].atom[6].bonded_neighbor[1] = 7; // CD--OE1
  AA[p].atom[6].bonded_neighbor[2] = 8; // CD--OE2
  AA[p].atom[7].nbonded_neighbors = 1; // OE1
  AA[p].atom[7].bonded_neighbor[0] = 6; // OE1--CD
  AA[p].atom[8].nbonded_neighbors = 1; // OE2
  AA[p].atom[8].bonded_neighbor[0] = 6; // OE2--CD
  AA[p].atom[9].nbonded_neighbors = 1; // H
  AA[p].atom[9].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[10].nbonded_neighbors = 1; // HA
  AA[p].atom[10].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[11].nbonded_neighbors = 1; //2HB
  AA[p].atom[11].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[12].nbonded_neighbors = 1; //3HB
  AA[p].atom[12].bonded_neighbor[0] = 4; //3HB--CB
  AA[p].atom[13].nbonded_neighbors = 1; //2HG
  AA[p].atom[13].bonded_neighbor[0] = 5; //2HG--CG
  AA[p].atom[14].nbonded_neighbors = 1; //3HG
  AA[p].atom[14].bonded_neighbor[0] = 5; //3HG--CG
  AA[p].atom[15].nbonded_neighbors = 1; // OXT
  AA[p].atom[15].bonded_neighbor[0] = 2; // OXT-C


  for ( k = 0; k < 10; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 10; k < 15; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template for building  HA
  AA[p].atom[10].ta[0] = 1; //  CA
  AA[p].atom[10].ta[1] = 0; //  N
  AA[p].atom[10].ta[2] = 2; //  C

  //bk   template for building 2HB
  AA[p].atom[11].ta[0] = 4; //  CB
  AA[p].atom[11].ta[1] = 1; //  CA
  AA[p].atom[11].ta[2] = 5; //  CG

  //bk   template for building 3HB
  AA[p].atom[12].ta[0] = 4; //  CB
  AA[p].atom[12].ta[1] = 1; //  CA
  AA[p].atom[12].ta[2] = 5; //  CG

  //bk   template for building 2HG
  AA[p].atom[13].ta[0] = 5; //  CG
  AA[p].atom[13].ta[1] = 4; //  CB
  AA[p].atom[13].ta[2] = 6; //  CD

  //bk   template for building 3HG
  AA[p].atom[14].ta[0] = 5; //  CG
  AA[p].atom[14].ta[1] = 4; //  CB
  AA[p].atom[14].ta[2] = 6; //  CD


  //bk   chi angles required to build atoms GLU
  //bk   chi angles needed for building  CG
  AA[p].chi_required[0] [5] = true;

  //bk   chi angles needed for building  CD
  AA[p].chi_required[0] [6] = true;
  AA[p].chi_required[1] [6] = true;

  //bk   chi angles needed for building  OE1
  AA[p].chi_required[0] [7] = true;
  AA[p].chi_required[1] [7] = true;
  AA[p].chi_required[2] [7] = true;

  //bk   chi angles needed for building  OE2
  AA[p].chi_required[0] [8] = true;
  AA[p].chi_required[1] [8] = true;
  AA[p].chi_required[2] [8] = true;

  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [11] = true;

  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [12] = true;

  //bk   chi angles needed for building 2HG
  AA[p].chi_required[0] [13] = true;
  AA[p].chi_required[1] [13] = true;

  //bk   chi angles needed for building 3HG
  AA[p].chi_required[0] [14] = true;
  AA[p].chi_required[1] [14] = true;




  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   CG
  //bk   four atoms that define chi angle  2
  AA[p].chi_atoms[1] [0] = 1; //   CA
  AA[p].chi_atoms[1] [1] = 4; //   CB
  AA[p].chi_atoms[1] [2] = 5; //   CG
  AA[p].chi_atoms[1] [3] = 6; //   CD
  //bk   four atoms that define chi angle  3
  AA[p].chi_atoms[2] [0] = 4; //   CB
  AA[p].chi_atoms[2] [1] = 5; //   CG
  AA[p].chi_atoms[2] [2] = 6; //   CD
  AA[p].chi_atoms[2] [3] = 7; //   OE1


 AA[p].atom[0].icoor[0]  = 0.09036078 ;	//   N
  AA[p].atom[0].icoor[1]  = 1.20944    ;	//   N
  AA[p].atom[0].icoor[2]  = -0.09877955;	//   N
  AA[p].atom[1].icoor[0]  = -0.6644503 ;	//   CA
AA[p].atom[1].icoor[1]  = 1.887128   ;	//   CA
AA[p].atom[1].icoor[2]  = -1.167784  ;	//   CA
AA[p].atom[2].icoor[0]  = -0.6415405 ;	//   C
AA[p].atom[2].icoor[1]  = 3.403009   ;	//   C
AA[p].atom[2].icoor[2]  = -0.9416856 ;	//   C
AA[p].atom[3].icoor[0]  =  0.5018099 ;	//   O
AA[p].atom[3].icoor[1]  = 3.910433   ;	//   O
AA[p].atom[3].icoor[2]  = -0.8834405 ;	//   O
AA[p].atom[4].icoor[0]  = -2.078524  ;	//   CB
AA[p].atom[4].icoor[1]  = 1.299661   ;	//   CB
AA[p].atom[4].icoor[2]  = -1.288965  ;	//   CB
AA[p].atom[5].icoor[0]  = -2.807634  ;	//   CG
AA[p].atom[5].icoor[1]  = 1.790975   ;	//   CG
AA[p].atom[5].icoor[2]  = -2.546677  ;	//   CG
AA[p].atom[6].icoor[0]  = -4.223983  ;	//   CD
AA[p].atom[6].icoor[1]  = 1.236637   ;	//   CD
AA[p].atom[6].icoor[2]  = -2.597517  ;	//   CD
AA[p].atom[7].icoor[0]  = -5.044554  ;	//   OE1
AA[p].atom[7].icoor[1]  = 1.733033   ;	//   OE1
AA[p].atom[7].icoor[2]  = -1.797006  ;	//   OE1
AA[p].atom[8].icoor[0]  = -4.450723  ;	//   OE2
AA[p].atom[8].icoor[1]  = 0.3260741  ;	//   OE2
AA[p].atom[8].icoor[2]  = -3.421581  ;	//   OE2
AA[p].atom[9].icoor[0]  = -0.4043283 ;	//   H
AA[p].atom[9].icoor[1]  = 1.301303   ;	//   H
AA[p].atom[9].icoor[2]  = 0.7769806  ;	//   H
AA[p].atom[10].icoor[0] = -0.1401839 ;	//   HA
AA[p].atom[10].icoor[1] = 1.690991   ;	//   HA
AA[p].atom[10].icoor[2] = -2.10393   ;	//   HA
AA[p].atom[11].icoor[0] = -2.006071  ;	//  2HB
AA[p].atom[11].icoor[1] = 0.2127525  ;	//  2HB
AA[p].atom[11].icoor[2] = -1.352877  ;	//  2HB
AA[p].atom[12].icoor[0] = -2.670931  ;	//  3HB
AA[p].atom[12].icoor[1] = 1.54428    ;	//  3HB
AA[p].atom[12].icoor[2] = -0.4052964 ;	//  3HB
AA[p].atom[13].icoor[0] = -2.881274  ;	//  2HG
AA[p].atom[13].icoor[1] = 2.877187   ;	//  2HG
AA[p].atom[13].icoor[2] = -2.559357  ;	//  2HG
AA[p].atom[14].icoor[0] = -2.268831  ;	//  3HG
AA[p].atom[14].icoor[1] = 1.463546   ;	//  3HG
AA[p].atom[14].icoor[2] = -3.436255  ;	//  3HG
AA[p].atom[15].icoor[0] = -1.674247  ;	//  OXT
AA[p].atom[15].icoor[1] = 3.928742   ;	//  OXT
AA[p].atom[15].icoor[2] = -0.470977  ;	//  OXT


}

//FIN DE AMINOACIDO GLUTAMATO C-Terminal



void init_ALA(Convention opt)
{


  //AMINOACIDO ALANINA
  int k, p =ALA;

  strcpy( AA[p].aa_name3, "ALA" );
  AA[p].aa_name1 = 'A';
  AA[p].mass = 71.08;


  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = false;
  AA[p].aa_is_nonpolar = true;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = false;

  AA[p].natoms = 10;
  AA[p].nheavyatoms = 5;
  AA[p].nchi = 0;

  AA[p].chi_atoms = NULL;
  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

  AA[p].chi_types = NULL;
  AA[p].natoms_EEF1 = 5;

  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " H  " ); // 3H of N-term
  strcpy( AA[p].atom[6].atom_name, " HA " );
  strcpy( AA[p].atom[7].atom_name, "1HB " );
  strcpy( AA[p].atom[8].atom_name, "2HB " );
  strcpy( AA[p].atom[9].atom_name, "3HB " );

  strcpy( AA[p].atom[10].atom_name," OXT" );
  strcpy( AA[p].atom[11].atom_name,"2H  " );
  strcpy( AA[p].atom[12].atom_name,"3H  " );

  switch (opt)
  {
  case Rosseta:
  case Sybil:
    AA[p].atom[0].charge = -0.47; //    N  <--- Nt(-0.30)
    AA[p].atom[1].charge = 0.07; //    CA  <--- Nt(0.21)
    AA[p].atom[2].charge = 0.51; //    C  <--- Ct(0.34)
    AA[p].atom[3].charge = -0.51; //    O  <--- Ct(-0.67)
    AA[p].atom[4].charge = -0.27; //    CB
    AA[p].atom[5].charge = 0.31; //    H  <--- Nt(0.33)
    AA[p].atom[6].charge = 0.09; //    HA  <--- Nt(0.10)
    AA[p].atom[7].charge = 0.09; //   1HB
    AA[p].atom[8].charge = 0.09; //   2HB
    AA[p].atom[9].charge = 0.09; //   3HB

    AA[p].atom[10].charge = -0.67; //   0XT
    AA[p].atom[11].charge = 0.33; //  1H
    AA[p].atom[12].charge = 0.33; //  2H
    break;

   case ICM:
      AA[p].atom[0].charge = -0.356; //    N
      AA[p].atom[1].charge = 0.064; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = -0.09; //    CB
      AA[p].atom[5].charge = 0.176; //    H
      AA[p].atom[6].charge = 0.02; //    HA
      AA[p].atom[7].charge = 0.04; //   1HB
      AA[p].atom[8].charge = 0.04; //   2HB
      AA[p].atom[9].charge = 0.04; //   3HB

      AA[p].atom[10].charge = -0.67; //   0XT
      AA[p].atom[11].charge = 0.33; //  1H
      AA[p].atom[12].charge = 0.33; //  2H
      break;
   case EEF1:
      AA[p].atom[0].charge = -0.35; //    N
      AA[p].atom[1].charge =  0.10; //    CA
      AA[p].atom[2].charge =  0.55; //    C
      AA[p].atom[3].charge = -0.55; //    O
      AA[p].atom[4].charge =  0.00; //    CB
      AA[p].atom[5].charge =  0.25; //    H
      AA[p].atom[6].charge =  0.00; //    HA
      AA[p].atom[7].charge =  0.00; //   1HB
      AA[p].atom[8].charge =  0.00; //   2HB
      AA[p].atom[9].charge =  0.00; //   3HB
      break;
   }


  /*   Nt/Ct - Charge Table

   1.0 ALA-N:2H    0.33000 HC
   1.0 ALA-N:HA    0.10000 HB
   1.0 ALA-N:3H    0.33000 HC
   1.0 ALA-N:CB   -0.27000 CT3
   1.0 ALA-N:1HB   0.09000 HA
   1.0 ALA-N:2HB   0.09000 HA
   1.0 ALA-N:3HB   0.09000 HA
   1.0 ALA-N:C     0.51000 C
   1.0 ALA-N:O    -0.51000 O
   1.0 ALA-N:N    -0.30000 NH3  *
   1.0 ALA-N:1H    0.33000 HC   *
   1.0 ALA-N:CA    0.21000 CT1  *

   1.0 ALA-C:HA    0.09000 HB
   1.0 ALA-C:CB   -0.27000 CT3
   1.0 ALA-C:1HB   0.09000 HA
   1.0 ALA-C:2HB   0.09000 HA
   1.0 ALA-C:3HB   0.09000 HA
   1.0 ALA-C:C     0.34000 CC  *
   1.0 ALA-C:O    -0.67000 OC  *
   1.0 ALA-C:OXT  -0.67000 OC
   1.0 ALA-C:N    -0.47000 NH1
   1.0 ALA-C:HN    0.31000 H
   1.0 ALA-C:CA    0.07000 CT1
  */

  switch (opt){
    case EEF1:
	  AA[p].atom[0].fullatom_type =   12; // NH1    N
	  AA[p].atom[1].fullatom_type =    5; // CH1E   CA
	  AA[p].atom[2].fullatom_type =    4; // C      C
	  AA[p].atom[3].fullatom_type =   16; // O      O
	  AA[p].atom[4].fullatom_type =    7; // CH3E   CB
	  AA[p].atom[5].fullatom_type =    0; // H      H
	  AA[p].atom[6].fullatom_type =   -1; //     HA
	  AA[p].atom[7].fullatom_type =   -1; //     1HB
	  AA[p].atom[8].fullatom_type =   -1; //     2HB
	  AA[p].atom[9].fullatom_type =   -1; //     3HB

	  AA[p].atom[10].fullatom_type =  17; //OC      OXT
	  AA[p].atom[11].fullatom_type =   1; //HC      2H
	  AA[p].atom[12].fullatom_type =   1; //HC      3H
	  break;
    case Sybil:
    	  AA[p].atom[0].fullatom_type = 40; // Nbb    N
    	  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
    	  AA[p].atom[2].fullatom_type = 10; // CObb   C
    	  AA[p].atom[3].fullatom_type = 56; // OCbb   O
    	  AA[p].atom[4].fullatom_type = 6; // CH3    CB
    	  AA[p].atom[5].fullatom_type = 66; // HNbb   H
    	  AA[p].atom[6].fullatom_type = 66; // Hapo   HA
    	  AA[p].atom[7].fullatom_type = 66; // Hapo  1HB
    	  AA[p].atom[8].fullatom_type = 66; // Hapo  2HB
    	  AA[p].atom[9].fullatom_type = 66; // Hapo  3HB

    	  AA[p].atom[10].fullatom_type = 57; // OXT
    	  AA[p].atom[11].fullatom_type = 66; // Hpol
    	  AA[p].atom[12].fullatom_type = 66; // Hpol
    	  break;
    default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 5; // CH3    CB
	  AA[p].atom[5].fullatom_type = 25; // HNbb   H
	  AA[p].atom[6].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[7].fullatom_type = 23; // Hapo  1HB
	  AA[p].atom[8].fullatom_type = 23; // Hapo  2HB
	  AA[p].atom[9].fullatom_type = 23; // Hapo  3HB

	  AA[p].atom[10].fullatom_type = 20; // OXT
	  AA[p].atom[11].fullatom_type = 22; // Hpol
	  AA[p].atom[12].fullatom_type = 22; // Hpol
	  break;
  }

  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 5; // N--H (3H) HNbb
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 6; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 7; // CB--1HB
  AA[p].atom[4].bonded_neighbor[2] = 8; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 9; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 1; // H
  AA[p].atom[5].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[6].nbonded_neighbors = 1; // HA
  AA[p].atom[6].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[7].nbonded_neighbors = 1; // 1HB
  AA[p].atom[7].bonded_neighbor[0] = 4; // 1HB--CB
  AA[p].atom[8].nbonded_neighbors = 1; // 2HB
  AA[p].atom[8].bonded_neighbor[0] = 4; // 2HB--CB
  AA[p].atom[9].nbonded_neighbors = 1; // 3HB
  AA[p].atom[9].bonded_neighbor[0] = 4; // 3HB--CB

  AA[p].atom[10].nbonded_neighbors = 1; // OXT
  AA[p].atom[10].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[11].nbonded_neighbors = 1; // 1H
  AA[p].atom[11].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[12].nbonded_neighbors = 1; // 2H
  AA[p].atom[12].bonded_neighbor[0] = 0; // H-N

  for ( k = 0; k < 6; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 6; k < 10; k++ )
    AA[p].atom[k].hastemplate = true;

  //bk   template for building  HA
  AA[p].atom[6].ta[0] = 1; //  CA
  AA[p].atom[6].ta[1] = 0; //  N
  AA[p].atom[6].ta[2] = 2; //  C

  //bk   template for building 1HB
  AA[p].atom[7].ta[0] = 4; //  CB
  AA[p].atom[7].ta[1] = 1; //  CA
  AA[p].atom[7].ta[2] = 0; //  N

  //bk   template for building 2HB
  AA[p].atom[8].ta[0] = 4; //  CB
  AA[p].atom[8].ta[1] = 1; //  CA
  AA[p].atom[8].ta[2] = 0; //  N

  //bk   template for building 3HB
  AA[p].atom[9].ta[0] = 4; //  CB
  AA[p].atom[9].ta[1] = 1; //  CA
  AA[p].atom[9].ta[2] = 0; //  N

  AA[p].atom[0].icoor[0] = 0.0000; //   N
  AA[p].atom[0].icoor[1] = 1.5340; //   N
  AA[p].atom[0].icoor[2] = 3.3340; //   N
  AA[p].atom[1].icoor[0] = 0.0000; //   CA
  AA[p].atom[1].icoor[1] = 2.8360; //   CA
  AA[p].atom[1].icoor[2] = 3.9900; //   CA
  AA[p].atom[2].icoor[0] = 0.0000; //   C
  AA[p].atom[2].icoor[1] = 2.6880; //   C
  AA[p].atom[2].icoor[2] = 5.5060; //   C
  AA[p].atom[3].icoor[0] = 0.0020; //   O
  AA[p].atom[3].icoor[1] = 1.5660; //   O
  AA[p].atom[3].icoor[2] = 6.0130; //   O
  AA[p].atom[4].icoor[0] = -1.1960; //  CB
  AA[p].atom[4].icoor[1] = 3.6590; //   CB
  AA[p].atom[4].icoor[2] = 3.5340; //   CB
  AA[p].atom[5].icoor[0] = -0.0470; //  H
  AA[p].atom[5].icoor[1] = 0.7070; //   H
  AA[p].atom[5].icoor[2] = 3.9120; //   H
  AA[p].atom[6].icoor[0] = 0.9120; //   HA
  AA[p].atom[6].icoor[1] = 3.3660; //   HA
  AA[p].atom[6].icoor[2] = 3.7150; //   HA
  AA[p].atom[7].icoor[0] = -1.1810; //  1HB
  AA[p].atom[7].icoor[1] = 4.6280; //  1HB
  AA[p].atom[7].icoor[2] = 4.0330; //  1HB
  AA[p].atom[8].icoor[0] = -1.1470; //  2HB
  AA[p].atom[8].icoor[1] = 3.8060; //  2HB
  AA[p].atom[8].icoor[2] = 2.4550; //  2HB
  AA[p].atom[9].icoor[0] = -2.1160; //  3HB
  AA[p].atom[9].icoor[1] = 3.1350; //  3HB
  AA[p].atom[9].icoor[2] = 3.7880; //  3HB

  // atom number for backbone HN
  AA[p].HNpos=5;
  // atom number for backbone HA
  AA[p].HApos=6;

  // number of polar hydrogens
  AA[p].nH_polar_complete=1;
  // atom numbers for polar H
  AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
  AA[p].Hpos_polar_complete[0]=5;

  // number of aromatic hydrogens
  AA[p].nH_aromatic_complete=0;

  // number of apolar hydrogens
  AA[p].nH_apolar_complete=4;
  //atom number for apolar hydrogens
  AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
  AA[p].Hpos_apolar_complete[0]=6;
  AA[p].Hpos_apolar_complete[1]=7;
  AA[p].Hpos_apolar_complete[2]=8;
  AA[p].Hpos_apolar_complete[3]=9;

  // number of acceptors
  AA[p].nacceptors=1;
  //acceptor information
  AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
  (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;

  //Number of Hydrogens connection
  AA[p].nH_hydrogen_connexions=5;
  //atoms hydrogens are connected too
  AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
  AA[p].Hbase[0][0]=5; AA[p].Hbase[0][1]=0;
  AA[p].Hbase[1][0]=6; AA[p].Hbase[1][1]=1;
  AA[p].Hbase[2][0]=7; AA[p].Hbase[2][1]=4;
  AA[p].Hbase[3][0]=8; AA[p].Hbase[3][1]=4;
  AA[p].Hbase[4][0]=9; AA[p].Hbase[4][1]=4;

  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=5;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=6;

  AA[p].atom[4].numHydrogens_atm=3;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=7;
  AA[p].atom[4].hydrogens_atm[1]=8;
  AA[p].atom[4].hydrogens_atm[2]=9;


}

//FIN DE AMINOACIDO ALANINA
/**
* Initializes Cysteine
*
* @param opt: Rosseta or ICM
*/
void init_CYS(Convention opt)
{
  //AMINOACIDO CYSTEINA
  int k, p =CYS;

  strcpy( AA[p].aa_name3, "CYS" );
  AA[p].aa_name1 = 'C';
  AA[p].mass = 103.146;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = false;
  AA[p].aa_is_nonpolar = true;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = false;

  AA[p].natoms = 11;
  AA[p].nheavyatoms = 6;
  AA[p].nchi = 1;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );


  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;

  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

  AA[p].chi_types = ( int * ) new int[AA[p].nchi];
  AA[p].chi_types[0] = 2;
  AA[p].natoms_EEF1 = 6;

  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " SG " );
  strcpy( AA[p].atom[6].atom_name, " H  " );
  strcpy( AA[p].atom[7].atom_name, " HA " );
  strcpy( AA[p].atom[8].atom_name, "2HB " );
  strcpy( AA[p].atom[9].atom_name, "3HB " );
  strcpy( AA[p].atom[10].atom_name, " HG " );

  strcpy( AA[p].atom[11].atom_name," OXT" );
  strcpy( AA[p].atom[12].atom_name,"2H  " );
  strcpy( AA[p].atom[13].atom_name,"3H  " );
// OJO COMENTADO MOMENTANEAMENTE
  switch(opt)
  {
  case Rosseta:
  case Sybil:
    AA[p].atom[0].charge = -0.47; //   N  <--- Nt(-0.30)
    AA[p].atom[1].charge = 0.07; //   CA  <--- Nt(0.21)
    AA[p].atom[2].charge = 0.51; //   C  <--- Ct(0.34)
    AA[p].atom[3].charge = -0.51; //   O  <--- Ct(-0.67)
    AA[p].atom[4].charge = -0.11; //   CB
    AA[p].atom[5].charge = -0.23; //   SG
    AA[p].atom[6].charge = 0.31; //   H  <--- Nt(0.33)
    AA[p].atom[7].charge = 0.09; //   HA  <--- Nt(0.10)
    AA[p].atom[8].charge = 0.09; //   2HB
    AA[p].atom[9].charge = 0.09; //   3HB
    AA[p].atom[10].charge = 0.16; //   HG

    AA[p].atom[11].charge = -0.67; //   0XT
    AA[p].atom[12].charge = 0.33; //  1H
    AA[p].atom[13].charge = 0.33; //  2H
    break;
  case ICM:
    AA[p].atom[0].charge = -0.356; //   N
    AA[p].atom[1].charge = 0.064; //   CA
    AA[p].atom[2].charge = 0.45; //   C
    AA[p].atom[3].charge = -0.384; //   O
    AA[p].atom[4].charge = -0.105; //   CB
    AA[p].atom[5].charge = 0.015; //   SG
    AA[p].atom[6].charge = 0.176; //   H
    AA[p].atom[7].charge = 0.02; //   HA
    AA[p].atom[8].charge = 0.055; //   2HB
    AA[p].atom[9].charge = 0.055; //   3HB
    AA[p].atom[10].charge = 0.01; //   HG

    AA[p].atom[11].charge = -0.67; //   0XT
    AA[p].atom[12].charge = 0.33; //  1H
    AA[p].atom[13].charge = 0.33; //  2H
    break;
  case EEF1:
    AA[p].atom[0].charge = -0.35; //   N
    AA[p].atom[1].charge =  0.10; //   CA
    AA[p].atom[2].charge =  0.55; //   C
    AA[p].atom[3].charge = -0.55; //   O
    AA[p].atom[4].charge =  0.19; //   CB
    AA[p].atom[5].charge = -0.19; //   SG
    AA[p].atom[6].charge =  0.25; //   H
    AA[p].atom[7].charge =  0.00; //   HA
    AA[p].atom[8].charge =  0.00; //   2HB
    AA[p].atom[9].charge =  0.00; //   3HB
    AA[p].atom[10].charge = 0.00; //   HG
    break;

  }

//CARGAS INCORRECTAS
/*AA[p].atom[0].charge = 0.5; //   N
AA[p].atom[1].charge = 0.5; //   CA
AA[p].atom[2].charge = 0.5; //   C
AA[p].atom[3].charge = 0.5; //   O
AA[p].atom[4].charge = 0.5; //   CB
AA[p].atom[5].charge = 0.5; //   SG
AA[p].atom[6].charge = 0.5; //   H
AA[p].atom[7].charge = 0.5; //   HA
AA[p].atom[8].charge = 0.5; //   2HB
AA[p].atom[9].charge = 0.5; //   3HB
AA[p].atom[10].charge = 0.5; //   HG

AA[p].atom[11].charge = 0.5; //   0XT
AA[p].atom[12].charge = 0.5; //  1H
AA[p].atom[13].charge = 0.5; //  2H
*/

  /*   Nt/Ct - Charge Table

   1.0 CYS-N:2H    0.33000 HC
   1.0 CYS-N:HA    0.10000 HB
   1.0 CYS-N:3H    0.33000 HC
   1.0 CYS-N:CB   -0.11000 CT2
   1.0 CYS-N:HB1   0.09000 HA
   1.0 CYS-N:1HB   0.09000 HA
   1.0 CYS-N:SG   -0.23000 S
   1.0 CYS-N:HG1   0.16000 HS
   1.0 CYS-N:C     0.51000 C
   1.0 CYS-N:O    -0.51000 O
   1.0 CYS-N:N    -0.30000 NH3
   1.0 CYS-N:1H    0.33000 HC
   1.0 CYS-N:CA    0.21000 CT1

   1.0 CYS-C:HA    0.09000 HB
   1.0 CYS-C:CB   -0.11000 CT2
   1.0 CYS-C:HB1   0.09000 HA
   1.0 CYS-C:1HB   0.09000 HA
   1.0 CYS-C:SG   -0.23000 S
   1.0 CYS-C:HG1   0.16000 HS
   1.0 CYS-C:C     0.34000 CC
   1.0 CYS-C:O    -0.67000 OC
   1.0 CYS-C:OXT  -0.67000 OC
   1.0 CYS-C:N    -0.47000 NH1
   1.0 CYS-C:HN    0.31000 H
   1.0 CYS-C:CA    0.07000 CT1
*/

  switch(opt) {
  	case EEF1:
	  AA[p].atom[0].fullatom_type =   12; // NH1    N
	  AA[p].atom[1].fullatom_type =    5; // CH1E   CA
	  AA[p].atom[2].fullatom_type =    4; // C      C
	  AA[p].atom[3].fullatom_type =   16; // O      O
	  AA[p].atom[4].fullatom_type =    6; // CH2E   CB
	  AA[p].atom[5].fullatom_type =   21; // SH1E   SG
	  AA[p].atom[6].fullatom_type =    0; // H      H
	  AA[p].atom[7].fullatom_type =   -1; //      HA
	  AA[p].atom[8].fullatom_type =   -1; //      2HB
	  AA[p].atom[9].fullatom_type =   -1; //      3HB
	  AA[p].atom[10].fullatom_type =  -1; //      HG

	  AA[p].atom[11].fullatom_type =  17; //OC      OXT
	  AA[p].atom[12].fullatom_type =   1; //HC      Hpol
	  AA[p].atom[13].fullatom_type =   1; //HC      Hpol
	  break;
  	case Sybil:
  		  AA[p].atom[0].fullatom_type = 40; // Nbb    N
  		  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
  		  AA[p].atom[2].fullatom_type = 10; // CObb   C
  		  AA[p].atom[3].fullatom_type = 56; // OCbb   O
  		  AA[p].atom[4].fullatom_type = 8; // CH2    CB
  		  AA[p].atom[5].fullatom_type = 63; // S      SG
  		  AA[p].atom[6].fullatom_type = 66; // HNbb   H
  		  AA[p].atom[7].fullatom_type = 66; // Hapo   HA
  		  AA[p].atom[8].fullatom_type = 66; // Hapo   2HB
  		  AA[p].atom[9].fullatom_type = 66; // Hapo   3HB
  		  AA[p].atom[10].fullatom_type = 66; // Hpol   HG

  		  AA[p].atom[11].fullatom_type = 57; // OXT
  		  AA[p].atom[12].fullatom_type = 66; // Hpol
  		  AA[p].atom[13].fullatom_type = 66; // Hpol
  		  break;
      default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 16; // S      SG
	  AA[p].atom[6].fullatom_type = 25; // HNbb   H
	  AA[p].atom[7].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[8].fullatom_type = 23; // Hapo   2HB
	  AA[p].atom[9].fullatom_type = 23; // Hapo   3HB
	  AA[p].atom[10].fullatom_type = 22; // Hpol   HG

	  AA[p].atom[11].fullatom_type = 20; // OXT
	  AA[p].atom[12].fullatom_type = 22; // Hpol
	  AA[p].atom[13].fullatom_type = 22; // Hpol
	  break;
  }
  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 6; // N--H (3H) HNbb
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 7; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--SG
  AA[p].atom[4].bonded_neighbor[2] = 8; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 9; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 2; // SG
  AA[p].atom[5].bonded_neighbor[0] = 4; // SG--CB
  AA[p].atom[5].bonded_neighbor[1] = 10; // SG--HG1
  AA[p].atom[6].nbonded_neighbors = 1; // H
  AA[p].atom[6].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[7].nbonded_neighbors = 1; // HA
  AA[p].atom[7].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[8].nbonded_neighbors = 1; //2HB
  AA[p].atom[8].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[9].nbonded_neighbors = 1; //3HB
  AA[p].atom[9].bonded_neighbor[0] = 4; //3HB--CB
  AA[p].atom[10].nbonded_neighbors = 1; // HG
  AA[p].atom[10].bonded_neighbor[0] = 5; // HG--SG

  AA[p].atom[11].nbonded_neighbors = 1; // OXT
  AA[p].atom[11].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[12].nbonded_neighbors = 1; // 1H
  AA[p].atom[12].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[13].nbonded_neighbors = 1; // 2H
  AA[p].atom[13].bonded_neighbor[0] = 0; // H-N



  for ( k = 0; k < 7; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 7; k < 11; k++ )
    AA[p].atom[k].hastemplate = true;

  //bk   template for building  HA
  AA[p].atom[7].ta[0] = 1; //  CA
  AA[p].atom[7].ta[1] = 0; //  N
  AA[p].atom[7].ta[2] = 2; //  C

  //bk   template for building 2HB
  AA[p].atom[8].ta[0] = 4; //  CB
  AA[p].atom[8].ta[1] = 1; //  CA
  AA[p].atom[8].ta[2] = 5; //  SG

  //bk   template for building 3HB
  AA[p].atom[9].ta[0] = 4; //  CB
  AA[p].atom[9].ta[1] = 1; //  CA
  AA[p].atom[9].ta[2] = 5; //  SG

  //bk   template for building HG
  AA[p].atom[10].ta[0] = 5; //  SG
  AA[p].atom[10].ta[1] = 4; //  CB
  AA[p].atom[10].ta[2] = 1; //  CA


  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   SG

  //bk   chi angles required to build atoms CYS
  //bk   chi angles needed for building  SG
  AA[p].chi_required[0] [5] = true;

  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [8] = true;

  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [9] = true;

  //bk   chi angles needed for building  HG
  AA[p].chi_required[0] [10] = true;



  //bk   template coordinates for the amino acid
  AA[p].atom[0].icoor[0] = 0.0000; //   N
  AA[p].atom[0].icoor[1] = 3.8070; //   N
  AA[p].atom[0].icoor[2] = 6.2240; //   N
  AA[p].atom[1].icoor[0] = 0.0000; //   CA
  AA[p].atom[1].icoor[1] = 3.7830; //   CA
  AA[p].atom[1].icoor[2] = 7.6880; //   CA
  AA[p].atom[2].icoor[0] = 0.0000; //   C
  AA[p].atom[2].icoor[1] = 5.1680; //   C
  AA[p].atom[2].icoor[2] = 8.3290; //   C
  AA[p].atom[3].icoor[0] = 0.0050; //   O
  AA[p].atom[3].icoor[1] = 6.1860; //   O
  AA[p].atom[3].icoor[2] = 7.6360; //   O
  AA[p].atom[4].icoor[0] = 1.3120; //   CB
  AA[p].atom[4].icoor[1] = 3.0690; //   CB
  AA[p].atom[4].icoor[2] = 8.0140; //   CB
  AA[p].atom[5].icoor[0] = 1.6350; //   SG
  AA[p].atom[5].icoor[1] = 2.8650; //   SG
  AA[p].atom[5].icoor[2] = 9.7820; //   SG
  AA[p].atom[6].icoor[0] = -0.0250; //   H
  AA[p].atom[6].icoor[1] = 4.7040; //   H
  AA[p].atom[6].icoor[2] = 5.7600; //   H
  AA[p].atom[7].icoor[0] = -0.8240; //   HA
  AA[p].atom[7].icoor[1] = 3.2090; //   HA
  AA[p].atom[7].icoor[2] = 8.1120; //   HA
  AA[p].atom[8].icoor[0] = 1.3110; //   2HB
  AA[p].atom[8].icoor[1] = 2.0650; //   2HB
  AA[p].atom[8].icoor[2] = 7.5890; //   2HB
  AA[p].atom[9].icoor[0] = 2.1560; //   3HB
  AA[p].atom[9].icoor[1] = 3.6320; //   3HB
  AA[p].atom[9].icoor[2] = 7.6160; //   3HB
  AA[p].atom[10].icoor[0] = 2.7990; //   HG
  AA[p].atom[10].icoor[1] = 2.2390; //   HG
  AA[p].atom[10].icoor[2] = 9.6390; //   HG

  // atom number for backbone HN
  AA[p].HNpos=6;
  // atom number for backbone HA
  AA[p].HApos=7;

  // number of polar hydrogens
  AA[p].nH_polar_complete=2;
  // atom numbers for polar H
  AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
  AA[p].Hpos_polar_complete[0]=6;
  AA[p].Hpos_polar_complete[1]=10;

  // number of aromatic hydrogens
  AA[p].nH_aromatic_complete=0;

  // number of apolar hydrogens
  AA[p].nH_apolar_complete=3;
  //atom number for apolar hydrogens
  AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
  AA[p].Hpos_apolar_complete[0]=7;
  AA[p].Hpos_apolar_complete[1]=8;
  AA[p].Hpos_apolar_complete[2]=9;

  // number of acceptors
  AA[p].nacceptors=1;
  //acceptor information
  AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
  (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;

  //Number of Hydrogens connection
  AA[p].nH_hydrogen_connexions=5;
  //atoms hydrogens are connected too
  AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
  AA[p].Hbase[0][0]=6; AA[p].Hbase[0][1]=0;
  AA[p].Hbase[1][0]=7; AA[p].Hbase[1][1]=1;
  AA[p].Hbase[2][0]=8; AA[p].Hbase[2][1]=4;
  AA[p].Hbase[3][0]=9; AA[p].Hbase[3][1]=4;
  AA[p].Hbase[4][0]=10; AA[p].Hbase[4][1]=5;

  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=6;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=7;

  AA[p].atom[4].numHydrogens_atm=2;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=8;
  AA[p].atom[4].hydrogens_atm[1]=9;

  AA[p].atom[5].numHydrogens_atm=1;
  AA[p].atom[5].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[5].numHydrogens_atm);
  AA[p].atom[5].hydrogens_atm[0]=10;

}

//FIN DE AMINOACIDO CYSTEINA
/**
* Initializes Aspartate
*
* @param opt: Rosseta or ICM
*/
void init_ASP(Convention opt)
{
  //AMINOACIDO ASPARTATO
  int k, p = ASP;

  strcpy( AA[p].aa_name3, "ASP" );
  AA[p].aa_name1 = 'D';
  AA[p].mass = 114.083;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = true;
  AA[p].aa_is_nonpolar = false;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = true;

  AA[p].natoms = 12;
  AA[p].nheavyatoms = 8;
  AA[p].nchi = 2;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );

  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;
  /* for ( int chino = 0; chino < AA[p].nchi; chino++ ) for ( int i = 0 ; i < AA[p].natoms; i++)
  if (AA[p].chi_required[chino][i]) printf( "--->%d %d true\n", chino, i); else  printf("--->%d %d false\n", chino, i); */

  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms +3) * sizeof( t_aa_atom_p ) );
  AA[p].chi_types = ( int * ) new int[AA[p].nchi];
  AA[p].chi_types[0] = 2;
  AA[p].chi_types[1] = 4;
  AA[p].natoms_EEF1 = 8;


  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " CG " );
  strcpy( AA[p].atom[6].atom_name, " OD1" );
  strcpy( AA[p].atom[7].atom_name, " OD2" );
  strcpy( AA[p].atom[8].atom_name, " H  " );
  strcpy( AA[p].atom[9].atom_name, " HA " );
  strcpy( AA[p].atom[10].atom_name, "2HB " );
  strcpy( AA[p].atom[11].atom_name, "3HB " );

  strcpy( AA[p].atom[12].atom_name," OXT" );
  strcpy( AA[p].atom[13].atom_name,"2H  " );
  strcpy( AA[p].atom[14].atom_name,"3H  " );


  switch(opt)
  {
    case Rosseta:
    case Sybil:
      AA[p].atom[0].charge = -0.47; //    N  <--- Nt(-0.30)
      AA[p].atom[1].charge = 0.07; //    CA  <--- Nt(0.21)
      AA[p].atom[2].charge = 0.51; //    C  <--- Ct(0.34)
      AA[p].atom[3].charge = -0.51; //    O  <--- Ct(-0.67)
      AA[p].atom[4].charge = -0.28; //    CB
      AA[p].atom[5].charge = 0.62; //    CG
      AA[p].atom[6].charge = -0.76; //    OD1
      AA[p].atom[7].charge = -0.76; //    OD2
      AA[p].atom[8].charge = 0.31; //    H  <--- Nt(0.33)
      AA[p].atom[9].charge = 0.09; //    HA  <--- Nt(0.10)
      AA[p].atom[10].charge = 0.09; //    2HB
      AA[p].atom[11].charge = 0.09; //    3HB

      AA[p].atom[12].charge = -0.67; //   0XT
      AA[p].atom[13].charge = 0.33; //  1H
      AA[p].atom[14].charge = 0.33; //  2H
      break;

    case ICM:
      AA[p].atom[0].charge = -0.356; //    N
      AA[p].atom[1].charge = 0.06; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = -0.17; //    CB
      AA[p].atom[5].charge = 0.5; //    CG
      AA[p].atom[6].charge = -0.57; //    OD1
      AA[p].atom[7].charge = -0.57; //    OD2
      AA[p].atom[8].charge = 0.176; //    H
      AA[p].atom[9].charge = 0.024; //    HA
      AA[p].atom[10].charge = 0.02; //    2HB
      AA[p].atom[11].charge = 0.02; //    3HB

      AA[p].atom[12].charge = -0.67; //   0XT
      AA[p].atom[13].charge = 0.33; //  1H
      AA[p].atom[14].charge = 0.33; //  2H
      break;
    case EEF1:
      AA[p].atom[0].charge = -0.35; //    N
      AA[p].atom[1].charge =  0.10; //    CA
      AA[p].atom[2].charge =  0.55; //    C
      AA[p].atom[3].charge = -0.55; //    O
      AA[p].atom[4].charge = -0.15; //    CB
      AA[p].atom[5].charge =  1.35; //    CG
      AA[p].atom[6].charge = -0.60; //    OD1
      AA[p].atom[7].charge = -0.60; //    OD2
      AA[p].atom[8].charge =  0.25; //    H
      AA[p].atom[9].charge =  0.00; //    HA
      AA[p].atom[10].charge = 0.00; //    2HB
      AA[p].atom[11].charge = 0.00; //    3HB
      break;

 }

  /*   Nt/Ct - Charge Table

   1.0 ASP-N:2H    0.33000 HC
   1.0 ASP-N:HA    0.10000 HB
   1.0 ASP-N:3H    0.33000 HC
   1.0 ASP-N:CB   -0.28000 CT2
   1.0 ASP-N:HB1   0.09000 HA
   1.0 ASP-N:1HB   0.09000 HA
   1.0 ASP-N:CG    0.62000 CC
   1.0 ASP-N:OD1  -0.76000 OC
   1.0 ASP-N:OD2  -0.76000 OC
   1.0 ASP-N:C     0.51000 C
   1.0 ASP-N:O    -0.51000 O
   1.0 ASP-N:N    -0.30000 NH3
   1.0 ASP-N:1H    0.33000 HC
   1.0 ASP-N:CA    0.21000 CT1

   1.0 ASP-C:HA    0.09000 HB
   1.0 ASP-C:CB   -0.28000 CT2
   1.0 ASP-C:HB1   0.09000 HA
   1.0 ASP-C:1HB   0.09000 HA
   1.0 ASP-C:CG    0.62000 CC
   1.0 ASP-C:OD1  -0.76000 OC
   1.0 ASP-C:OD2  -0.76000 OC
   1.0 ASP-C:C     0.34000 CC
   1.0 ASP-C:O    -0.67000 OC
   1.0 ASP-C:OXT  -0.67000 OC
   1.0 ASP-C:N    -0.47000 NH1
   1.0 ASP-C:HN    0.31000 H
   1.0 ASP-C:CA    0.07000 CT1
  */

  switch(opt){
    case EEF1:
	  AA[p].atom[0].fullatom_type =   12; // NH1    N
	  AA[p].atom[1].fullatom_type =    5; // CH1E   CA
	  AA[p].atom[2].fullatom_type =    4; // C      C
	  AA[p].atom[3].fullatom_type =   16; // O      O
	  AA[p].atom[4].fullatom_type =    6; // CH2E   CB
	  AA[p].atom[5].fullatom_type =    4; // C      CG
	  AA[p].atom[6].fullatom_type =   17; // OC     OD1
	  AA[p].atom[7].fullatom_type =   17; // OC     OD2
	  AA[p].atom[8].fullatom_type =    0; // H      H
	  AA[p].atom[9].fullatom_type =   -1; //      HA
	  AA[p].atom[10].fullatom_type =  -1; //      2HB
	  AA[p].atom[11].fullatom_type =  -1; //      3HB

	  AA[p].atom[12].fullatom_type =  17; //OC      OXT
	  AA[p].atom[13].fullatom_type =   1; //HC	   2H
	  AA[p].atom[14].fullatom_type =   1; //HC	   3H
	  break;
    case Sybil:
      AA[p].atom[0].fullatom_type = 40; // Nbb    N
   	  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
   	  AA[p].atom[2].fullatom_type = 10; // CObb   C
   	  AA[p].atom[3].fullatom_type = 56; // OCbb   O
   	  AA[p].atom[4].fullatom_type = 8; // CH2    CB
   	  AA[p].atom[5].fullatom_type = 15; // COO    CG
   	  AA[p].atom[6].fullatom_type = 57; // OOC    OD1
   	  AA[p].atom[7].fullatom_type = 57; // OOC    OD2
   	  AA[p].atom[8].fullatom_type = 66; // HNbb   H
   	  AA[p].atom[9].fullatom_type = 66; // Hapo   HA
   	  AA[p].atom[10].fullatom_type = 66; // Hapo  2HB
   	  AA[p].atom[11].fullatom_type = 66; // Hapo  3HB

   	  AA[p].atom[12].fullatom_type = 57; // OXT
   	  AA[p].atom[13].fullatom_type = 66; // Hpol
   	  AA[p].atom[14].fullatom_type = 66; // Hpol
   	  break;
    default:
   	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 2; // COO    CG
	  AA[p].atom[6].fullatom_type = 15; // OOC    OD1
	  AA[p].atom[7].fullatom_type = 15; // OOC    OD2
	  AA[p].atom[8].fullatom_type = 25; // HNbb   H
	  AA[p].atom[9].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[10].fullatom_type = 23; // Hapo  2HB
	  AA[p].atom[11].fullatom_type = 23; // Hapo  3HB

	  AA[p].atom[12].fullatom_type = 20; // OXT
	  AA[p].atom[13].fullatom_type = 22; // Hpol
	  AA[p].atom[14].fullatom_type = 22; // Hpol
	  break;
  }


  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 8; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 9; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--CG
  AA[p].atom[4].bonded_neighbor[2] = 10; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 11; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 3; // CG
  AA[p].atom[5].bonded_neighbor[0] = 4; // CG--CB
  AA[p].atom[5].bonded_neighbor[1] = 6; // CG--OD1
  AA[p].atom[5].bonded_neighbor[2] = 7; // CG--OD2
  AA[p].atom[6].nbonded_neighbors = 1; // OD1
  AA[p].atom[6].bonded_neighbor[0] = 5; // OD1--CG
  AA[p].atom[7].nbonded_neighbors = 1; // OD2
  AA[p].atom[7].bonded_neighbor[0] = 5; // OD2--CG
  AA[p].atom[8].nbonded_neighbors = 1; // H
  AA[p].atom[8].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[9].nbonded_neighbors = 1; // HA
  AA[p].atom[9].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[10].nbonded_neighbors = 1; //2HB
  AA[p].atom[10].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[11].nbonded_neighbors = 1; //3HB
  AA[p].atom[11].bonded_neighbor[0] = 4; //3HB--CB

  AA[p].atom[12].nbonded_neighbors = 1; // OXT
  AA[p].atom[12].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[13].nbonded_neighbors = 1; // 1H
  AA[p].atom[13].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[14].nbonded_neighbors = 1; // 2H
  AA[p].atom[14].bonded_neighbor[0] = 0; // H-N


  for ( k = 0; k < 9; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 9; k < 12; k++ )
    AA[p].atom[k].hastemplate = true;

  //bk   template for building  HA
  AA[p].atom[9].ta[0] = 1; //  CA
  AA[p].atom[9].ta[1] = 0; //  N
  AA[p].atom[9].ta[2] = 2; //  C

  //bk   template for building 2HB
  AA[p].atom[10].ta[0] = 4; //  CB
  AA[p].atom[10].ta[1] = 1; //  CA
  AA[p].atom[10].ta[2] = 5; //  CG

  //bk   template for building 3HB
  AA[p].atom[11].ta[0] = 4; //  CB
  AA[p].atom[11].ta[1] = 1; //  CA
  AA[p].atom[11].ta[2] = 5; //  CG



  //bk   chi angles required to build atoms ASP
  //bk   chi angles needed for building  CG
  AA[p].chi_required[0] [5] = true;

  //bk   chi angles needed for building  OD1
  AA[p].chi_required[0] [6] = true;
  AA[p].chi_required[1] [6] = true;

  //bk   chi angles needed for building  OD2
  AA[p].chi_required[0] [7] = true;
  AA[p].chi_required[1] [7] = true;

  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [10] = true;

  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [11] = true;



  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   CG

  //bk   four atoms that define chi angle  2
  AA[p].chi_atoms[1] [0] = 1; //   CA
  AA[p].chi_atoms[1] [1] = 4; //   CB
  AA[p].chi_atoms[1] [2] = 5; //   CG
  AA[p].chi_atoms[1] [3] = 6; //   OD1


  AA[p].atom[0].icoor[0] = 0.0000; //   N
  AA[p].atom[0].icoor[1] = 5.2000; //   N
  AA[p].atom[0].icoor[2] = 9.6580; //   N
  AA[p].atom[1].icoor[0] = 0.0000; //   CA
  AA[p].atom[1].icoor[1] = 6.4570; //   CA
  AA[p].atom[1].icoor[2] = 10.3940; //   CA
  AA[p].atom[2].icoor[0] = 0.0000; //   C
  AA[p].atom[2].icoor[1] = 6.2160; //   C
  AA[p].atom[2].icoor[2] = 11.8980; //   C
  AA[p].atom[3].icoor[0] = 0.0010; //   O
  AA[p].atom[3].icoor[1] = 5.0720; //   O
  AA[p].atom[3].icoor[2] = 12.3530; //   O
  AA[p].atom[4].icoor[0] = -1.2080; //   CB
  AA[p].atom[4].icoor[1] = 7.3100; //   CB
  AA[p].atom[4].icoor[2] = 9.9990; //   CB
  AA[p].atom[5].icoor[0] = -1.2950; //   CG
  AA[p].atom[5].icoor[1] = 8.6490; //   CG
  AA[p].atom[5].icoor[2] = 10.7190; //   CG
  AA[p].atom[6].icoor[0] = -2.2210; //   OD1
  AA[p].atom[6].icoor[1] = 9.3810; //   OD1
  AA[p].atom[6].icoor[2] = 10.4610; //   OD1
  AA[p].atom[7].icoor[0] = -0.5300; //   OD2
  AA[p].atom[7].icoor[1] = 8.8620; //   OD2
  AA[p].atom[7].icoor[2] = 11.6290; //   OD2
  AA[p].atom[8].icoor[0] = -0.0400; //   H
  AA[p].atom[8].icoor[1] = 4.3370; //   H
  AA[p].atom[8].icoor[2] = 10.1830; //   H
  AA[p].atom[9].icoor[0] = 0.9090; //   HA
  AA[p].atom[9].icoor[1] = 7.0160; //   HA
  AA[p].atom[9].icoor[2] = 10.1700; //   HA
  AA[p].atom[10].icoor[0] = -1.2940; //   2HB
  AA[p].atom[10].icoor[1] = 7.4700; //   2HB
  AA[p].atom[10].icoor[2] = 8.9240; //   2HB
  AA[p].atom[11].icoor[0] = -2.0180; //   3HB
  AA[p].atom[11].icoor[1] = 6.6630; //   3HB
  AA[p].atom[11].icoor[2] = 10.3370; //   3HB

  // atom number for backbone HN
  AA[p].HNpos=8;
  // atom number for backbone HA
  AA[p].HApos=9;

  // number of polar hydrogens
  AA[p].nH_polar_complete=1;
  // atom numbers for polar H
  AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
  AA[p].Hpos_polar_complete[0]=8;

  // number of aromatic hydrogens
  AA[p].nH_aromatic_complete=0;

  // number of apolar hydrogens
  AA[p].nH_apolar_complete=3;
  //atom number for apolar hydrogens
  AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
  AA[p].Hpos_apolar_complete[0]=9;
  AA[p].Hpos_apolar_complete[1]=10;
  AA[p].Hpos_apolar_complete[2]=11;

  // number of acceptors
  AA[p].nacceptors=3;
  //acceptor information
  AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
  (AA[p].aBase[1])[0]=6;(AA[p].aBase[1])[1]=5;
  (AA[p].aBase[2])[0]=7;(AA[p].aBase[2])[1]=5;
  (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;
  (AA[p].accpt_pos[1])[0]=1;(AA[p].accpt_pos[1])[1]=6;
  (AA[p].accpt_pos[2])[0]=2;(AA[p].accpt_pos[2])[1]=7;

  //Number of Hydrogens connection
  AA[p].nH_hydrogen_connexions=4;
  //atoms hydrogens are connected too
  AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
  AA[p].Hbase[0][0]=8; AA[p].Hbase[0][1]=0;
  AA[p].Hbase[1][0]=9; AA[p].Hbase[1][1]=1;
  AA[p].Hbase[2][0]=10; AA[p].Hbase[2][1]=4;
  AA[p].Hbase[3][0]=11; AA[p].Hbase[3][1]=4;


  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=8;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=9;

  AA[p].atom[4].numHydrogens_atm=2;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=10;
  AA[p].atom[4].hydrogens_atm[0]=11;

}

//FIN DE AMINOACIDO ASPARTATO

/**
* Initializes Glutamate
*
* @param opt: Rosseta or ICM
*/
void init_GLU(Convention opt)
{

  //AMINOACIDO GLUTAMATO
  int k, p = GLU;


  strcpy( AA[p].aa_name3, "GLU" );
  AA[p].aa_name1 = 'E';
  AA[p].mass = 128.11;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = true;
  AA[p].aa_is_nonpolar = false;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = true;

  AA[p].natoms = 15;
  AA[p].nheavyatoms = 9;
  AA[p].nchi = 3;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );

  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;

  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

// Bug 23/11/2009
//  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );
  AA[p].chi_types = ( int * ) new int[AA[p].nchi];
  AA[p].chi_types[0] = 2;
  AA[p].chi_types[1] = 2;
  AA[p].chi_types[2] = 4;
  AA[p].natoms_EEF1 = 9;


  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " CG " );
  strcpy( AA[p].atom[6].atom_name, " CD " );
  strcpy( AA[p].atom[7].atom_name, " OE1" );
  strcpy( AA[p].atom[8].atom_name, " OE2" );
  strcpy( AA[p].atom[9].atom_name, " H  " );
  strcpy( AA[p].atom[10].atom_name, " HA " );
  strcpy( AA[p].atom[11].atom_name, "2HB " );
  strcpy( AA[p].atom[12].atom_name, "3HB " );
  strcpy( AA[p].atom[13].atom_name, "2HG " );
  strcpy( AA[p].atom[14].atom_name, "3HG " );

  strcpy( AA[p].atom[15].atom_name," OXT" );
  strcpy( AA[p].atom[16].atom_name,"2H  " );
  strcpy( AA[p].atom[17].atom_name,"3H  " );


  switch (opt)
  {
    case Rosseta:
    case Sybil:
      AA[p].atom[0].charge = -0.47; //    N
      AA[p].atom[1].charge = 0.07; //    CA
      AA[p].atom[2].charge = 0.51; //    C
      AA[p].atom[3].charge = -0.51; //    O
      AA[p].atom[4].charge = -0.18; //    CB
      AA[p].atom[5].charge = -0.28; //    CG
      AA[p].atom[6].charge = 0.62; //    CD
      AA[p].atom[7].charge = -0.76; //    OE1
      AA[p].atom[8].charge = -0.76; //    OE2
      AA[p].atom[9].charge = 0.31; //    H
      AA[p].atom[10].charge = 0.09; //    HA
      AA[p].atom[11].charge = 0.09; //    2HB
      AA[p].atom[12].charge = 0.09; //    3HB
      AA[p].atom[13].charge = 0.09; //    2HG
      AA[p].atom[14].charge = 0.09; //    3HG

      AA[p].atom[15].charge = -0.67; //   0XT
      AA[p].atom[16].charge = 0.33; //  1H
      AA[p].atom[17].charge = 0.33; //  2H
      break;
    case ICM:
      AA[p].atom[0].charge = -0.356; //    N
      AA[p].atom[1].charge = 0.064; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = -0.12; //    CB
      AA[p].atom[5].charge = -0.17; //    CG
      AA[p].atom[6].charge = 0.5; //    CD
      AA[p].atom[7].charge = -0.57; //    OE1
      AA[p].atom[8].charge = -0.57; //    OE2
      AA[p].atom[9].charge = 0.176; //    H
      AA[p].atom[10].charge = 0.02; //    HA
      AA[p].atom[11].charge = 0.02; //    2HB
      AA[p].atom[12].charge = 0.02; //    3HB
      AA[p].atom[13].charge = -0.04; //    2HG
      AA[p].atom[14].charge = -0.04; //    3HG

      AA[p].atom[15].charge = -0.67; //   0XT
      AA[p].atom[16].charge = 0.33; //  1H
      AA[p].atom[17].charge = 0.33; //  2H
      break;
   case EEF1:
      AA[p].atom[0].charge = -0.35; //    N
      AA[p].atom[1].charge =  0.10; //    CA
      AA[p].atom[2].charge =  0.55; //    C
      AA[p].atom[3].charge = -0.55; //    O
      AA[p].atom[4].charge =  0.00; //    CB
      AA[p].atom[5].charge = -0.15; //    CG
      AA[p].atom[6].charge =  1.35; //    CD
      AA[p].atom[7].charge = -0.60; //    OE1
      AA[p].atom[8].charge = -0.60; //    OE2
      AA[p].atom[9].charge =  0.25; //    H
      AA[p].atom[10].charge = 0.00; //    HA
      AA[p].atom[11].charge = 0.00; //    2HB
      AA[p].atom[12].charge = 0.00; //    3HB
      AA[p].atom[13].charge = 0.00; //    2HG
      AA[p].atom[14].charge = 0.00; //    3HG
      break;

  }
/*
   1.0 GLU-N:2H    0.33000 HC
   1.0 GLU-N:HA    0.10000 HB
   1.0 GLU-N:3H    0.33000 HC
   1.0 GLU-N:CB   -0.18000 CT2
   1.0 GLU-N:HB1   0.09000 HA
   1.0 GLU-N:1HB   0.09000 HA
   1.0 GLU-N:CG   -0.28000 CT2
   1.0 GLU-N:HG1   0.09000 HA
   1.0 GLU-N:1HG   0.09000 HA
   1.0 GLU-N:CD    0.62000 CC
   1.0 GLU-N:OE1  -0.76000 OC
   1.0 GLU-N:OE2  -0.76000 OC
   1.0 GLU-N:C     0.51000 C
   1.0 GLU-N:N    -0.30000 NH3
   1.0 GLU-N:O    -0.51000 O
   1.0 GLU-N:1H    0.33000 HC
   1.0 GLU-N:CA    0.21000 CT1

   1.0 GLU-C:HA    0.09000 HB
   1.0 GLU-C:CB   -0.18000 CT2
   1.0 GLU-C:HB1   0.09000 HA
   1.0 GLU-C:1HB   0.09000 HA
   1.0 GLU-C:CG   -0.28000 CT2
   1.0 GLU-C:HG1   0.09000 HA
   1.0 GLU-C:1HG   0.09000 HA
   1.0 GLU-C:CD    0.62000 CC
   1.0 GLU-C:O    -0.67000 OC
   1.0 GLU-C:OE1  -0.76000 OC
   1.0 GLU-C:OXT  -0.67000 OC
   1.0 GLU-C:OE2  -0.76000 OC
   1.0 GLU-C:C     0.34000 CC
   1.0 GLU-C:N    -0.47000 NH1
   1.0 GLU-C:HN    0.31000 H
   1.0 GLU-C:CA    0.07000 CT1
*/

  switch ( opt ) {
  	case EEF1:
	  AA[p].atom[0].fullatom_type =  12; // NH1    N
	  AA[p].atom[1].fullatom_type =   5; // CH1E   CA
	  AA[p].atom[2].fullatom_type =   4; // C      C
	  AA[p].atom[3].fullatom_type =  16; // O      O
	  AA[p].atom[4].fullatom_type =   6; // CH2E   CB
	  AA[p].atom[5].fullatom_type =   6; // CH2E   CG
	  AA[p].atom[6].fullatom_type =   4; // C      CD
	  AA[p].atom[7].fullatom_type =  17; // OC     OE1
	  AA[p].atom[8].fullatom_type =  17; // OC     OE2
	  AA[p].atom[9].fullatom_type =   0; // H      H
	  AA[p].atom[10].fullatom_type = -1; //    HA
	  AA[p].atom[11].fullatom_type = -1; //   2HB
	  AA[p].atom[12].fullatom_type = -1; //   3HB
	  AA[p].atom[13].fullatom_type = -1; //   2HG
	  AA[p].atom[14].fullatom_type = -1; //   3HG

	  AA[p].atom[15].fullatom_type = 17; //OC      OXT
  	  AA[p].atom[16].fullatom_type =  1; //HC      Hpol
  	  AA[p].atom[17].fullatom_type =  1; //HC      Hpol
	  break;
  	case Sybil:
	  AA[p].atom[0].fullatom_type = 40; // Nbb    N
	  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 10; // CObb   C
	  AA[p].atom[3].fullatom_type = 56; // OCbb   O
	  AA[p].atom[4].fullatom_type = 8; // CH2    CB
	  AA[p].atom[5].fullatom_type = 8; // CH2    CG
	  AA[p].atom[6].fullatom_type = 15; // COO    CD
	  AA[p].atom[7].fullatom_type = 57; // OOC    OE1
	  AA[p].atom[8].fullatom_type = 57; // OOC    OE2
	  AA[p].atom[9].fullatom_type = 66; // HNbb   H
	  AA[p].atom[10].fullatom_type = 66; // Hapo   HA
	  AA[p].atom[11].fullatom_type = 66; // Hapo  2HB
	  AA[p].atom[12].fullatom_type = 66; // Hapo  3HB
	  AA[p].atom[13].fullatom_type = 66; // Hapo  2HG
	  AA[p].atom[14].fullatom_type = 66; // Hapo  3HG

	  AA[p].atom[15].fullatom_type = 57; // OXT
  	  AA[p].atom[16].fullatom_type = 66; // Hpol
  	  AA[p].atom[17].fullatom_type = 66; // Hpol
	  break;
  	default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 4; // CH2    CG
	  AA[p].atom[6].fullatom_type = 2; // COO    CD
	  AA[p].atom[7].fullatom_type = 15; // OOC    OE1
	  AA[p].atom[8].fullatom_type = 15; // OOC    OE2
	  AA[p].atom[9].fullatom_type = 25; // HNbb   H
	  AA[p].atom[10].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[11].fullatom_type = 23; // Hapo  2HB
	  AA[p].atom[12].fullatom_type = 23; // Hapo  3HB
	  AA[p].atom[13].fullatom_type = 23; // Hapo  2HG
	  AA[p].atom[14].fullatom_type = 23; // Hapo  3HG

	  AA[p].atom[15].fullatom_type = 20; // OXT
  	  AA[p].atom[16].fullatom_type = 22; // Hpol
  	  AA[p].atom[17].fullatom_type = 22; // Hpol
	  break;
  }



  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 9; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 10; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--CG
  AA[p].atom[4].bonded_neighbor[2] = 11; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 12; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 4; // CG
  AA[p].atom[5].bonded_neighbor[0] = 4; // CG--CB
  AA[p].atom[5].bonded_neighbor[1] = 6; // CG--CD
  AA[p].atom[5].bonded_neighbor[2] = 13; // CG--2HG
  AA[p].atom[5].bonded_neighbor[3] = 14; // CG--3HG
  AA[p].atom[6].nbonded_neighbors = 3; // CD
  AA[p].atom[6].bonded_neighbor[0] = 5; // CD--CG
  AA[p].atom[6].bonded_neighbor[1] = 7; // CD--OE1
  AA[p].atom[6].bonded_neighbor[2] = 8; // CD--OE2
  AA[p].atom[7].nbonded_neighbors = 1; // OE1
  AA[p].atom[7].bonded_neighbor[0] = 6; // OE1--CD
  AA[p].atom[8].nbonded_neighbors = 1; // OE2
  AA[p].atom[8].bonded_neighbor[0] = 6; // OE2--CD
  AA[p].atom[9].nbonded_neighbors = 1; // H
  AA[p].atom[9].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[10].nbonded_neighbors = 1; // HA
  AA[p].atom[10].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[11].nbonded_neighbors = 1; //2HB
  AA[p].atom[11].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[12].nbonded_neighbors = 1; //3HB
  AA[p].atom[12].bonded_neighbor[0] = 4; //3HB--CB
  AA[p].atom[13].nbonded_neighbors = 1; //2HG
  AA[p].atom[13].bonded_neighbor[0] = 5; //2HG--CG
  AA[p].atom[14].nbonded_neighbors = 1; //3HG
  AA[p].atom[14].bonded_neighbor[0] = 5; //3HG--CG

  AA[p].atom[15].nbonded_neighbors = 1; // OXT
  AA[p].atom[15].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[16].nbonded_neighbors = 1; // 1H
  AA[p].atom[16].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[17].nbonded_neighbors = 1; // 2H
  AA[p].atom[17].bonded_neighbor[0] = 0; // H-N


  for ( k = 0; k < 10; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 10; k < 15; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template for building  HA
  AA[p].atom[10].ta[0] = 1; //  CA
  AA[p].atom[10].ta[1] = 0; //  N
  AA[p].atom[10].ta[2] = 2; //  C

  //bk   template for building 2HB
  AA[p].atom[11].ta[0] = 4; //  CB
  AA[p].atom[11].ta[1] = 1; //  CA
  AA[p].atom[11].ta[2] = 5; //  CG

  //bk   template for building 3HB
  AA[p].atom[12].ta[0] = 4; //  CB
  AA[p].atom[12].ta[1] = 1; //  CA
  AA[p].atom[12].ta[2] = 5; //  CG

  //bk   template for building 2HG
  AA[p].atom[13].ta[0] = 5; //  CG
  AA[p].atom[13].ta[1] = 4; //  CB
  AA[p].atom[13].ta[2] = 6; //  CD

  //bk   template for building 3HG
  AA[p].atom[14].ta[0] = 5; //  CG
  AA[p].atom[14].ta[1] = 4; //  CB
  AA[p].atom[14].ta[2] = 6; //  CD


  //bk   chi angles required to build atoms GLU
  //bk   chi angles needed for building  CG
  AA[p].chi_required[0] [5] = true;

  //bk   chi angles needed for building  CD
  AA[p].chi_required[0] [6] = true;
  AA[p].chi_required[1] [6] = true;

  //bk   chi angles needed for building  OE1
  AA[p].chi_required[0] [7] = true;
  AA[p].chi_required[1] [7] = true;
  AA[p].chi_required[2] [7] = true;

  //bk   chi angles needed for building  OE2
  AA[p].chi_required[0] [8] = true;
  AA[p].chi_required[1] [8] = true;
  AA[p].chi_required[2] [8] = true;

  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [11] = true;

  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [12] = true;

  //bk   chi angles needed for building 2HG
  AA[p].chi_required[0] [13] = true;
  AA[p].chi_required[1] [13] = true;

  //bk   chi angles needed for building 3HG
  AA[p].chi_required[0] [14] = true;
  AA[p].chi_required[1] [14] = true;




  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   CG
  //bk   four atoms that define chi angle  2
  AA[p].chi_atoms[1] [0] = 1; //   CA
  AA[p].chi_atoms[1] [1] = 4; //   CB
  AA[p].chi_atoms[1] [2] = 5; //   CG
  AA[p].chi_atoms[1] [3] = 6; //   CD
  //bk   four atoms that define chi angle  3
  AA[p].chi_atoms[2] [0] = 4; //   CB
  AA[p].chi_atoms[2] [1] = 5; //   CG
  AA[p].chi_atoms[2] [2] = 6; //   CD
  AA[p].chi_atoms[2] [3] = 7; //   OE1


  AA[p].atom[0].icoor[0] = 0.0000; //   N
  AA[p].atom[0].icoor[1] = 7.3000; //   N
  AA[p].atom[0].icoor[2] = 12.6660; //   N
  AA[p].atom[1].icoor[0] = 0.0000; //   CA
  AA[p].atom[1].icoor[1] = 7.2080; //   CA
  AA[p].atom[1].icoor[2] = 14.1210; //   CA
  AA[p].atom[2].icoor[0] = 0.0000; //   C
  AA[p].atom[2].icoor[1] = 8.5910; //   C
  AA[p].atom[2].icoor[2] = 14.7610; //   C
  AA[p].atom[3].icoor[0] = -0.0030; //   O
  AA[p].atom[3].icoor[1] = 9.6070; //   O
  AA[p].atom[3].icoor[2] = 14.0660; //   O
  AA[p].atom[4].icoor[0] = 1.2080; //   CB
  AA[p].atom[4].icoor[1] = 6.4040; //   CB
  AA[p].atom[4].icoor[2] = 14.6070; //   CB
  AA[p].atom[5].icoor[0] = 1.2950; //   CG
  AA[p].atom[5].icoor[1] = 6.2520; //   CG
  AA[p].atom[5].icoor[2] = 16.1190; //   CG
  AA[p].atom[6].icoor[0] = 2.5040; //   CD
  AA[p].atom[6].icoor[1] = 5.4530; //   CD
  AA[p].atom[6].icoor[2] = 16.5190; //   CD
  AA[p].atom[7].icoor[0] = 2.6990; //   OE1
  AA[p].atom[7].icoor[1] = 5.2540; //   OE1
  AA[p].atom[7].icoor[2] = 17.6940; //   OE1
  AA[p].atom[8].icoor[0] = 3.1710; //   OE2
  AA[p].atom[8].icoor[1] = 4.9460; //   OE2
  AA[p].atom[8].icoor[2] = 15.6480; //   OE2
  AA[p].atom[9].icoor[0] = 0.0400; //   H
  AA[p].atom[9].icoor[1] = 8.2170; //   H
  AA[p].atom[9].icoor[2] = 12.2440; //   H
  AA[p].atom[10].icoor[0] = -0.9090; //   HA
  AA[p].atom[10].icoor[1] = 6.7100; //   HA
  AA[p].atom[10].icoor[2] = 14.4600; //   HA
  AA[p].atom[11].icoor[0] = 1.1400; //  2HB
  AA[p].atom[11].icoor[1] = 5.4180; //  2HB
  AA[p].atom[11].icoor[2] = 14.1470; //  2HB
  AA[p].atom[12].icoor[0] = 2.0990; //  3HB
  AA[p].atom[12].icoor[1] = 6.9150; //  3HB
  AA[p].atom[12].icoor[2] = 14.2420; //  3HB
  AA[p].atom[13].icoor[0] = 1.2970; //  2HG
  AA[p].atom[13].icoor[1] = 7.2030; //  2HG
  AA[p].atom[13].icoor[2] = 16.6510; //  2HG
  AA[p].atom[14].icoor[0] = 0.3930; //  3HG
  AA[p].atom[14].icoor[1] = 5.6970; //  3HG
  AA[p].atom[14].icoor[2] = 16.3730; //  3HG

// atom number for backbone HN
 AA[p].HNpos=9;
 // atom number for backbone HA
 AA[p].HApos=10;

 // number of polar hydrogens
 AA[p].nH_polar_complete=1;
 // atom numbers for polar H
 AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
 AA[p].Hpos_polar_complete[0]=9;

 // number of aromatic hydrogens
 AA[p].nH_aromatic_complete=0;

 // number of apolar hydrogens
 AA[p].nH_apolar_complete=5;
 //atom number for apolar hydrogens
 AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
 AA[p].Hpos_apolar_complete[0]=10;
 AA[p].Hpos_apolar_complete[1]=11;
 AA[p].Hpos_apolar_complete[2]=12;
 AA[p].Hpos_apolar_complete[3]=13;
 AA[p].Hpos_apolar_complete[4]=14;

 // number of acceptors
 AA[p].nacceptors=3;
 //acceptor information
 AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
 AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
 (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
 (AA[p].aBase[1])[0]=7;(AA[p].aBase[1])[1]=6;
 (AA[p].aBase[2])[0]=8;(AA[p].aBase[2])[1]=6;
 (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;
 (AA[p].accpt_pos[1])[0]=1;(AA[p].accpt_pos[1])[1]=7;
 (AA[p].accpt_pos[2])[0]=2;(AA[p].accpt_pos[2])[1]=8;

 //Number of Hydrogens connection
 AA[p].nH_hydrogen_connexions=6;
 //atoms hydrogens are connected too
 AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
 AA[p].Hbase[0][0]=9; AA[p].Hbase[0][1]=0;
 AA[p].Hbase[1][0]=10; AA[p].Hbase[1][1]=1;
 AA[p].Hbase[2][0]=11; AA[p].Hbase[2][1]=4;
 AA[p].Hbase[3][0]=12; AA[p].Hbase[3][1]=4;
 AA[p].Hbase[4][0]=13; AA[p].Hbase[4][1]=5;
 AA[p].Hbase[5][0]=14; AA[p].Hbase[5][1]=5;

  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=9;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=10;

  AA[p].atom[4].numHydrogens_atm=2;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=11;
  AA[p].atom[4].hydrogens_atm[1]=12;

  AA[p].atom[5].numHydrogens_atm=2;
  AA[p].atom[5].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[5].numHydrogens_atm);
  AA[p].atom[5].hydrogens_atm[0]=13;
  AA[p].atom[5].hydrogens_atm[1]=14;

}

//FIN DE AMINOACIDO GLUTAMATO

/**
* Initializes Phenylalanine
*
* @param opt: Rosseta or ICM
*/
void init_PHE(Convention opt)
{

  //AMINOACIDO Phenylalanine
  int k, p = PHE;

  strcpy( AA[p].aa_name3, "PHE" );
  AA[p].aa_name1 = 'F';
  AA[p].mass = 147.178;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = false;
  AA[p].aa_is_nonpolar = true;
  AA[p].aa_is_aromatic = true;
  AA[p].aa_is_charged = false;

  AA[p].natoms = 20;
  AA[p].nheavyatoms = 11;
  AA[p].nchi = 2;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );
  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;


  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

  AA[p].chi_types = ( int * ) new int[AA[p].nchi];
  AA[p].chi_types[0] = 2;
  AA[p].chi_types[1] = 4;
  AA[p].natoms_EEF1 = 11;


  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " CG " );
  strcpy( AA[p].atom[6].atom_name, " CD1" );
  strcpy( AA[p].atom[7].atom_name, " CD2" );
  strcpy( AA[p].atom[8].atom_name, " CE1" );
  strcpy( AA[p].atom[9].atom_name, " CE2" );
  strcpy( AA[p].atom[10].atom_name, " CZ " );
  strcpy( AA[p].atom[11].atom_name, " H  " );
  strcpy( AA[p].atom[12].atom_name, " HD1" );
  strcpy( AA[p].atom[13].atom_name, " HE1" );
  strcpy( AA[p].atom[14].atom_name, " HZ " );
  strcpy( AA[p].atom[15].atom_name, " HE2" );
  strcpy( AA[p].atom[16].atom_name, " HD2" );
  strcpy( AA[p].atom[17].atom_name, " HA " );
  strcpy( AA[p].atom[18].atom_name, "2HB " );
  strcpy( AA[p].atom[19].atom_name, "3HB " );

  strcpy( AA[p].atom[20].atom_name," OXT" );
  strcpy( AA[p].atom[21].atom_name,"2H  " );
  strcpy( AA[p].atom[22].atom_name,"3H  " );


  switch(opt)
  {
    case Rosseta:
    case Sybil:
     AA[p].atom[0].charge = -0.47; //    N
      AA[p].atom[1].charge = 0.07; //    CA  <--- Nt(0.21)
      AA[p].atom[2].charge = 0.51; //    C
      AA[p].atom[3].charge = -0.51; //    O
      AA[p].atom[4].charge = -0.18; //    CB
      AA[p].atom[5].charge = 0.00; //    CG
      AA[p].atom[6].charge = -0.115; //    CD1
      AA[p].atom[7].charge = -0.115; //    CD2
      AA[p].atom[8].charge = -0.115; //    CE1
      AA[p].atom[9].charge = -0.115; //    CE2
      AA[p].atom[10].charge = -0.115; //    CZ
      AA[p].atom[11].charge = 0.31; //    H
      AA[p].atom[12].charge = 0.115; //    HD1
      AA[p].atom[13].charge = 0.115; //    HE1
      AA[p].atom[14].charge = 0.115; //    HZ
      AA[p].atom[15].charge = 0.115; //    HE2
      AA[p].atom[16].charge = 0.115; //    HD2
      AA[p].atom[17].charge = 0.09; //    HA
      AA[p].atom[18].charge = 0.09; //    2HB
      AA[p].atom[19].charge = 0.09; //    3HB

      AA[p].atom[20].charge = -0.67; //   0XT
      AA[p].atom[21].charge = 0.33; //  1H
      AA[p].atom[22].charge = 0.33; //  2H
      break;
    case ICM:
      AA[p].atom[0].charge = -0.356; //    N
      AA[p].atom[1].charge = 0.064; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = -0.04; //    CB
      AA[p].atom[5].charge = 0.02; //    CG
      AA[p].atom[6].charge = -0.015; //    CD1
      AA[p].atom[7].charge = -0.015; //    CD2
      AA[p].atom[8].charge = -0.015; //    CE1
      AA[p].atom[9].charge = -0.015; //    CE2
      AA[p].atom[10].charge = -0.015; //    CZ
      AA[p].atom[11].charge = 0.176; //    H
      AA[p].atom[12].charge = 0.01; //    HD1
      AA[p].atom[13].charge = 0.01; //    HE1
      AA[p].atom[14].charge = 0.005; //    HZ
      AA[p].atom[15].charge = 0.010; //    HE2
      AA[p].atom[16].charge = 0.010; //    HD2
      AA[p].atom[17].charge = 0.02; //    HA
      AA[p].atom[18].charge = 0.025; //    2HB
      AA[p].atom[19].charge = 0.025; //    3HB
      break;
  case EEF1:
      AA[p].atom[0].charge = -0.35; //    N
      AA[p].atom[1].charge =  0.10; //    CA
      AA[p].atom[2].charge =  0.55; //    C
      AA[p].atom[3].charge = -0.55; //    O
      AA[p].atom[4].charge =  0.00; //    CB
      AA[p].atom[5].charge =  0.00; //    CG
      AA[p].atom[6].charge =  0.00; //    CD1
      AA[p].atom[7].charge =  0.00; //    CD2
      AA[p].atom[8].charge =  0.00; //    CE1
      AA[p].atom[9].charge =  0.00; //    CE2
      AA[p].atom[10].charge = 0.00; //    CZ
      AA[p].atom[11].charge = 0.25; //    H
      AA[p].atom[12].charge = 0.00; //    HD1
      AA[p].atom[13].charge = 0.00; //    HE1
      AA[p].atom[14].charge = 0.00; //    HZ
      AA[p].atom[15].charge = 0.00; //    HE2
      AA[p].atom[16].charge = 0.00; //    HD2
      AA[p].atom[17].charge = 0.00; //    HA
      AA[p].atom[18].charge = 0.00; //    2HB
      AA[p].atom[19].charge = 0.00; //    3HB
      break;

  }

/*
   1.0 PHE-N:2H    0.33000 HC
   1.0 PHE-N:CE2  -0.11500 CA
   1.0 PHE-N:HA    0.10000 HB
   1.0 PHE-N:3H    0.33000 HC
   1.0 PHE-N:HE2   0.11500 HP
   1.0 PHE-N:CB   -0.18000 CT2
   1.0 PHE-N:C     0.51000 C
   1.0 PHE-N:HB1   0.09000 HA
   1.0 PHE-N:1HB   0.09000 HA
   1.0 PHE-N:CG    0.00000 CA
   1.0 PHE-N:CD1  -0.11500 CA
   1.0 PHE-N:HD1   0.11500 HP
   1.0 PHE-N:O    -0.51000 O
   1.0 PHE-N:CE1  -0.11500 CA
   1.0 PHE-N:HE1   0.11500 HP
   1.0 PHE-N:CZ   -0.11500 CA
   1.0 PHE-N:HZ    0.11500 HP
   1.0 PHE-N:N    -0.30000 NH3
   1.0 PHE-N:CD2  -0.11500 CA
   1.0 PHE-N:1H    0.33000 HC
   1.0 PHE-N:HD2   0.11500 HP
   1.0 PHE-N:CA    0.21000 CT1

   1.0 PHE-C:CE2  -0.11500 CA
   1.0 PHE-C:HA    0.09000 HB
   1.0 PHE-C:HE2   0.11500 HP
   1.0 PHE-C:CB   -0.18000 CT2
   1.0 PHE-C:C     0.34000 CC
   1.0 PHE-C:HB1   0.09000 HA
   1.0 PHE-C:1HB   0.09000 HA
   1.0 PHE-C:CG    0.00000 CA
   1.0 PHE-C:CD1  -0.11500 CA
   1.0 PHE-C:HD1   0.11500 HP
   1.0 PHE-C:CE1  -0.11500 CA
   1.0 PHE-C:O    -0.67000 OC
   1.0 PHE-C:HE1   0.11500 HP
   1.0 PHE-C:OXT  -0.67000 OC
   1.0 PHE-C:CZ   -0.11500 CA
   1.0 PHE-C:HZ    0.11500 HP
   1.0 PHE-C:N    -0.47000 NH1
   1.0 PHE-C:CD2  -0.11500 CA
   1.0 PHE-C:HN    0.31000 H
   1.0 PHE-C:HD2   0.11500 HP
   1.0 PHE-C:CA    0.07000 CT1
*/

  switch ( opt ) {
  	case EEF1:
	  AA[p].atom[0].fullatom_type =  12; // NH1    N
	  AA[p].atom[1].fullatom_type =   5; // CH1E   CA
	  AA[p].atom[2].fullatom_type =   4; // C      C
	  AA[p].atom[3].fullatom_type =  16; // O      O
	  AA[p].atom[4].fullatom_type =   6; // CH2E   CB
	  AA[p].atom[5].fullatom_type =  24; // CR     CG
	  AA[p].atom[6].fullatom_type =   8; // CR1E   CD1
	  AA[p].atom[7].fullatom_type =   8; // CR1E   CD2
	  AA[p].atom[8].fullatom_type =   8; // CR1E   CE1
	  AA[p].atom[9].fullatom_type =   8; // CR1E   CE2
	  AA[p].atom[10].fullatom_type =  8; // CR1E   CZ
	  AA[p].atom[11].fullatom_type =  0; // H      H
	  AA[p].atom[12].fullatom_type = -1; // Haro   HD1
	  AA[p].atom[13].fullatom_type = -1; // Haro   HE1
	  AA[p].atom[14].fullatom_type = -1; // Haro   HZ
	  AA[p].atom[15].fullatom_type = -1; // Haro   HE2
	  AA[p].atom[16].fullatom_type = -1; // Haro   HD2
	  AA[p].atom[17].fullatom_type = -1; // Hapo   HA
	  AA[p].atom[18].fullatom_type = -1; // Hapo   2HB
	  AA[p].atom[19].fullatom_type = -1; // Hapo   3HB

	  AA[p].atom[20].fullatom_type = 17; // OC     OXT
	  AA[p].atom[21].fullatom_type =  1; // HC     Hpol
	  AA[p].atom[22].fullatom_type =  1; // HC     Hpol
	  break;
  	case Sybil:
  		  AA[p].atom[0].fullatom_type = 40; // Nbb    N
  		  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
  		  AA[p].atom[2].fullatom_type = 10; // CObb   C
  		  AA[p].atom[3].fullatom_type = 56; // OCbb   O
  		  AA[p].atom[4].fullatom_type = 8; // CH2    CB
  		  AA[p].atom[5].fullatom_type = 12; // aroC   CG
  		  AA[p].atom[6].fullatom_type = 12; // aroC   CD1
  		  AA[p].atom[7].fullatom_type = 12; // aroC   CD2
  		  AA[p].atom[8].fullatom_type = 12; // aroC   CE1
  		  AA[p].atom[9].fullatom_type = 12; // aroC   CE2
  		  AA[p].atom[10].fullatom_type = 12; // aroC   CZ
  		  AA[p].atom[11].fullatom_type = 66; // HNbb   H
  		  AA[p].atom[12].fullatom_type = 66; // Haro   HD1
  		  AA[p].atom[13].fullatom_type = 66; // Haro   HE1
  		  AA[p].atom[14].fullatom_type = 66; // Haro   HZ
  		  AA[p].atom[15].fullatom_type = 66; // Haro   HE2
  		  AA[p].atom[16].fullatom_type = 66; // Haro   HD2
  		  AA[p].atom[17].fullatom_type = 66; // Hapo   HA
  		  AA[p].atom[18].fullatom_type = 66; // Hapo  2HB
  		  AA[p].atom[19].fullatom_type = 66; // Hapo  3HB

  		  AA[p].atom[20].fullatom_type = 57; // OXT
  		  AA[p].atom[21].fullatom_type = 66; // Hpol
  		  AA[p].atom[22].fullatom_type = 66; // Hpol
  		  break;
  	default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 6; // aroC   CG
	  AA[p].atom[6].fullatom_type = 6; // aroC   CD1
	  AA[p].atom[7].fullatom_type = 6; // aroC   CD2
	  AA[p].atom[8].fullatom_type = 6; // aroC   CE1
	  AA[p].atom[9].fullatom_type = 6; // aroC   CE2
	  AA[p].atom[10].fullatom_type = 6; // aroC   CZ
	  AA[p].atom[11].fullatom_type = 25; // HNbb   H
	  AA[p].atom[12].fullatom_type = 24; // Haro   HD1
	  AA[p].atom[13].fullatom_type = 24; // Haro   HE1
	  AA[p].atom[14].fullatom_type = 24; // Haro   HZ
	  AA[p].atom[15].fullatom_type = 24; // Haro   HE2
	  AA[p].atom[16].fullatom_type = 24; // Haro   HD2
	  AA[p].atom[17].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[18].fullatom_type = 23; // Hapo  2HB
	  AA[p].atom[19].fullatom_type = 23; // Hapo  3HB

	  AA[p].atom[20].fullatom_type = 20; // OXT
	  AA[p].atom[21].fullatom_type = 22; // Hpol
	  AA[p].atom[22].fullatom_type = 22; // Hpol
	  break;
  }

  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 11; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 17; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--CG
  AA[p].atom[4].bonded_neighbor[2] = 18; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 19; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 3; // CG
  AA[p].atom[5].bonded_neighbor[0] = 4; // CG--CB
  AA[p].atom[5].bonded_neighbor[1] = 6; // CG--CD1
  AA[p].atom[5].bonded_neighbor[2] = 7; // CG--CD2
  AA[p].atom[6].nbonded_neighbors = 3; // CD1
  AA[p].atom[6].bonded_neighbor[0] = 5; // CD1--CG
  AA[p].atom[6].bonded_neighbor[1] = 8; // CD1--CE1
  AA[p].atom[6].bonded_neighbor[2] = 12; // CD1--HD1
  AA[p].atom[7].nbonded_neighbors = 3; // CD2
  AA[p].atom[7].bonded_neighbor[0] = 5; // CD2--CG
  AA[p].atom[7].bonded_neighbor[1] = 9; // CD2--CE2
  AA[p].atom[7].bonded_neighbor[2] = 16; // CD2--HD2
  AA[p].atom[8].nbonded_neighbors = 3; // CE1
  AA[p].atom[8].bonded_neighbor[0] = 6; // CE1--CD1
  AA[p].atom[8].bonded_neighbor[1] = 10; // CE1--CZ
  AA[p].atom[8].bonded_neighbor[2] = 13; // CE1--HE1
  AA[p].atom[9].nbonded_neighbors = 3; // CE2
  AA[p].atom[9].bonded_neighbor[0] = 7; // CE2--CD2
  AA[p].atom[9].bonded_neighbor[1] = 10; // CE2--CZ
  AA[p].atom[9].bonded_neighbor[2] = 15; // CE2--HE2
  AA[p].atom[10].nbonded_neighbors = 3; // CZ
  AA[p].atom[10].bonded_neighbor[0] = 8; // CZ--CE1
  AA[p].atom[10].bonded_neighbor[1] = 9; // CZ--CE2
  AA[p].atom[10].bonded_neighbor[2] = 14; // CZ--HZ
  AA[p].atom[11].nbonded_neighbors = 1; // H
  AA[p].atom[11].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[12].nbonded_neighbors = 1; // HD1
  AA[p].atom[12].bonded_neighbor[0] = 6; // HD1--CD1
  AA[p].atom[13].nbonded_neighbors = 1; // HE1
  AA[p].atom[13].bonded_neighbor[0] = 8; // HE1--CE1
  AA[p].atom[14].nbonded_neighbors = 1; // HZ
  AA[p].atom[14].bonded_neighbor[0] = 10; // HZ--CZ
  AA[p].atom[15].nbonded_neighbors = 1; // HE2
  AA[p].atom[15].bonded_neighbor[0] = 9; // HE2--CE2
  AA[p].atom[16].nbonded_neighbors = 1; // HD2
  AA[p].atom[16].bonded_neighbor[0] = 7; // HD2--CD2
  AA[p].atom[17].nbonded_neighbors = 1; // HA
  AA[p].atom[17].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[18].nbonded_neighbors = 1; //2HB
  AA[p].atom[18].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[19].nbonded_neighbors = 1; //3HB
  AA[p].atom[19].bonded_neighbor[0] = 4; //3HB--CB

  AA[p].atom[20].nbonded_neighbors = 1; // OXT
  AA[p].atom[20].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[21].nbonded_neighbors = 1; // 1H
  AA[p].atom[21].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[22].nbonded_neighbors = 1; // 2H
  AA[p].atom[22].bonded_neighbor[0] = 0; // H-N



  for ( k = 0; k < 12; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 12; k < 20; k++ )
    AA[p].atom[k].hastemplate = true;




  //bk   template for building  HD1
  AA[p].atom[12].ta[0] = 6; //   CD1
  AA[p].atom[12].ta[1] = 5; //   CG
  AA[p].atom[12].ta[2] = 8; //   CE1

  //bk   template for building  HE1
  AA[p].atom[13].ta[0] = 8; //   CE1
  AA[p].atom[13].ta[1] = 6; //   CD1
  AA[p].atom[13].ta[2] = 10; //   CZ

  //bk   template for building  HZ
  AA[p].atom[14].ta[0] = 10; //   CZ
  AA[p].atom[14].ta[1] = 8; //   CE1
  AA[p].atom[14].ta[2] = 9; //   CE2

  //bk   template for building  HE2
  AA[p].atom[15].ta[0] = 9; //   CE2
  AA[p].atom[15].ta[1] = 7; //   CD2
  AA[p].atom[15].ta[2] = 10; //   CZ

  //bk   template for building  HD2
  AA[p].atom[16].ta[0] = 7; //   CD2
  AA[p].atom[16].ta[1] = 5; //   CG
  AA[p].atom[16].ta[2] = 9; //   CE2

  //bk   template for building  HA
  AA[p].atom[17].ta[0] = 1; //   CA
  AA[p].atom[17].ta[1] = 0; //   N
  AA[p].atom[17].ta[2] = 2; //   C

  //bk   template for building 2HB
  AA[p].atom[18].ta[0] = 4; //   CB
  AA[p].atom[18].ta[1] = 1; //   CA
  AA[p].atom[18].ta[2] = 5; //   CG

  //bk   template for building 3HB
  AA[p].atom[19].ta[0] = 4; //   CB
  AA[p].atom[19].ta[1] = 1; //   CA
  AA[p].atom[19].ta[2] = 5; //   CG


  //bk   chi angles required to build atoms PHE
  //bk   chi angles needed for building  CG
  AA[p].chi_required[0] [5] = true;

  //bk   chi angles needed for building  CD1
  AA[p].chi_required[0] [6] = true;
  AA[p].chi_required[1] [6] = true;

  //bk   chi angles needed for building  CD2
  AA[p].chi_required[0] [7] = true;
  AA[p].chi_required[1] [7] = true;

  //bk   chi angles needed for building  CE1
  AA[p].chi_required[0] [8] = true;
  AA[p].chi_required[1] [8] = true;

  //bk   chi angles needed for building  CE2
  AA[p].chi_required[0] [9] = true;
  AA[p].chi_required[1] [9] = true;

  //bk   chi angles needed for building  CZ
  AA[p].chi_required[0] [10] = true;
  AA[p].chi_required[1] [10] = true;

  //bk   chi angles needed for building  HD1
  AA[p].chi_required[0] [12] = true;
  AA[p].chi_required[1] [12] = true;



  //bk   chi angles needed for building  HE1
  AA[p].chi_required[0] [13] = true;
  AA[p].chi_required[1] [13] = true;


  //bk   chi angles needed for building  HZ
  AA[p].chi_required[0] [14] = true;
  AA[p].chi_required[1] [14] = true;

  //bk   chi angles needed for building  HE2
  AA[p].chi_required[0] [15] = true;
  AA[p].chi_required[1] [15] = true;

  //bk   chi angles needed for building  HD2
  AA[p].chi_required[0] [16] = true;
  AA[p].chi_required[1] [16] = true;

  //bk   chi angles needed for building  HA

  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [18] = true;

  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [19] = true;



  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   CG
  //bk   four atoms that define chi angle  2
  AA[p].chi_atoms[1] [0] = 1; //   CA
  AA[p].chi_atoms[1] [1] = 4; //   CB
  AA[p].chi_atoms[1] [2] = 5; //   CG
  AA[p].chi_atoms[1] [3] = 6; //   CD1



  AA[p].atom[0].icoor[0] = 0.0000; //   N
  AA[p].atom[0].icoor[1] = 8.6220; //   N
  AA[p].atom[0].icoor[2] = 16.0890; //   N
  AA[p].atom[1].icoor[0] = 0.0000; //   CA
  AA[p].atom[1].icoor[1] = 9.8810; //   CA
  AA[p].atom[1].icoor[2] = 16.8260; //   CA
  AA[p].atom[2].icoor[0] = 0.0000; //   C
  AA[p].atom[2].icoor[1] = 9.6390; //   C
  AA[p].atom[2].icoor[2] = 18.3300; //   C
  AA[p].atom[3].icoor[0] = 0.0020; //   O
  AA[p].atom[3].icoor[1] = 8.4940; //   O
  AA[p].atom[3].icoor[2] = 18.7830; //   O
  AA[p].atom[4].icoor[0] = -1.2080; //   CB
  AA[p].atom[4].icoor[1] = 10.7320; //   CB
  AA[p].atom[4].icoor[2] = 16.4300; //   CB
  AA[p].atom[5].icoor[0] = -1.2820; //   CG
  AA[p].atom[5].icoor[1] = 12.0500; //   CG
  AA[p].atom[5].icoor[2] = 17.1470; //   CG
  AA[p].atom[6].icoor[0] = -2.3200; //   CD1
  AA[p].atom[6].icoor[1] = 12.9350; //   CD1
  AA[p].atom[6].icoor[2] = 16.8960; //   CD1
  AA[p].atom[7].icoor[0] = -0.3140; //   CD2
  AA[p].atom[7].icoor[1] = 12.4080; //   CD2
  AA[p].atom[7].icoor[2] = 18.0730; //   CD2
  AA[p].atom[8].icoor[0] = -2.3890; //   CE1
  AA[p].atom[8].icoor[1] = 14.1480; //   CE1
  AA[p].atom[8].icoor[2] = 17.5550; //   CE1
  AA[p].atom[9].icoor[0] = -0.3810; //   CE2
  AA[p].atom[9].icoor[1] = 13.6190; //   CE2
  AA[p].atom[9].icoor[2] = 18.7340; //   CE2
  AA[p].atom[10].icoor[0] = -1.4200; //   CZ
  AA[p].atom[10].icoor[1] = 14.4900; //   CZ
  AA[p].atom[10].icoor[2] = 18.4740; //   CZ
  AA[p].atom[11].icoor[0] = -0.0400; //   H
  AA[p].atom[11].icoor[1] = 7.7600; //   H
  AA[p].atom[11].icoor[2] = 16.6140; //   H
  AA[p].atom[12].icoor[0] = -3.0870; //   HD1
  AA[p].atom[12].icoor[1] = 12.6640; //   HD1
  AA[p].atom[12].icoor[2] = 16.1700; //   HD1
  AA[p].atom[13].icoor[0] = -3.2100; //   HE1
  AA[p].atom[13].icoor[1] = 14.8340; //   HE1
  AA[p].atom[13].icoor[2] = 17.3480; //   HE1
  AA[p].atom[14].icoor[0] = -1.4730; //   HZ
  AA[p].atom[14].icoor[1] = 15.4460; //   HZ
  AA[p].atom[14].icoor[2] = 18.9930; //   HZ
  AA[p].atom[15].icoor[0] = 0.3860; //   HE2
  AA[p].atom[15].icoor[1] = 13.8880; //   HE2
  AA[p].atom[15].icoor[2] = 19.4600; //   HE2
  AA[p].atom[16].icoor[0] = 0.5070; //   HD2
  AA[p].atom[16].icoor[1] = 11.7200; //   HD2
  AA[p].atom[16].icoor[2] = 18.2780; //   HD2
  AA[p].atom[17].icoor[0] = 0.9090; //   HA
  AA[p].atom[17].icoor[1] = 10.4410; //   HA
  AA[p].atom[17].icoor[2] = 16.6020; //   HA
  AA[p].atom[18].icoor[0] = -1.1750; //  2HB
  AA[p].atom[18].icoor[1] = 10.9580; //  2HB
  AA[p].atom[18].icoor[2] = 15.3650; //  2HB
  AA[p].atom[19].icoor[0] = -2.1310; //  3HB
  AA[p].atom[19].icoor[1] = 10.2010; //  3HB
  AA[p].atom[19].icoor[2] = 16.6590; //  3HB


// atom number for backbone HN
 AA[p].HNpos=11;
 // atom number for backbone HA
 AA[p].HApos=17;

 // number of polar hydrogens
 AA[p].nH_polar_complete=1;
 // atom numbers for polar H
 AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
 AA[p].Hpos_polar_complete[0]=11;

 // number of aromatic hydrogens
 AA[p].nH_aromatic_complete=5;
 AA[p].Hpos_aromatic_complete=(int*)malloc(sizeof(int)*AA[p].nH_aromatic_complete);
 AA[p].Hpos_aromatic_complete[0]=12;
 AA[p].Hpos_aromatic_complete[1]=13;
 AA[p].Hpos_aromatic_complete[2]=14;
 AA[p].Hpos_aromatic_complete[3]=15;
 AA[p].Hpos_aromatic_complete[4]=16;

 // number of apolar hydrogens
 AA[p].nH_apolar_complete=3;
 //atom number for apolar hydrogens
 AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
 AA[p].Hpos_apolar_complete[0]=17;
 AA[p].Hpos_apolar_complete[1]=18;
 AA[p].Hpos_apolar_complete[2]=19;

 // number of acceptors
 AA[p].nacceptors=1;
 //acceptor information
 AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
 AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
 (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
 (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;

 //Number of Hydrogens connection
 AA[p].nH_hydrogen_connexions=9;
 //atoms hydrogens are connected too
 AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
 AA[p].Hbase[0][0]=11; AA[p].Hbase[0][1]=0;
 AA[p].Hbase[1][0]=12; AA[p].Hbase[1][1]=6;
 AA[p].Hbase[2][0]=13; AA[p].Hbase[2][1]=8;
 AA[p].Hbase[3][0]=14; AA[p].Hbase[3][1]=10;
 AA[p].Hbase[4][0]=15; AA[p].Hbase[4][1]=9;
 AA[p].Hbase[5][0]=16; AA[p].Hbase[5][1]=7;
 AA[p].Hbase[6][0]=17; AA[p].Hbase[6][1]=1;
 AA[p].Hbase[7][0]=18; AA[p].Hbase[7][1]=4;
 AA[p].Hbase[8][0]=19; AA[p].Hbase[8][1]=4;

  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=11;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=17;

  AA[p].atom[4].numHydrogens_atm=2;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=18;
  AA[p].atom[4].hydrogens_atm[1]=18;

  AA[p].atom[6].numHydrogens_atm=1;
  AA[p].atom[6].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[6].numHydrogens_atm);
  AA[p].atom[6].hydrogens_atm[0]=12;

  AA[p].atom[7].numHydrogens_atm=1;
  AA[p].atom[7].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[7].numHydrogens_atm);
  AA[p].atom[7].hydrogens_atm[0]=16;

  AA[p].atom[8].numHydrogens_atm=1;
  AA[p].atom[8].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[8].numHydrogens_atm);
  AA[p].atom[8].hydrogens_atm[0]=13;

  AA[p].atom[9].numHydrogens_atm=1;
  AA[p].atom[9].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[9].numHydrogens_atm);
  AA[p].atom[9].hydrogens_atm[0]=15;

  AA[p].atom[10].numHydrogens_atm=1;
  AA[p].atom[10].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[10].numHydrogens_atm);
  AA[p].atom[10].hydrogens_atm[0]=14;

  //FIN DE AMINOACIDO Phenylalanine

}

/**
* Initializes Glycine
*
* @param opt: Rosseta or ICM
*/
void init_GLY(Convention opt)
{
  //AMINOACIDO GLYCINE
  int k, p =GLY;

  strcpy( AA[p].aa_name3, "GLY" );
  AA[p].aa_name1 = 'G';
  AA[p].mass = 57.053;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = false;
  AA[p].aa_is_nonpolar = true;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = false;

  AA[p].natoms = 7;
  AA[p].nheavyatoms = 4;
  AA[p].nchi = 0;

  AA[p].chi_atoms = NULL;
  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

  AA[p].chi_types = NULL;
  AA[p].natoms_EEF1 = 4;

  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " H  " );
  strcpy( AA[p].atom[5].atom_name, "2HA " );
  strcpy( AA[p].atom[6].atom_name, "3HA " );

  strcpy( AA[p].atom[7].atom_name," OXT" );
  strcpy( AA[p].atom[8].atom_name,"2H  " );
  strcpy( AA[p].atom[9].atom_name,"3H  " );


  switch(opt)
  {
    case Rosseta:
    case Sybil:
    AA[p].atom[0].charge = -0.47; //    N <---Nt(-0.30)
      AA[p].atom[1].charge = -0.02; //    CA <---Nt(0.13)
      AA[p].atom[2].charge = 0.51; //    C <---Ct(0.34)
      AA[p].atom[3].charge = -0.51; //    O <---Ct(0.67)
      AA[p].atom[4].charge = 0.31; //    H <---Nt(0.33)
      AA[p].atom[5].charge = 0.09; //    2HA
      AA[p].atom[6].charge = 0.09; //    3HA

      AA[p].atom[7].charge = -0.67; //   0XT
      AA[p].atom[8].charge = 0.33; //  1H
      AA[p].atom[9].charge = 0.33; //  2H
      break;

    case ICM:
      AA[p].atom[0].charge = -0.344; //    N
      AA[p].atom[1].charge = -0.008; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = 0.176; //    H
      AA[p].atom[5].charge = 0.055; //    2HA
      AA[p].atom[6].charge = 0.055; //    3HA

      AA[p].atom[7].charge = -0.67; //   0XT
      AA[p].atom[8].charge = 0.33; //  1H
      AA[p].atom[9].charge = 0.33; //  2H
      break;
    case EEF1:
      AA[p].atom[0].charge = -0.35; //    N
      AA[p].atom[1].charge =  0.10; //    CA
      AA[p].atom[2].charge =  0.55; //    C
      AA[p].atom[3].charge = -0.55; //    O
      AA[p].atom[4].charge =  0.25; //    H
      AA[p].atom[5].charge =  0.00; //    2HA
      AA[p].atom[6].charge =  0.00; //    3HA
      break;
  }
  /*
   1.0 GLY-N:2H    0.33000 HC
   1.0 GLY-N:HA1   0.09000 HB
   1.0 GLY-N:3H    0.33000 HC   *
   1.0 GLY-N:1HA   0.09000 HB
   1.0 GLY-N:C     0.51000 C
   1.0 GLY-N:O    -0.51000 O
   1.0 GLY-N:N    -0.30000 NH3  *
   1.0 GLY-N:1H    0.33000 HC
   1.0 GLY-N:CA    0.13000 CT2  *

   1.0 GLY-C:HA1   0.09000 HB
   1.0 GLY-C:1HA   0.09000 HB
   1.0 GLY-C:C     0.34000 CC   *
   1.0 GLY-C:O    -0.67000 OC   *
   1.0 GLY-C:OXT  -0.67000 OC
   1.0 GLY-C:N    -0.47000 NH1
   1.0 GLY-C:HN    0.31000 H
   1.0 GLY-C:CA   -0.02000 CT2
*/

  switch (opt) {
  	case EEF1:
	  AA[p].atom[0].fullatom_type = 12; // NH1    N
	  AA[p].atom[1].fullatom_type =  6; // CH2E   CA
	  AA[p].atom[2].fullatom_type =  4; // C      C
	  AA[p].atom[3].fullatom_type = 16; // O      O
	  AA[p].atom[4].fullatom_type =  0; // H      H
	  AA[p].atom[5].fullatom_type = -1; // HA     2HA
	  AA[p].atom[6].fullatom_type = -1; // HA     3HA

	  AA[p].atom[7].fullatom_type = 17; // OC     OXT
	  AA[p].atom[8].fullatom_type =  1; // Hpol   2H
	  AA[p].atom[9].fullatom_type =  1; // Hpol   3H
	  break;
  	case Sybil:
  		  AA[p].atom[0].fullatom_type = 40; // Nbb    N
  		  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
  		  AA[p].atom[2].fullatom_type = 10; // CObb   C
  		  AA[p].atom[3].fullatom_type = 56; // OCbb   O
  		  AA[p].atom[4].fullatom_type = 66; // HNbb   H
  		  AA[p].atom[5].fullatom_type = 66; // Hapo   2HA
  		  AA[p].atom[6].fullatom_type = 66; // Hapo   3HA

  		  AA[p].atom[7].fullatom_type = 57; // OXT
  		  AA[p].atom[8].fullatom_type = 66; // Hpol
  		  AA[p].atom[9].fullatom_type = 66; // Hpol
  		  break;
  	default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 25; // HNbb   H
	  AA[p].atom[5].fullatom_type = 23; // Hapo   2HA
	  AA[p].atom[6].fullatom_type = 23; // Hapo   3HA

	  AA[p].atom[7].fullatom_type = 20; // OXT
	  AA[p].atom[8].fullatom_type = 22; // Hpol
	  AA[p].atom[9].fullatom_type = 22; // Hpol
	  break;
  }


  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 4; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 5; // CA--2HA
  AA[p].atom[1].bonded_neighbor[3] = 6; // CA--3HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 1; // H
  AA[p].atom[4].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[5].nbonded_neighbors = 1; // 2HA
  AA[p].atom[5].bonded_neighbor[0] = 1; // 2HA--CA
  AA[p].atom[6].nbonded_neighbors = 1; // 3HA
  AA[p].atom[6].bonded_neighbor[0] = 1; // 3HA--CA

  AA[p].atom[7].nbonded_neighbors = 1; // OXT
  AA[p].atom[7].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[8].nbonded_neighbors = 1; // 1H
  AA[p].atom[8].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[9].nbonded_neighbors = 1; // 2H
  AA[p].atom[9].bonded_neighbor[0] = 0; // H-N


  for ( k = 0; k < 5; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 5; k < 7; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template for building  2HA
  AA[p].atom[5].ta[0] = 1; //   CA
  AA[p].atom[5].ta[1] = 0; //   N
  AA[p].atom[5].ta[2] = 2; //   C

  //bk   template for building  3HA
  AA[p].atom[6].ta[0] = 1; //   CA
  AA[p].atom[6].ta[1] = 0; //   N
  AA[p].atom[6].ta[2] = 3; //   C



  AA[p].atom[0].icoor[0] = 0.0000; //   N
  AA[p].atom[0].icoor[1] = 10.7220; //   N
  AA[p].atom[0].icoor[2] = 19.0990; //   N
  AA[p].atom[1].icoor[0] = 0.0000; //   CA
  AA[p].atom[1].icoor[1] = 10.6280; //   CA
  AA[p].atom[1].icoor[2] = 20.5540; //   CA
  AA[p].atom[2].icoor[0] = 0.0000; //   C
  AA[p].atom[2].icoor[1] = 12.0090; //   C
  AA[p].atom[2].icoor[2] = 21.1960; //   C
  AA[p].atom[3].icoor[0] = -0.0030; //   O
  AA[p].atom[3].icoor[1] = 13.0390; //   O
  AA[p].atom[3].icoor[2] = 20.5210; //   O
  AA[p].atom[4].icoor[0] = 0.0000; //   H
  AA[p].atom[4].icoor[1] = 11.6400; //   H
  AA[p].atom[4].icoor[2] = 18.6790; //   H
  AA[p].atom[5].icoor[0] = -0.8900; //  2HA
  AA[p].atom[5].icoor[1] = 10.0870; //  2HA
  AA[p].atom[5].icoor[2] = 20.8760; //  2HA
  AA[p].atom[6].icoor[0] = 0.8890; //  3HA
  AA[p].atom[6].icoor[1] = 10.0870; //  3HA
  AA[p].atom[6].icoor[2] = 20.8760; //  3HA

  // atom number for backbone HN
  AA[p].HNpos=4;
 // atom number for backbone HA
   AA[p].HApos=6;

 // number of polar hydrogens
 AA[p].nH_polar_complete=1;
 // atom numbers for polar H
 AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
 AA[p].Hpos_polar_complete[0]=4;

 // number of aromatic hydrogens
 AA[p].nH_aromatic_complete=0;

 // number of apolar hydrogens
 AA[p].nH_apolar_complete=2;
 //atom number for apolar hydrogens
 AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
 AA[p].Hpos_apolar_complete[0]=5;
 AA[p].Hpos_apolar_complete[1]=6;

 // number of acceptors
 AA[p].nacceptors=1;
 //acceptor information
 AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
 AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
 (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
 (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;

 //Number of Hydrogens connection
 AA[p].nH_hydrogen_connexions=3;
 //atoms hydrogens are connected too
 AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
 AA[p].Hbase[0][0]=4; AA[p].Hbase[0][1]=0;
 AA[p].Hbase[1][0]=5; AA[p].Hbase[1][1]=1;
 AA[p].Hbase[2][0]=6; AA[p].Hbase[2][1]=1;

  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=4;

  AA[p].atom[1].numHydrogens_atm=2;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=5;
  AA[p].atom[1].hydrogens_atm[1]=6;

}

//FIN DE AMINOACIDO GLYCINE

/**
* Initializes Histidine
*
* @param opt: Rosseta or ICM
*/
void init_HIS(Convention opt)
{
  //AMINOACIDO HISTIDINE
  int k, p =HIS;


  strcpy( AA[p].aa_name3, "HIS" ); // <--- HIS = HSD (Protonated at NE2, neutral)
  AA[p].aa_name1 = 'H';            // Note: HID (Protonated at ND1, neutral)
  AA[p].mass = 137.143;            // Note: HIP (Di-Protonated at ND1 & NE2, +)

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = true;
  AA[p].aa_is_nonpolar = false;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = false;

  AA[p].natoms = 17;
  AA[p].nheavyatoms = 10;
  AA[p].nchi = 2;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );
  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;

  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

  AA[p].chi_types = ( int * ) new int[AA[p].nchi];
  AA[p].chi_types[0] = 2;
  AA[p].chi_types[1] = 4;
  AA[p].natoms_EEF1 = 11;


  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " CG " );
  strcpy( AA[p].atom[6].atom_name, " ND1" );
  strcpy( AA[p].atom[7].atom_name, " CD2" );
  strcpy( AA[p].atom[8].atom_name, " CE1" );
  strcpy( AA[p].atom[9].atom_name, " NE2" ); // <--- Protonated here!
  strcpy( AA[p].atom[10].atom_name," H  " );
  strcpy( AA[p].atom[11].atom_name," HE2" );// <--- The proton
  strcpy( AA[p].atom[12].atom_name," HA " );
  strcpy( AA[p].atom[13].atom_name,"2HB " );
  strcpy( AA[p].atom[14].atom_name,"3HB " );
  strcpy( AA[p].atom[15].atom_name," HE1" );
  strcpy( AA[p].atom[16].atom_name," HD2" );

  strcpy( AA[p].atom[17].atom_name," OXT" );
  strcpy( AA[p].atom[18].atom_name,"2H  " );
  strcpy( AA[p].atom[19].atom_name,"3H  " );


  switch(opt)
  {
    case Rosseta:
    case Sybil:
     AA[p].atom[0].charge = -0.47; //    N
      AA[p].atom[1].charge = 0.07; //    CA
      AA[p].atom[2].charge = 0.51; //    C
      AA[p].atom[3].charge = -0.51; //    O
      AA[p].atom[4].charge = -0.08; //    CB
      AA[p].atom[5].charge = 0.22; //    CG
      AA[p].atom[6].charge = -0.70; //    ND1
      AA[p].atom[7].charge = -0.05; //    CD2
      AA[p].atom[8].charge = 0.25; //    CE1  (¿HE1?)
      AA[p].atom[9].charge = -0.36; //    NE2
      AA[p].atom[10].charge = 0.31; //    H
      AA[p].atom[11].charge = 0.32; //    HE2
      AA[p].atom[12].charge = 0.09; //    HA
      AA[p].atom[13].charge = 0.09; //    2HB
      AA[p].atom[14].charge = 0.09; //    3HB
      AA[p].atom[15].charge = 0.13; //    HE1
      AA[p].atom[16].charge = 0.09; //    HD2

      AA[p].atom[17].charge = -0.67; //   0XT
      AA[p].atom[18].charge = 0.33; //  1H
      AA[p].atom[19].charge = 0.33; //  2H
      break;

    case ICM:
      AA[p].atom[0].charge = -0.356; //    N
      AA[p].atom[1].charge = 0.064; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = -0.05; //    CB
      AA[p].atom[5].charge = 0.15; //    CG
      AA[p].atom[6].charge = -0.2; //    ND1
      AA[p].atom[7].charge = 0.1; //    CD2
      AA[p].atom[8].charge = 0.275; //    CE1
      AA[p].atom[9].charge = -0.2; //    NE2
      AA[p].atom[10].charge = 0.176; //    H
      AA[p].atom[11].charge = 0.265; //    HE2
      AA[p].atom[12].charge = 0.02; //    HA
      AA[p].atom[13].charge = 0.065; //    2HB
      AA[p].atom[14].charge = 0.065; //    3HB
      AA[p].atom[15].charge = 0.15; //    HE1
      AA[p].atom[16].charge = 0.135; //    HD2

      AA[p].atom[17].charge = -0.67; //   0XT
      AA[p].atom[18].charge = 0.33; //  1H
      AA[p].atom[19].charge = 0.33; //  2H
      break;

   case EEF1:
      AA[p].atom[0].charge = -0.35; //    N
      AA[p].atom[1].charge =  0.10; //    CA
      AA[p].atom[2].charge =  0.55; //    C
      AA[p].atom[3].charge = -0.55; //    O
      AA[p].atom[4].charge =  0.00; //    CB
      AA[p].atom[5].charge =  0.10; //    CG
      AA[p].atom[6].charge = -0.40; //    ND1
      AA[p].atom[7].charge =  0.10; //    CD2
      AA[p].atom[8].charge =  0.30; //    CE1
      AA[p].atom[9].charge = -0.40; //    NE2
      AA[p].atom[10].charge = 0.25; //    H
      AA[p].atom[11].charge = 0.00; //    HE2
      AA[p].atom[12].charge = 0.00; //    HA
      AA[p].atom[13].charge = 0.00; //    2HB
      AA[p].atom[14].charge = 0.00; //    3HB
      AA[p].atom[15].charge = 0.00; //    HE1
      AA[p].atom[16].charge = 0.00; //    HD2
      break;

  }
  /*
   1.0 HSE-N:2H    0.33000 HC
   1.0 HSE-N:O    -0.51000 O
   1.0 HSE-N:HA    0.10000 HB
   1.0 HSE-N:3H    0.33000 HC
   1.0 HSE-N:CB   -0.08000 CT2
   1.0 HSE-N:HB1   0.09000 HA
   1.0 HSE-N:HB2   0.09000 HA
   1.0 HSE-N:ND1  -0.70000 NR2
   1.0 HSE-N:CG    0.22000 CPH1
   1.0 HSE-N:CE1   0.25000 CPH2
   1.0 HSE-N:HE1   0.13000 HR1
   1.0 HSE-N:NE2  -0.36000 NR1
   1.0 HSE-N:HE2   0.32000 H
   1.0 HSE-N:CD2  -0.05000 CPH1
   1.0 HSE-N:N    -0.30000 NH3
   1.0 HSE-N:HD2   0.09000 HR3
   1.0 HSE-N:1H    0.33000 HC
   1.0 HSE-N:C     0.51000 C
   1.0 HSE-N:CA    0.21000 CT1

   1.0 HSE-C:HA    0.09000 HB
   1.0 HSE-C:CB   -0.08000 CT2
   1.0 HSE-C:HB1   0.09000 HA
   1.0 HSE-C:HB2   0.09000 HA
   1.0 HSE-C:ND1  -0.70000 NR2
   1.0 HSE-C:CG    0.22000 CPH1
   1.0 HSE-C:CE1   0.25000 CPH2
   1.0 HSE-C:HE1   0.13000 HR1
   1.0 HSE-C:O    -0.67000 OC
   1.0 HSE-C:NE2  -0.36000 NR1
   1.0 HSE-C:OXT  -0.67000 OC
   1.0 HSE-C:HE2   0.32000 H
   1.0 HSE-C:CD2  -0.05000 CPH1
   1.0 HSE-C:N    -0.47000 NH1
   1.0 HSE-C:HD2   0.09000 HR3
   1.0 HSE-C:HN    0.31000 H
   1.0 HSE-C:C     0.34000 CC
   1.0 HSE-C:CA    0.07000 CT1
*/

  switch (opt) {
  	case EEF1:
	  AA[p].atom[0].fullatom_type =  12; // NH1    N
	  AA[p].atom[1].fullatom_type =   5; // CH1E   CA
	  AA[p].atom[2].fullatom_type =   4; // C      C
	  AA[p].atom[3].fullatom_type =  16; // O      O
	  AA[p].atom[4].fullatom_type =   6; // CH2E   CB
	  AA[p].atom[5].fullatom_type =  24; // CR     CG
	  AA[p].atom[6].fullatom_type =  10; // NR     ND1
	  AA[p].atom[7].fullatom_type =   8; // CR1E   CD2
	  AA[p].atom[8].fullatom_type =   8; // CR1E   CE1
	  AA[p].atom[9].fullatom_type =  12; // NH1    NE2
	  AA[p].atom[10].fullatom_type =  0; // H      H
	  AA[p].atom[11].fullatom_type =  0; // H      HE2
	  AA[p].atom[12].fullatom_type = -1; // Hapo   HA
	  AA[p].atom[13].fullatom_type = -1; // Hapo   2HB
	  AA[p].atom[14].fullatom_type = -1; // Hapo   3HB
	  AA[p].atom[15].fullatom_type = -1; // Hapo   HE1
	  AA[p].atom[16].fullatom_type = -1; // Hapo   HD2

	  AA[p].atom[17].fullatom_type = 17; // OC     OXT
	  AA[p].atom[18].fullatom_type =  1; // Hpol   2H
	  AA[p].atom[19].fullatom_type =  1; // Hpol   3H
	  break;
  	case Sybil:
  	  AA[p].atom[0].fullatom_type = 40; // Nbb    N
  	  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
  	  AA[p].atom[2].fullatom_type = 10; // CObb   C
  	  AA[p].atom[3].fullatom_type = 56; // OCbb   O
  	  AA[p].atom[4].fullatom_type = 8; // CH2    CB
  	  AA[p].atom[5].fullatom_type = 11; // aroC   CG
  	  AA[p].atom[6].fullatom_type = 44; // Nhis   ND1
  	  AA[p].atom[7].fullatom_type = 11; // aroC   CD2
  	  AA[p].atom[8].fullatom_type = 11; // aroC   CE1
  	  AA[p].atom[9].fullatom_type = 42; // Ntrp   NE2
  	  AA[p].atom[10].fullatom_type = 66; // HNbb   H
  	  AA[p].atom[11].fullatom_type = 66; // Hpol   HE2
  	  AA[p].atom[12].fullatom_type = 66; // Hapo   HA
  	  AA[p].atom[13].fullatom_type = 66; // Hapo  2HB
  	  AA[p].atom[14].fullatom_type = 66; // Hapo  3HB
  	  AA[p].atom[15].fullatom_type = 66; // Hapo   HE1
  	  AA[p].atom[16].fullatom_type = 66; // Hapo   HD2

  	  AA[p].atom[17].fullatom_type = 57; // OXT
  	  AA[p].atom[18].fullatom_type = 66; // Hpol
  	  AA[p].atom[19].fullatom_type = 66; // Hpol
  	  break;
  	default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 6; // aroC   CG
	  AA[p].atom[6].fullatom_type = 8; // Nhis   ND1
	  AA[p].atom[7].fullatom_type = 6; // aroC   CD2
	  AA[p].atom[8].fullatom_type = 6; // aroC   CE1
	  AA[p].atom[9].fullatom_type = 7; // Ntrp   NE2
	  AA[p].atom[10].fullatom_type = 25; // HNbb   H
	  AA[p].atom[11].fullatom_type = 22; // Hpol   HE2
	  AA[p].atom[12].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[13].fullatom_type = 23; // Hapo  2HB
	  AA[p].atom[14].fullatom_type = 23; // Hapo  3HB
	  AA[p].atom[15].fullatom_type = 23; // Hapo   HE1
	  AA[p].atom[16].fullatom_type = 23; // Hapo   HD2

	  AA[p].atom[17].fullatom_type = 20; // OXT
	  AA[p].atom[18].fullatom_type = 22; // Hpol
	  AA[p].atom[19].fullatom_type = 22; // Hpol
	  break;
  }

  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 10; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 12; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--CG
  AA[p].atom[4].bonded_neighbor[2] = 13; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 14; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 3; // CG
  AA[p].atom[5].bonded_neighbor[0] = 4; // CG--CB
  AA[p].atom[5].bonded_neighbor[1] = 6; // CG--ND1
  AA[p].atom[5].bonded_neighbor[2] = 7; // CG--CD2
  AA[p].atom[6].nbonded_neighbors = 2; // ND1
  AA[p].atom[6].bonded_neighbor[0] = 5; // ND1--CG
  AA[p].atom[6].bonded_neighbor[1] = 8; // ND1--CE1
  AA[p].atom[7].nbonded_neighbors = 3; // CD2
  AA[p].atom[7].bonded_neighbor[0] = 5; // CD2--CG
  AA[p].atom[7].bonded_neighbor[1] = 9; // CD2--NE2
  AA[p].atom[7].bonded_neighbor[2] = 16; // CD2--HD2
  AA[p].atom[8].nbonded_neighbors = 3; // CE1
  AA[p].atom[8].bonded_neighbor[0] = 6; // CE1--ND1
  AA[p].atom[8].bonded_neighbor[1] = 9; // CE1--NE2
  AA[p].atom[8].bonded_neighbor[2] = 15; // CE1--HE1
  AA[p].atom[9].nbonded_neighbors = 3; // NE2
  AA[p].atom[9].bonded_neighbor[0] = 7; // NE2--CD2
  AA[p].atom[9].bonded_neighbor[1] = 8; // NE2--CE1
  AA[p].atom[9].bonded_neighbor[2] = 11; // NE2--HE2
  AA[p].atom[10].nbonded_neighbors = 1; // H
  AA[p].atom[10].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[11].nbonded_neighbors = 1; // HE2
  AA[p].atom[11].bonded_neighbor[0] = 9; // HE2--NE2
  AA[p].atom[12].nbonded_neighbors = 1; // HA
  AA[p].atom[12].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[13].nbonded_neighbors = 1; //2HB
  AA[p].atom[13].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[14].nbonded_neighbors = 1; //3HB
  AA[p].atom[14].bonded_neighbor[0] = 4; //3HB--CB
  AA[p].atom[15].nbonded_neighbors = 1; // HE1
  AA[p].atom[15].bonded_neighbor[0] = 8; // HE1--CE1
  AA[p].atom[16].nbonded_neighbors = 1; // HD2
  AA[p].atom[16].bonded_neighbor[0] = 7; // HD2--CD2

  AA[p].atom[17].nbonded_neighbors = 1; // OXT
  AA[p].atom[17].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[18].nbonded_neighbors = 1; // 1H
  AA[p].atom[18].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[19].nbonded_neighbors = 1; // 2H
  AA[p].atom[19].bonded_neighbor[0] = 0; // H-N


  for ( k = 0; k < 11; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 11; k < 17; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template for building  HE2
  AA[p].atom[11].ta[0] = 9; //   NE2
  AA[p].atom[11].ta[1] = 7; //   CD2
  AA[p].atom[11].ta[2] = 8; //   CE1

  //bk   template for building  HA
  AA[p].atom[12].ta[0] = 1; //   CA
  AA[p].atom[12].ta[1] = 0; //   N
  AA[p].atom[12].ta[2] = 2; //   C

  //bk   template for building  2HB
  AA[p].atom[13].ta[0] = 4; //   CB
  AA[p].atom[13].ta[1] = 1; //   CA
  AA[p].atom[13].ta[2] = 5; //   CG

  //bk   template for building  3HB
  AA[p].atom[14].ta[0] = 4; //   CB
  AA[p].atom[14].ta[1] = 1; //   CA
  AA[p].atom[14].ta[2] = 5; //   CG

  //bk   template for building  HE1
  AA[p].atom[15].ta[0] = 8; //   CE1
  AA[p].atom[15].ta[1] = 6; //   ND1
  AA[p].atom[15].ta[2] = 9; //  NE2

  //bk   template for building  HD2
  AA[p].atom[16].ta[0] = 7; //   CD2
  AA[p].atom[16].ta[1] = 5; //   CG
  AA[p].atom[16].ta[2] = 9; //   NE2




  //bk   chi angles required to build atoms HIS
  //bk   chi angles needed for building  CG:HE1
  AA[p].chi_required[0] [5] = true;

  //bk   chi angles needed for building  ND1
  AA[p].chi_required[0] [6] = true;
  AA[p].chi_required[1] [6] = true;

  //bk   chi angles needed for building  CD2
  AA[p].chi_required[0] [7] = true;
  AA[p].chi_required[1] [7] = true;

  //bk   chi angles needed for building  CE1
  AA[p].chi_required[0] [8] = true;
  AA[p].chi_required[1] [8] = true;

  //bk   chi angles needed for building  NE2
  AA[p].chi_required[0] [9] = true;
  AA[p].chi_required[1] [9] = true;

  //bk   chi angles needed for building  HE2
  AA[p].chi_required[0] [11] = true;
  AA[p].chi_required[1] [11] = true;

  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [13] = true;

  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [14] = true;

  //bk   chi angles needed for building  HE1
  AA[p].chi_required[0] [15] = true;
  AA[p].chi_required[1] [15] = true;

  //bk   chi angles needed for building  HD2
  AA[p].chi_required[0] [16] = true;
  AA[p].chi_required[1] [16] = true;

  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   CG
  //bk   four atoms that define chi angle  2
  AA[p].chi_atoms[1] [0] = 1; //   CA
  AA[p].chi_atoms[1] [1] = 4; //   CB
  AA[p].chi_atoms[1] [2] = 5; //   CG:HE1
  AA[p].chi_atoms[1] [3] = 6; //   ND1



  AA[p].atom[0].icoor[0] = 0.0000; //   N
  AA[p].atom[0].icoor[1] = 12.0140; //   N
  AA[p].atom[0].icoor[2] = 22.5220; //   N
  AA[p].atom[1].icoor[0] = 0.0000; //   CA
  AA[p].atom[1].icoor[1] = 13.2300; //   CA
  AA[p].atom[1].icoor[2] = 23.3160; //   CA
  AA[p].atom[2].icoor[0] = 0.0000; //   C
  AA[p].atom[2].icoor[1] = 12.6360; //   C
  AA[p].atom[2].icoor[2] = 24.7310; //   C
  AA[p].atom[3].icoor[0] = 0.1540; //   O
  AA[p].atom[3].icoor[1] = 11.4260; //   O:HE1
  AA[p].atom[3].icoor[2] = 24.9130; //   O
  AA[p].atom[4].icoor[0] = -1.2080; //   CB
  AA[p].atom[4].icoor[1] = 14.1050; //   CB
  AA[p].atom[4].icoor[2] = 22.9660; //   CB
  AA[p].atom[5].icoor[0] = -1.2740; //   CG
  AA[p].atom[5].icoor[1] = 15.3800; //   CG
  AA[p].atom[5].icoor[2] = 23.7480; //   CG
  AA[p].atom[6].icoor[0] = -2.2840; //   ND1
  AA[p].atom[6].icoor[1] = 16.3040; //   ND1
  AA[p].atom[6].icoor[2] = 23.5800; //   ND1
  AA[p].atom[7].icoor[0] = -0.4560; //   CD2
  AA[p].atom[7].icoor[1] = 15.8850; //   CD2
  AA[p].atom[7].icoor[2] = 24.7010; //   CD2
  AA[p].atom[8].icoor[0] = -2.0840; //   CE1
  AA[p].atom[8].icoor[1] = 17.3230; //   CE1:HE1
  AA[p].atom[8].icoor[2] = 24.3980; //   CE1
  AA[p].atom[9].icoor[0] = -0.9820; //   NE2
  AA[p].atom[9].icoor[1] = 17.0930; //   NE2
  AA[p].atom[9].icoor[2] = 25.0880; //   NE2
  AA[p].atom[10].icoor[0] = -0.0410; //   H
  AA[p].atom[10].icoor[1] = 11.1270; //   H
  AA[p].atom[10].icoor[2] = 23.0030; //   H
  AA[p].atom[11].icoor[0] = -0.5830; //   HE2
  AA[p].atom[11].icoor[1] = 17.7010; //   HE2
  AA[p].atom[11].icoor[2] = 25.7890; //   HE2
  AA[p].atom[12].icoor[0] = 0.9100; //   HA
  AA[p].atom[12].icoor[1] = 13.7960; //   HA
  AA[p].atom[12].icoor[2] = 23.1190; //   HA
  AA[p].atom[13].icoor[0] = -1.1760; //  2HB
  AA[p].atom[13].icoor[1] = 14.3870; //  2HB:HE1
  AA[p].atom[13].icoor[2] = 21.9130; //  2HB
  AA[p].atom[14].icoor[0] = -2.1330; //  3HB
  AA[p].atom[14].icoor[1] = 13.5660; //  3HB
  AA[p].atom[14].icoor[2] = 23.1690; //  3HB
  AA[p].atom[15].icoor[0] = -2.7780; //   HE1
  AA[p].atom[15].icoor[1] = 18.1630; //   HE1
  AA[p].atom[15].icoor[2] = 24.4180; //   HE1
  AA[p].atom[16].icoor[0] = 0.4620; //   HD2
  AA[p].atom[16].icoor[1] = 15.5170; //   HD2
  AA[p].atom[16].icoor[2] = 25.1600; //   HD2

  // atom number for backbone HN
  AA[p].HNpos=10;
  // atom number for backbone HA
  AA[p].HApos=12;

  // number of polar hydrogens
  AA[p].nH_polar_complete=2;
  // atom numbers for polar H
  AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
  AA[p].Hpos_polar_complete[0]=10;
  AA[p].Hpos_polar_complete[1]=11;

  // number of aromatic hydrogens
  AA[p].nH_aromatic_complete=0;

  // number of apolar hydrogens
  AA[p].nH_apolar_complete=5;
  //atom number for apolar hydrogens
  AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
  AA[p].Hpos_apolar_complete[0]=12;
  AA[p].Hpos_apolar_complete[1]=13;
  AA[p].Hpos_apolar_complete[2]=14;
  AA[p].Hpos_apolar_complete[3]=15;
  AA[p].Hpos_apolar_complete[4]=16;

  // number of acceptors
  AA[p].nacceptors=2;
  //acceptor information
  AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
  (AA[p].aBase[1])[0]=6;(AA[p].aBase[1])[1]=5;
  (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;
  (AA[p].accpt_pos[1])[0]=1;(AA[p].accpt_pos[1])[1]=6;
  //OJO. HAY UNA RELACION ABASE2(6,8)

  //Number of Hydrogens connection
  AA[p].nH_hydrogen_connexions=7;
  //atoms hydrogens are connected too
  AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
  AA[p].Hbase[0][0]=10; AA[p].Hbase[0][1]=0;
  AA[p].Hbase[1][0]=11; AA[p].Hbase[1][1]=9;
  AA[p].Hbase[2][0]=12; AA[p].Hbase[2][1]=1;
  AA[p].Hbase[3][0]=13; AA[p].Hbase[3][1]=4;
  AA[p].Hbase[4][0]=14; AA[p].Hbase[4][1]=4;
  AA[p].Hbase[5][0]=15; AA[p].Hbase[5][1]=8;
  AA[p].Hbase[6][0]=16; AA[p].Hbase[6][1]=7;

  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=10;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=12;

  AA[p].atom[4].numHydrogens_atm=2;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=13;
  AA[p].atom[4].hydrogens_atm[1]=14;

  AA[p].atom[7].numHydrogens_atm=1;
  AA[p].atom[7].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[7].numHydrogens_atm);
  AA[p].atom[7].hydrogens_atm[0]=16;

  AA[p].atom[8].numHydrogens_atm=1;
  AA[p].atom[8].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[8].numHydrogens_atm);
  AA[p].atom[8].hydrogens_atm[0]=15;

  AA[p].atom[9].numHydrogens_atm=1;
  AA[p].atom[9].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[9].numHydrogens_atm);
  AA[p].atom[9].hydrogens_atm[0]=11;

  //FIN DE AMINOACIDO HISTIDINE
}

/**
* Initializes Isoleucine
*
* @param opt: Rosseta or ICM
*/
void init_ILE(Convention opt)
{
  //AMINOACIDO ISOLEUCINE
  int k, p =ILE;


  strcpy( AA[p].aa_name3, "ILE" );
  AA[p].aa_name1 = 'I';
  AA[p].mass = 113.161;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = false;
  AA[p].aa_is_nonpolar = true;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = false;

  AA[p].natoms = 19;
  AA[p].nheavyatoms = 8;
  AA[p].nchi = 2;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );

  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;

  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

  AA[p].chi_types = ( int * ) new int[AA[p].nchi];
  AA[p].chi_types[0] = 2;
  AA[p].chi_types[1] = 2;
  AA[p].natoms_EEF1 = 8;


  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " CG1" );
  strcpy( AA[p].atom[6].atom_name, " CG2" );
  strcpy( AA[p].atom[7].atom_name, " CD1" );
  strcpy( AA[p].atom[8].atom_name, " H  " );
  strcpy( AA[p].atom[9].atom_name, " HA " );
  strcpy( AA[p].atom[10].atom_name, " HB " );
  strcpy( AA[p].atom[11].atom_name, "1HG2" );
  strcpy( AA[p].atom[12].atom_name, "2HG2" );
  strcpy( AA[p].atom[13].atom_name, "3HG2" );
  strcpy( AA[p].atom[14].atom_name, "2HG1" );
  strcpy( AA[p].atom[15].atom_name, "3HG1" );
  strcpy( AA[p].atom[16].atom_name, "1HD1" );
  strcpy( AA[p].atom[17].atom_name, "2HD1" );
  strcpy( AA[p].atom[18].atom_name, "3HD1" );

  strcpy( AA[p].atom[19].atom_name," OXT" );
  strcpy( AA[p].atom[20].atom_name,"2H  " );
  strcpy( AA[p].atom[21].atom_name,"3H  " );

  switch(opt)
  {
    case Rosseta:
    case Sybil:
    AA[p].atom[0].charge = -0.47; //    N
      AA[p].atom[1].charge = 0.07; //    CA
      AA[p].atom[2].charge = 0.51; //    C
      AA[p].atom[3].charge = -0.51; //    O
      AA[p].atom[4].charge = -0.09; //    CB
      AA[p].atom[5].charge = -0.18; //    CG1
      AA[p].atom[6].charge = -0.27; //    CG2
      AA[p].atom[7].charge = -0.27; //    CD1
      AA[p].atom[8].charge = 0.31; //    H
      AA[p].atom[9].charge = 0.09; //    HA
      AA[p].atom[10].charge = 0.09; //    HB
      AA[p].atom[11].charge = 0.09; //    1HG2
      AA[p].atom[12].charge = 0.09; //    2HG2
      AA[p].atom[13].charge = 0.09; //    3HG2
      AA[p].atom[14].charge = 0.09; //    2HG1
      AA[p].atom[15].charge = 0.09; //    3HG1
      AA[p].atom[16].charge = 0.09; //    1HD1
      AA[p].atom[17].charge = 0.09; //    2HD1
      AA[p].atom[18].charge = 0.09; //    3HD1

      AA[p].atom[19].charge = -0.67; //   0XT
      AA[p].atom[20].charge = 0.33; //  1H
      AA[p].atom[21].charge = 0.33; //  2H
      break;

    case ICM:
      AA[p].atom[0].charge = -0.356; //    N
      AA[p].atom[1].charge = 0.064; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = -0.005; //    CB
      AA[p].atom[5].charge = -0.02; //    CG1
      AA[p].atom[6].charge = -0.075; //    CG2
      AA[p].atom[7].charge = -0.075; //    CD1
      AA[p].atom[8].charge = 0.176; //    H
      AA[p].atom[9].charge = 0.02; //    HA
      AA[p].atom[10].charge = 0.025; //    HB
      AA[p].atom[11].charge = 0.025; //    1HG2
      AA[p].atom[12].charge = 0.025; //    2HG2
      AA[p].atom[13].charge = 0.025; //    3HG2
      AA[p].atom[14].charge = 0.015; //    2HG1
      AA[p].atom[15].charge = 0.015; //    3HG1
      AA[p].atom[16].charge = 0.025; //    1HD1
      AA[p].atom[17].charge = 0.025; //    2HD1
      AA[p].atom[18].charge = 0.025; //    3HD1

      AA[p].atom[19].charge = -0.67; //   0XT
      AA[p].atom[20].charge = 0.33; //  1H
      AA[p].atom[21].charge = 0.33; //  2H
      break;
   case EEF1:
      AA[p].atom[0].charge = -0.35; //    N
      AA[p].atom[1].charge =  0.10; //    CA
      AA[p].atom[2].charge =  0.55; //    C
      AA[p].atom[3].charge = -0.55; //    O
      AA[p].atom[4].charge =  0.00; //    CB
      AA[p].atom[5].charge =  0.00; //    CG1
      AA[p].atom[6].charge =  0.00; //    CG2
      AA[p].atom[7].charge =  0.00; //    CD1
      AA[p].atom[8].charge =  0.25; //    H
      AA[p].atom[9].charge =  0.00; //    HA
      AA[p].atom[10].charge = 0.00; //    HB
      AA[p].atom[11].charge = 0.00; //    1HG2
      AA[p].atom[12].charge = 0.00; //    2HG2
      AA[p].atom[13].charge = 0.00; //    3HG2
      AA[p].atom[14].charge = 0.00; //    2HG1
      AA[p].atom[15].charge = 0.00; //    3HG1
      AA[p].atom[16].charge = 0.00; //    1HD1
      AA[p].atom[17].charge = 0.00; //    2HD1
      AA[p].atom[18].charge = 0.00; //    3HD1
      break;

  }
  /*
   1.0 ILE-N:2H    0.33000 HC
   1.0 ILE-N:HD3   0.09000 HA
   1.0 ILE-N:HA    0.10000 HB
   1.0 ILE-N:3H    0.33000 HC
   1.0 ILE-N:C     0.51000 C
   1.0 ILE-N:CB   -0.09000 CT1
   1.0 ILE-N:O    -0.51000 O
   1.0 ILE-N:HB    0.09000 HA
   1.0 ILE-N:CG2  -0.27000 CT3
   1.0 ILE-N:1HG2  0.09000 HA
   1.0 ILE-N:2HG2  0.09000 HA
   1.0 ILE-N:3HG2  0.09000 HA
   1.0 ILE-N:CG1  -0.18000 CT2
   1.0 ILE-N:HG11  0.09000 HA
   1.0 ILE-N:1HG1  0.09000 HA
   1.0 ILE-N:CD1  -0.27000 CT3
   1.0 ILE-N:N    -0.30000 NH3
   1.0 ILE-N:HD1   0.09000 HA
   1.0 ILE-N:1H    0.33000 HC
   1.0 ILE-N:HD2   0.09000 HA
   1.0 ILE-N:CA    0.21000 CT1

   1.0 ILE-C:HD3   0.09000 HA
   1.0 ILE-C:HA    0.09000 HB
   1.0 ILE-C:C     0.34000 CC
   1.0 ILE-C:CB   -0.09000 CT1
   1.0 ILE-C:HB    0.09000 HA
   1.0 ILE-C:CG2  -0.27000 CT3
   1.0 ILE-C:1HG2  0.09000 HA
   1.0 ILE-C:2HG2  0.09000 HA
   1.0 ILE-C:3HG2  0.09000 HA
   1.0 ILE-C:CG1  -0.18000 CT2
   1.0 ILE-C:O    -0.67000 OC
   1.0 ILE-C:HG11  0.09000 HA
   1.0 ILE-C:OXT  -0.67000 OC
   1.0 ILE-C:1HG1  0.09000 HA
   1.0 ILE-C:CD1  -0.27000 CT3
   1.0 ILE-C:N    -0.47000 NH1
   1.0 ILE-C:HD1   0.09000 HA
   1.0 ILE-C:HN    0.31000 H
   1.0 ILE-C:HD2   0.09000 HA
   1.0 ILE-C:CA    0.07000 CT1
 */

  switch (opt) {
  	case EEF1:
	  AA[p].atom[0].fullatom_type = 12; // NH1    N
	  AA[p].atom[1].fullatom_type =  5; // CH1E   CA
	  AA[p].atom[2].fullatom_type =  4; // C      C
	  AA[p].atom[3].fullatom_type = 16; // O      O
	  AA[p].atom[4].fullatom_type =  5; // CH1E   CB
	  AA[p].atom[5].fullatom_type =  6; // CH2E   CG1
	  AA[p].atom[6].fullatom_type =  7; // CH3E   CG2
	  AA[p].atom[7].fullatom_type =  7; // CH3E   CD1
	  AA[p].atom[8].fullatom_type =  0; // H      H
	  AA[p].atom[9].fullatom_type =  -1; // Hapo  HA
	  AA[p].atom[10].fullatom_type = -1; // Hapo  HB
	  AA[p].atom[11].fullatom_type = -1; // Hapo  1HG2
	  AA[p].atom[12].fullatom_type = -1; // Hapo  2HG2
	  AA[p].atom[13].fullatom_type = -1; // Hapo  3HG2
	  AA[p].atom[14].fullatom_type = -1; // Hapo  2HG1
	  AA[p].atom[15].fullatom_type = -1; // Hapo  3HG1
	  AA[p].atom[16].fullatom_type = -1; // Hapo  1HD1
	  AA[p].atom[17].fullatom_type = -1; // Hapo  2HD1
	  AA[p].atom[18].fullatom_type = -1; // Hapo  3HD1

	  AA[p].atom[19].fullatom_type = 17; // OC    OXT
	  AA[p].atom[20].fullatom_type =  1; // Hpol  2H
	  AA[p].atom[21].fullatom_type =  1; // Hpol  3H
	  break;
  	case Sybil:
	  AA[p].atom[0].fullatom_type = 40; // Nbb    N
	  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 10; // CObb   C
	  AA[p].atom[3].fullatom_type = 56; // OCbb   O
	  AA[p].atom[4].fullatom_type = 9; // CH1    CB
	  AA[p].atom[5].fullatom_type = 8; // CH2    CG1
	  AA[p].atom[6].fullatom_type = 6; // CH3    CG2
	  AA[p].atom[7].fullatom_type = 6; // CH3    CD1
	  AA[p].atom[8].fullatom_type = 66; // HNbb   H
	  AA[p].atom[9].fullatom_type = 66; // Hapo   HA
	  AA[p].atom[10].fullatom_type = 66; // Hapo   HB
	  AA[p].atom[11].fullatom_type = 66; // Hapo  1HG2
	  AA[p].atom[12].fullatom_type = 66; // Hapo  2HG2
	  AA[p].atom[13].fullatom_type = 66; // Hapo  3HG2
	  AA[p].atom[14].fullatom_type = 66; // Hapo  2HG1
	  AA[p].atom[15].fullatom_type = 66; // Hapo  3HG1
	  AA[p].atom[16].fullatom_type = 66; // Hapo  1HD1
	  AA[p].atom[17].fullatom_type = 66; // Hapo  2HD1
	  AA[p].atom[18].fullatom_type = 66; // Hapo  3HD1

	  AA[p].atom[19].fullatom_type = 57; // OXT
	  AA[p].atom[20].fullatom_type = 66; // Hpol
	  AA[p].atom[21].fullatom_type = 66; // Hpol
	  break;
  	default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 3; // CH1    CB
	  AA[p].atom[5].fullatom_type = 4; // CH2    CG1
	  AA[p].atom[6].fullatom_type = 5; // CH3    CG2
	  AA[p].atom[7].fullatom_type = 5; // CH3    CD1
	  AA[p].atom[8].fullatom_type = 25; // HNbb   H
	  AA[p].atom[9].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[10].fullatom_type = 23; // Hapo   HB
	  AA[p].atom[11].fullatom_type = 23; // Hapo  1HG2
	  AA[p].atom[12].fullatom_type = 23; // Hapo  2HG2
	  AA[p].atom[13].fullatom_type = 23; // Hapo  3HG2
	  AA[p].atom[14].fullatom_type = 23; // Hapo  2HG1
	  AA[p].atom[15].fullatom_type = 23; // Hapo  3HG1
	  AA[p].atom[16].fullatom_type = 23; // Hapo  1HD1
	  AA[p].atom[17].fullatom_type = 23; // Hapo  2HD1
	  AA[p].atom[18].fullatom_type = 23; // Hapo  3HD1

	  AA[p].atom[19].fullatom_type = 20; // OXT
	  AA[p].atom[20].fullatom_type = 22; // Hpol
	  AA[p].atom[21].fullatom_type = 22; // Hpol
	  break;
  }

  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 8; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 9; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--CG1
  AA[p].atom[4].bonded_neighbor[2] = 6; // CB--CG2
  AA[p].atom[4].bonded_neighbor[3] = 10; // CB--HB
  AA[p].atom[5].nbonded_neighbors = 4; // CG1
  AA[p].atom[5].bonded_neighbor[0] = 4; // CG1--CB
  AA[p].atom[5].bonded_neighbor[1] = 7; // CG1--CD1
  AA[p].atom[5].bonded_neighbor[2] = 14; // CG1--2HG1
  AA[p].atom[5].bonded_neighbor[3] = 15; // CG1--3HG1
  AA[p].atom[6].nbonded_neighbors = 4; // CG2
  AA[p].atom[6].bonded_neighbor[0] = 4; // CG2--CB
  AA[p].atom[6].bonded_neighbor[1] = 11; // CG2--1HG2
  AA[p].atom[6].bonded_neighbor[2] = 12; // CG2--2HG2
  AA[p].atom[6].bonded_neighbor[3] = 13; // CG2--3HG2
  AA[p].atom[7].nbonded_neighbors = 4; // CD1
  AA[p].atom[7].bonded_neighbor[0] = 5; // CD1--CG1
  AA[p].atom[7].bonded_neighbor[1] = 16; // CD1--1HD1
  AA[p].atom[7].bonded_neighbor[2] = 17; // CD1--2HD1
  AA[p].atom[7].bonded_neighbor[3] = 18; // CD1--3HD1
  AA[p].atom[8].nbonded_neighbors = 1; // H
  AA[p].atom[8].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[9].nbonded_neighbors = 1; // HA
  AA[p].atom[9].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[10].nbonded_neighbors = 1; // HB
  AA[p].atom[10].bonded_neighbor[0] = 4; // HB--CB
  AA[p].atom[11].nbonded_neighbors = 1; //1HG2
  AA[p].atom[11].bonded_neighbor[0] = 6; //1HG2--CG2
  AA[p].atom[12].nbonded_neighbors = 1; //2HG2
  AA[p].atom[12].bonded_neighbor[0] = 6; //2HG2--CG2
  AA[p].atom[13].nbonded_neighbors = 1; //3HG2
  AA[p].atom[13].bonded_neighbor[0] = 6; //3HG2--CG2
  AA[p].atom[14].nbonded_neighbors = 1; //2HG1
  AA[p].atom[14].bonded_neighbor[0] = 5; //2HG1--CG1
  AA[p].atom[15].nbonded_neighbors = 1; //3HG1
  AA[p].atom[15].bonded_neighbor[0] = 5; //3HG1--CG1
  AA[p].atom[16].nbonded_neighbors = 1; //1HD1
  AA[p].atom[16].bonded_neighbor[0] = 7; //1HD1--CD1
  AA[p].atom[17].nbonded_neighbors = 1; //2HD1
  AA[p].atom[17].bonded_neighbor[0] = 7; //2HD1--CD1
  AA[p].atom[18].nbonded_neighbors = 1; //3HD1
  AA[p].atom[18].bonded_neighbor[0] = 7; //3HD1--CD1

  AA[p].atom[19].nbonded_neighbors = 1; // OXT
  AA[p].atom[19].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[20].nbonded_neighbors = 1; // 1H
  AA[p].atom[20].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[21].nbonded_neighbors = 1; // 2H
  AA[p].atom[21].bonded_neighbor[0] = 0; // H-N


  for ( k = 0; k < 9; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 9; k < 19; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template for building  HA
  AA[p].atom[9].ta[0] = 1; //   CA
  AA[p].atom[9].ta[1] = 0; //   N
  AA[p].atom[9].ta[2] = 2; //   C

  //bk   template for building  HB
  AA[p].atom[10].ta[0] = 4; //   CB
  AA[p].atom[10].ta[1] = 1; //   CA
  AA[p].atom[10].ta[2] = 5; //   CG1

  //bk   template for building 1HG2
  AA[p].atom[11].ta[0] = 6; //   CG2
  AA[p].atom[11].ta[1] = 4; //   CB
  AA[p].atom[11].ta[2] = 1; //   CA

  //bk   template for building 2HG2
  AA[p].atom[12].ta[0] = 6; //   CG2
  AA[p].atom[12].ta[1] = 4; //   CB
  AA[p].atom[12].ta[2] = 1; //   CA

  //bk   template for building 3HG2
  AA[p].atom[13].ta[0] = 6; //   CG2
  AA[p].atom[13].ta[1] = 4; //   CB
  AA[p].atom[13].ta[2] = 1; //   CA

  //bk   template for building 2HG1
  AA[p].atom[14].ta[0] = 5; //   CG1
  AA[p].atom[14].ta[1] = 4; //   CB
  AA[p].atom[14].ta[2] = 7; //   CD1

  //bk   template for building 3HG1
  AA[p].atom[15].ta[0] = 5; //   CG1
  AA[p].atom[15].ta[1] = 4; //   CB
  AA[p].atom[15].ta[2] = 7; //   CD1

  //bk   template for building 1HD1
  AA[p].atom[16].ta[0] = 7; //   CD1
  AA[p].atom[16].ta[1] = 5; //   CG1
  AA[p].atom[16].ta[2] = 4; //   CB

  //bk   template for building 2HD1
  AA[p].atom[17].ta[0] = 7; //   CD1
  AA[p].atom[17].ta[1] = 5; //   CG1
  AA[p].atom[17].ta[2] = 4; //   CB

  //bk   template for building 3HD1
  AA[p].atom[18].ta[0] = 7; //   CD1
  AA[p].atom[18].ta[1] = 5; //   CG1
  AA[p].atom[18].ta[2] = 4; //   CB


  //bk   chi angles required to build atoms ILE
  //bk   chi angles needed for building  CG1
  AA[p].chi_required[0] [5] = true;

  //bk   chi angles needed for building  CG2
  AA[p].chi_required[0] [6] = true;

  //bk   chi angles needed for building  CD1
  AA[p].chi_required[0] [7] = true;
  AA[p].chi_required[1] [7] = true;

  //bk   chi angles needed for building  HB
  AA[p].chi_required[0] [10] = true;

  //bk   chi angles needed for building 1HG2
  AA[p].chi_required[0] [11] = true;

  //bk   chi angles needed for building 2HG2
  AA[p].chi_required[0] [12] = true;

  //bk   chi angles needed for building 3HG2
  AA[p].chi_required[0] [13] = true;

  //bk   chi angles needed for building 2HG1
  AA[p].chi_required[0] [14] = true;
  AA[p].chi_required[1] [14] = true;

  //bk   chi angles needed for building 3HG1
  AA[p].chi_required[0] [15] = true;
  AA[p].chi_required[1] [15] = true;

  //bk   chi angles needed for building 1HD1
  AA[p].chi_required[0] [16] = true;
  AA[p].chi_required[1] [16] = true;

  //bk   chi angles needed for building 2HD1
  AA[p].chi_required[0] [17] = true;
  AA[p].chi_required[1] [17] = true;

  //bk   chi angles needed for building 3HD1
  AA[p].chi_required[0] [18] = true;
  AA[p].chi_required[1] [18] = true;




  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   CG1
  //bk   four atoms that define chi angle  2
  AA[p].chi_atoms[1] [0] = 1; //   CA
  AA[p].chi_atoms[1] [1] = 4; //   CB
  AA[p].chi_atoms[1] [2] = 5; //   CG1
  AA[p].chi_atoms[1] [3] = 7; //   CD1



  AA[p].atom[0].icoor[0] = 0.0000; //   N
  AA[p].atom[0].icoor[1] = 13.5100; //   N
  AA[p].atom[0].icoor[2] = 25.7330; //   N
  AA[p].atom[1].icoor[0] = 0.0000; //   CA
  AA[p].atom[1].icoor[1] = 13.0790; //   CA
  AA[p].atom[1].icoor[2] = 27.1250; //   CA
  AA[p].atom[2].icoor[0] = 0.0000; //   C
  AA[p].atom[2].icoor[1] = 14.2670; //   C
  AA[p].atom[2].icoor[2] = 28.0870; //   C
  AA[p].atom[3].icoor[0] = 0.0680; //   O
  AA[p].atom[3].icoor[1] = 15.4230; //   O
  AA[p].atom[3].icoor[2] = 27.6610; //   O
  AA[p].atom[4].icoor[0] = 1.2140; //   CB
  AA[p].atom[4].icoor[1] = 12.1860; //   CB
  AA[p].atom[4].icoor[2] = 27.4400; //   CB
  AA[p].atom[5].icoor[0] = 1.1890; //   CG1
  AA[p].atom[5].icoor[1] = 11.7490; //   CG1
  AA[p].atom[5].icoor[2] = 28.9070; //   CG1
  AA[p].atom[6].icoor[0] = 2.5090; //   CG2
  AA[p].atom[6].icoor[1] = 12.9170; //   CG2
  AA[p].atom[6].icoor[2] = 27.1210; //   CG2
  AA[p].atom[7].icoor[0] = 2.3510; //   CD1
  AA[p].atom[7].icoor[1] = 10.8660; //   CD1
  AA[p].atom[7].icoor[2] = 29.3010; //   CD1
  AA[p].atom[8].icoor[0] = 0.0310; //   H
  AA[p].atom[8].icoor[1] = 14.5010; //   H
  AA[p].atom[8].icoor[2] = 25.5390; //   H
  AA[p].atom[9].icoor[0] = -0.9170; //   HA
  AA[p].atom[9].icoor[1] = 12.5400; //   HA
  AA[p].atom[9].icoor[2] = 27.3600; //   HA
  AA[p].atom[10].icoor[0] = 1.1520; //   HB
  AA[p].atom[10].icoor[1] = 11.2780; //   HB
  AA[p].atom[10].icoor[2] = 26.8410; //   HB
  AA[p].atom[11].icoor[0] = 3.3570; //  1HG2
  AA[p].atom[11].icoor[1] = 12.2720; //  1HG2
  AA[p].atom[11].icoor[2] = 27.3500; //  1HG2
  AA[p].atom[12].icoor[0] = 2.5280; //  2HG2
  AA[p].atom[12].icoor[1] = 13.1790; //  2HG2
  AA[p].atom[12].icoor[2] = 26.0640; //  2HG2
  AA[p].atom[13].icoor[0] = 2.5720; //  3HG2
  AA[p].atom[13].icoor[1] = 13.8250; //  3HG2
  AA[p].atom[13].icoor[2] = 27.7210; //  3HG2
  AA[p].atom[14].icoor[0] = 1.1970; //  2HG1
  AA[p].atom[14].icoor[1] = 12.6530; //  2HG1
  AA[p].atom[14].icoor[2] = 29.5150; //  2HG1
  AA[p].atom[15].icoor[0] = 0.2540; //  3HG1
  AA[p].atom[15].icoor[1] = 11.2120; //  3HG1
  AA[p].atom[15].icoor[2] = 29.0690; //  3HG1
  AA[p].atom[16].icoor[0] = 2.2640; //  1HD1
  AA[p].atom[16].icoor[1] = 10.5970; //  1HD1
  AA[p].atom[16].icoor[2] = 30.3540; //  1HD1
  AA[p].atom[17].icoor[0] = 2.3430; //  2HD1
  AA[p].atom[17].icoor[1] = 9.9600; //  2HD1
  AA[p].atom[17].icoor[2] = 28.6940; //  2HD1
  AA[p].atom[18].icoor[0] = 3.2860; //  3HD1
  AA[p].atom[18].icoor[1] = 11.4010; //  3HD1
  AA[p].atom[18].icoor[2] = 29.1410; //  3HD1

  // atom number for backbone HN
  AA[p].HNpos=8;
  // atom number for backbone HA
  AA[p].HApos=9;

  // number of polar hydrogens
  AA[p].nH_polar_complete=1;
  // atom numbers for polar H
  AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
  AA[p].Hpos_polar_complete[0]=8;

  // number of aromatic hydrogens
  AA[p].nH_aromatic_complete=0;

  // number of apolar hydrogens
  AA[p].nH_apolar_complete=10;
  //atom number for apolar hydrogens
  AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
  AA[p].Hpos_apolar_complete[0]=9;
  AA[p].Hpos_apolar_complete[1]=10;
  AA[p].Hpos_apolar_complete[2]=11;
  AA[p].Hpos_apolar_complete[3]=12;
  AA[p].Hpos_apolar_complete[4]=13;
  AA[p].Hpos_apolar_complete[5]=14;
  AA[p].Hpos_apolar_complete[6]=15;
  AA[p].Hpos_apolar_complete[7]=16;
  AA[p].Hpos_apolar_complete[8]=17;
  AA[p].Hpos_apolar_complete[9]=18;

  // number of acceptors
  AA[p].nacceptors=1;
//acceptor information
AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
(AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
(AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;


  //Number of Hydrogens connection
  AA[p].nH_hydrogen_connexions=11;
  //atoms hydrogens are connected too
  AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
  AA[p].Hbase[0][0]=10; AA[p].Hbase[0][1]=0;
  AA[p].Hbase[1][0]=11; AA[p].Hbase[1][1]=9;
  AA[p].Hbase[2][0]=12; AA[p].Hbase[2][1]=1;
  AA[p].Hbase[3][0]=13; AA[p].Hbase[3][1]=4;
  AA[p].Hbase[4][0]=14; AA[p].Hbase[4][1]=4;
  AA[p].Hbase[5][0]=15; AA[p].Hbase[5][1]=8;
  AA[p].Hbase[6][0]=16; AA[p].Hbase[6][1]=7;
  AA[p].Hbase[7][0]=13; AA[p].Hbase[7][1]=4;
  AA[p].Hbase[8][0]=14; AA[p].Hbase[8][1]=4;
  AA[p].Hbase[9][0]=15; AA[p].Hbase[9][1]=8;
  AA[p].Hbase[10][0]=16; AA[p].Hbase[10][1]=7;

  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=8;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=9;

  AA[p].atom[4].numHydrogens_atm=1;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=10;

  AA[p].atom[5].numHydrogens_atm=2;
  AA[p].atom[5].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[5].numHydrogens_atm);
  AA[p].atom[5].hydrogens_atm[0]=14;
  AA[p].atom[5].hydrogens_atm[1]=15;

  AA[p].atom[6].numHydrogens_atm=3;
  AA[p].atom[6].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[6].numHydrogens_atm);
  AA[p].atom[6].hydrogens_atm[0]=11;
  AA[p].atom[6].hydrogens_atm[1]=12;
  AA[p].atom[6].hydrogens_atm[2]=13;

  AA[p].atom[7].numHydrogens_atm=3;
  AA[p].atom[7].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[7].numHydrogens_atm);
  AA[p].atom[7].hydrogens_atm[0]=16;
  AA[p].atom[7].hydrogens_atm[1]=17;
  AA[p].atom[7].hydrogens_atm[2]=18;

}

//FIN DE AMINOACIDO  ISOLEUCINE

/**
* Initializes Lysine
*
* @param opt: Rosseta or ICM
*/
void init_LYS(Convention opt)
{
  //AMINOACIDO LYSINA
  int k, p =LYS;


  strcpy( AA[p].aa_name3, "LYS" );
  AA[p].aa_name1 = 'K';
  AA[p].mass = 129.184;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = true;
  AA[p].aa_is_nonpolar = false;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = true;

  AA[p].natoms = 22;
  AA[p].nheavyatoms = 9;
  AA[p].nchi = 4;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );
  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;


  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

  AA[p].chi_types = ( int * ) new int[AA[p].nchi];
  AA[p].chi_types[0] = 2;
  AA[p].chi_types[1] = 2;
  AA[p].chi_types[2] = 2;
  AA[p].chi_types[3] = 2;
  AA[p].natoms_EEF1 = 12;


  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " CG " );
  strcpy( AA[p].atom[6].atom_name, " CD " );
  strcpy( AA[p].atom[7].atom_name, " CE " );
  strcpy( AA[p].atom[8].atom_name, " NZ " );
  strcpy( AA[p].atom[9].atom_name, " H  " );
  strcpy( AA[p].atom[10].atom_name, "1HZ " );
  strcpy( AA[p].atom[11].atom_name, "2HZ " );
  strcpy( AA[p].atom[12].atom_name, "3HZ " );
  strcpy( AA[p].atom[13].atom_name, " HA " );
  strcpy( AA[p].atom[14].atom_name, "2HB " );
  strcpy( AA[p].atom[15].atom_name, "3HB " );
  strcpy( AA[p].atom[16].atom_name, "2HG " );
  strcpy( AA[p].atom[17].atom_name, "3HG " );
  strcpy( AA[p].atom[18].atom_name, "2HD " );
  strcpy( AA[p].atom[19].atom_name, "3HD " );
  strcpy( AA[p].atom[20].atom_name, "2HE " );
  strcpy( AA[p].atom[21].atom_name, "3HE " );

  strcpy( AA[p].atom[22].atom_name," OXT" );
  strcpy( AA[p].atom[23].atom_name,"2H  " );
  strcpy( AA[p].atom[24].atom_name,"3H  " );

  switch(opt)
  {
    case Rosseta:
    case Sybil:
     AA[p].atom[0].charge = -0.47; //    N
      AA[p].atom[1].charge = 0.07; //    CA
      AA[p].atom[2].charge = 0.51; //    C
      AA[p].atom[3].charge = -0.51; //    O
      AA[p].atom[4].charge = -0.18; //    CB
      AA[p].atom[5].charge = -0.18; //    CG
      AA[p].atom[6].charge = -0.18; //    CD
      AA[p].atom[7].charge = 0.21; //    CE
      AA[p].atom[8].charge = -0.30; //    NZ
      AA[p].atom[9].charge = 0.31; //    H
      AA[p].atom[10].charge = 0.33; //    1HZ
      AA[p].atom[11].charge = 0.33; //    2HZ
      AA[p].atom[12].charge = 0.33; //    3HZ
      AA[p].atom[13].charge = 0.09; //    HA
      AA[p].atom[14].charge = 0.09; //    2HB
      AA[p].atom[15].charge = 0.09; //    3HB
      AA[p].atom[16].charge = 0.09; //    2HG
      AA[p].atom[17].charge = 0.09; //    3HG
      AA[p].atom[18].charge = 0.09; //    2HD
      AA[p].atom[19].charge = 0.09; //    3HD
      AA[p].atom[20].charge = 0.05; //    2HE
      AA[p].atom[21].charge = 0.05; //    3HE

      AA[p].atom[22].charge = -0.67; //   0XT
      AA[p].atom[23].charge = 0.33; //  1H
      AA[p].atom[24].charge = 0.33; //  2H
      break;

    case ICM:
      AA[p].atom[0].charge = -0.356; //    N
      AA[p].atom[1].charge = 0.064; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = -0.03; //    CB
      AA[p].atom[5].charge = -0.025; //    CG
      AA[p].atom[6].charge = -0.115; //    CD
      AA[p].atom[7].charge = 0.05; //    CE
      AA[p].atom[8].charge = -0.32; //    NZ
      AA[p].atom[9].charge = 0.176; //    H
      AA[p].atom[10].charge = 0.32; //    1HZ
      AA[p].atom[11].charge = 0.32; //    2HZ
      AA[p].atom[12].charge = 0.32; //    3HZ
      AA[p].atom[13].charge = 0.02; //    HA
      AA[p].atom[14].charge = 0.015; //    2HB
      AA[p].atom[15].charge = 0.015; //    3HB
      AA[p].atom[16].charge = 0.02; //    2HG
      AA[p].atom[17].charge = 0.02; //    3HG
      AA[p].atom[18].charge = 0.1; //    2HD
      AA[p].atom[19].charge = 0.1; //    3HD
      AA[p].atom[20].charge = 0.12; //    2HE
      AA[p].atom[21].charge = 0.12; //    3HE

      AA[p].atom[22].charge = -0.67; //   0XT
      AA[p].atom[23].charge = 0.33; //  1H
      AA[p].atom[24].charge = 0.33; //  2H
      break;
     case EEF1:
      AA[p].atom[0].charge = -0.35; //    N
      AA[p].atom[1].charge =  0.10; //    CA
      AA[p].atom[2].charge =  0.55; //    C
      AA[p].atom[3].charge = -0.55; //    O
      AA[p].atom[4].charge =  0.00; //    CB
      AA[p].atom[5].charge =  0.00; //    CG
      AA[p].atom[6].charge =  0.00; //    CD
      AA[p].atom[7].charge =  0.00; //    CE
      AA[p].atom[8].charge = -1.35; //    NZ
      AA[p].atom[9].charge =  0.25; //    H
      AA[p].atom[10].charge = 0.45; //    1HZ
      AA[p].atom[11].charge = 0.45; //    2HZ
      AA[p].atom[12].charge = 0.45; //    3HZ
      AA[p].atom[13].charge = 0.00; //    HA
      AA[p].atom[14].charge = 0.00; //    2HB
      AA[p].atom[15].charge = 0.00; //    3HB
      AA[p].atom[16].charge = 0.00; //    2HG
      AA[p].atom[17].charge = 0.00; //    3HG
      AA[p].atom[18].charge = 0.00; //    2HD
      AA[p].atom[19].charge = 0.00; //    3HD
      AA[p].atom[20].charge = 0.00; //    2HE
      AA[p].atom[21].charge = 0.00; //    3HE
      break;
  }


  /*
   1.0 LYS-N:2H    0.33000 HC
   1.0 LYS-N:NZ   -0.30000 NH3
   1.0 LYS-N:HA    0.10000 HB
   1.0 LYS-N:3H    0.33000 HC
   1.0 LYS-N:1HZ   0.33000 HC
   1.0 LYS-N:CB   -0.18000 CT2
   1.0 LYS-N:2HZ   0.33000 HC
   1.0 LYS-N:HB1   0.09000 HA
   1.0 LYS-N:1HB   0.09000 HA
   1.0 LYS-N:CG   -0.18000 CT2
   1.0 LYS-N:HG1   0.09000 HA
   1.0 LYS-N:1HG   0.09000 HA
   1.0 LYS-N:3HZ   0.33000 HC
   1.0 LYS-N:CD   -0.18000 CT2
   1.0 LYS-N:C     0.51000 C
   1.0 LYS-N:HD1   0.09000 HA
   1.0 LYS-N:O    -0.51000 O
   1.0 LYS-N:1HD   0.09000 HA
   1.0 LYS-N:CE    0.21000 CT2
   1.0 LYS-N:N    -0.30000 NH3
   1.0 LYS-N:HE1   0.05000 HA
   1.0 LYS-N:1H    0.33000 HC
   1.0 LYS-N:1HE   0.05000 HA
   1.0 LYS-N:CA    0.21000 CT1

   1.0 LYS-C:NZ   -0.30000 NH3
   1.0 LYS-C:HA    0.09000 HB
   1.0 LYS-C:1HZ   0.33000 HC
   1.0 LYS-C:CB   -0.18000 CT2
   1.0 LYS-C:2HZ   0.33000 HC
   1.0 LYS-C:HB1   0.09000 HA
   1.0 LYS-C:1HB   0.09000 HA
   1.0 LYS-C:CG   -0.18000 CT2
   1.0 LYS-C:HG1   0.09000 HA
   1.0 LYS-C:1HG   0.09000 HA
   1.0 LYS-C:3HZ   0.33000 HC
   1.0 LYS-C:CD   -0.18000 CT2
   1.0 LYS-C:O    -0.67000 OC
   1.0 LYS-C:C     0.34000 CC
   1.0 LYS-C:HD1   0.09000 HA
   1.0 LYS-C:OXT  -0.67000 OC
   1.0 LYS-C:1HD   0.09000 HA
   1.0 LYS-C:CE    0.21000 CT2
   1.0 LYS-C:N    -0.47000 NH1
   1.0 LYS-C:HE1   0.05000 HA
   1.0 LYS-C:HN    0.31000 H
   1.0 LYS-C:1HE   0.05000 HA
   1.0 LYS-C:CA    0.07000 CT1
*/

  switch (opt) {
  	case EEF1:
	  AA[p].atom[0].fullatom_type =  12; // NH1    N
	  AA[p].atom[1].fullatom_type =   5; // CH1E   CA
	  AA[p].atom[2].fullatom_type =   4; // C      C
	  AA[p].atom[3].fullatom_type =  16; // O      O
	  AA[p].atom[4].fullatom_type =   6; // CH2E   CB
	  AA[p].atom[5].fullatom_type =   6; // CH2E   CG
	  AA[p].atom[6].fullatom_type =   6; // CH2E   CD
	  AA[p].atom[7].fullatom_type =   6; // CH2E   CE
	  AA[p].atom[8].fullatom_type =  14; // NH3    NZ
	  AA[p].atom[9].fullatom_type =   0; // H      H
	  AA[p].atom[10].fullatom_type =  1; // HC     1HZ
	  AA[p].atom[11].fullatom_type =  1; // HC     2HZ
	  AA[p].atom[12].fullatom_type =  1; // HC     3HZ
	  AA[p].atom[13].fullatom_type = -1; // Hapo   HA
	  AA[p].atom[14].fullatom_type = -1; // Hapo   2HB
	  AA[p].atom[15].fullatom_type = -1; // Hapo   3HB
	  AA[p].atom[16].fullatom_type = -1; // Hapo   2HG
	  AA[p].atom[17].fullatom_type = -1; // Hapo   3HG
	  AA[p].atom[18].fullatom_type = -1; // Hapo   2HD
	  AA[p].atom[19].fullatom_type = -1; // Hapo   3HD
	  AA[p].atom[20].fullatom_type = -1; // Hapo   2HE
	  AA[p].atom[21].fullatom_type = -1; // Hapo   3HE

	  AA[p].atom[22].fullatom_type = 17; // OC     OXT
	  AA[p].atom[23].fullatom_type =  1; // HC     2H
	  AA[p].atom[24].fullatom_type =  1; // HC     3H
	  break;
  	case Sybil:
	  AA[p].atom[0].fullatom_type = 40; // Nbb    N
	  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 10; // CObb   C
	  AA[p].atom[3].fullatom_type = 56; // OCbb   O
	  AA[p].atom[4].fullatom_type = 8; // CH2    CB
	  AA[p].atom[5].fullatom_type = 8; // CH2    CG
	  AA[p].atom[6].fullatom_type = 8; // CH2    CD
	  AA[p].atom[7].fullatom_type = 5; // CH2    CE
	  AA[p].atom[8].fullatom_type = 35; // Nlys   NZ
	  AA[p].atom[9].fullatom_type = 66; // HNbb   H
	  AA[p].atom[10].fullatom_type = 66; // Hpol  1HZ
	  AA[p].atom[11].fullatom_type = 66; // Hpol  2HZ
	  AA[p].atom[12].fullatom_type = 66; // Hpol  3HZ
	  AA[p].atom[13].fullatom_type = 66; // Hapo   HA
	  AA[p].atom[14].fullatom_type = 66; // Hapo  2HB
	  AA[p].atom[15].fullatom_type = 66; // Hapo  3HB
	  AA[p].atom[16].fullatom_type = 66; // Hapo  2HG
	  AA[p].atom[17].fullatom_type = 66; // Hapo  3HG
	  AA[p].atom[18].fullatom_type = 66; // Hapo  2HD
	  AA[p].atom[19].fullatom_type = 66; // Hapo  3HD
	  AA[p].atom[20].fullatom_type = 66; // Hapo  2HE
	  AA[p].atom[21].fullatom_type = 66; // Hapo  3HE

	  AA[p].atom[22].fullatom_type = 57; // OXT
	  AA[p].atom[23].fullatom_type = 66; // Hpol
	  AA[p].atom[24].fullatom_type = 66; // Hpol
	  break;
  	default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 4; // CH2    CG
	  AA[p].atom[6].fullatom_type = 4; // CH2    CD
	  AA[p].atom[7].fullatom_type = 4; // CH2    CE
	  AA[p].atom[8].fullatom_type = 10; // Nlys   NZ
	  AA[p].atom[9].fullatom_type = 25; // HNbb   H
	  AA[p].atom[10].fullatom_type = 22; // Hpol  1HZ
	  AA[p].atom[11].fullatom_type = 22; // Hpol  2HZ
	  AA[p].atom[12].fullatom_type = 22; // Hpol  3HZ
	  AA[p].atom[13].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[14].fullatom_type = 23; // Hapo  2HB
	  AA[p].atom[15].fullatom_type = 23; // Hapo  3HB
	  AA[p].atom[16].fullatom_type = 23; // Hapo  2HG
	  AA[p].atom[17].fullatom_type = 23; // Hapo  3HG
	  AA[p].atom[18].fullatom_type = 23; // Hapo  2HD
	  AA[p].atom[19].fullatom_type = 23; // Hapo  3HD
	  AA[p].atom[20].fullatom_type = 23; // Hapo  2HE
	  AA[p].atom[21].fullatom_type = 23; // Hapo  3HE

	  AA[p].atom[22].fullatom_type = 20; // OXT
	  AA[p].atom[23].fullatom_type = 22; // Hpol
	  AA[p].atom[24].fullatom_type = 22; // Hpol
	  break;
  }

  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 9; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 13; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--CG
  AA[p].atom[4].bonded_neighbor[2] = 14; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 15; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 4; // CG
  AA[p].atom[5].bonded_neighbor[0] = 4; // CG--CB
  AA[p].atom[5].bonded_neighbor[1] = 6; // CG--CD
  AA[p].atom[5].bonded_neighbor[2] = 16; // CG--2HG
  AA[p].atom[5].bonded_neighbor[3] = 17; // CG--3HG
  AA[p].atom[6].nbonded_neighbors = 4; // CD
  AA[p].atom[6].bonded_neighbor[0] = 5; // CD--CG
  AA[p].atom[6].bonded_neighbor[1] = 7; // CD--CE
  AA[p].atom[6].bonded_neighbor[2] = 18; // CD--2HD
  AA[p].atom[6].bonded_neighbor[3] = 19; // CD--3HD
  AA[p].atom[7].nbonded_neighbors = 4; // CE
  AA[p].atom[7].bonded_neighbor[0] = 6; // CE--CD
  AA[p].atom[7].bonded_neighbor[1] = 8; // CE--NZ
  AA[p].atom[7].bonded_neighbor[2] = 20; // CE--2HE
  AA[p].atom[7].bonded_neighbor[3] = 21; // CE--3HE
  AA[p].atom[8].nbonded_neighbors = 4; // NZ
  AA[p].atom[8].bonded_neighbor[0] = 7; // NZ--CE
  AA[p].atom[8].bonded_neighbor[1] = 10; // NZ--1HZ
  AA[p].atom[8].bonded_neighbor[2] = 11; // NZ--2HZ
  AA[p].atom[8].bonded_neighbor[3] = 12; // NZ--3HZ
  AA[p].atom[9].nbonded_neighbors = 1; // H
  AA[p].atom[9].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[10].nbonded_neighbors = 1; //1HZ
  AA[p].atom[10].bonded_neighbor[0] = 8; //1HZ--NZ
  AA[p].atom[11].nbonded_neighbors = 1; //2HZ
  AA[p].atom[11].bonded_neighbor[0] = 8; //2HZ--NZ
  AA[p].atom[12].nbonded_neighbors = 1; //3HZ
  AA[p].atom[12].bonded_neighbor[0] = 8; //3HZ--NZ
  AA[p].atom[13].nbonded_neighbors = 1; // HA
  AA[p].atom[13].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[14].nbonded_neighbors = 1; //2HB
  AA[p].atom[14].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[15].nbonded_neighbors = 1; //3HB
  AA[p].atom[15].bonded_neighbor[0] = 4; //3HB--CB
  AA[p].atom[16].nbonded_neighbors = 1; //2HG
  AA[p].atom[16].bonded_neighbor[0] = 5; //2HG--CG
  AA[p].atom[17].nbonded_neighbors = 1; //3HG
  AA[p].atom[17].bonded_neighbor[0] = 5; //3HG--CG
  AA[p].atom[18].nbonded_neighbors = 1; //2HD
  AA[p].atom[18].bonded_neighbor[0] = 6; //2HD--CD
  AA[p].atom[19].nbonded_neighbors = 1; //3HD
  AA[p].atom[19].bonded_neighbor[0] = 6; //3HD--CD
  AA[p].atom[20].nbonded_neighbors = 1; //2HE
  AA[p].atom[20].bonded_neighbor[0] = 7; //2HE--CE
  AA[p].atom[21].nbonded_neighbors = 1; //3HE
  AA[p].atom[21].bonded_neighbor[0] = 7; //3HE--CE

  AA[p].atom[22].nbonded_neighbors = 1; // OXT
  AA[p].atom[22].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[23].nbonded_neighbors = 1; // 1H
  AA[p].atom[23].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[24].nbonded_neighbors = 1; // 2H
  AA[p].atom[24].bonded_neighbor[0] = 0; // H-N


  for ( k = 0; k < 10; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 10; k < 22; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template for building 1HZ
  AA[p].atom[10].ta[0] = 8; //   NZ
  AA[p].atom[10].ta[1] = 7; //   CE
  AA[p].atom[10].ta[2] = 6; //   CD

  //bk   template for building 2HZ
  AA[p].atom[11].ta[0] = 8; //   NZ
  AA[p].atom[11].ta[1] = 7; //   CE
  AA[p].atom[11].ta[2] = 6; //   CD

  //bk   template for building 3HZ
  AA[p].atom[12].ta[0] = 8; //   NZ
  AA[p].atom[12].ta[1] = 7; //   CE
  AA[p].atom[12].ta[2] = 6; //   CD

  //bk   template for building  HA
  AA[p].atom[13].ta[0] = 1; //   CA
  AA[p].atom[13].ta[1] = 0; //   N
  AA[p].atom[13].ta[2] = 2; //   C

  //bk   template for building 2HB
  AA[p].atom[14].ta[0] = 4; //   CB
  AA[p].atom[14].ta[1] = 1; //   CA
  AA[p].atom[14].ta[2] = 5; //   CG

  //bk   template for building 3HB
  AA[p].atom[15].ta[0] = 4; //   CB
  AA[p].atom[15].ta[1] = 1; //   CA
  AA[p].atom[15].ta[2] = 5; //   CG

  //bk   template for building 2HG
  AA[p].atom[16].ta[0] = 5; //   CG
  AA[p].atom[16].ta[1] = 4; //   CB
  AA[p].atom[16].ta[2] = 6; //   CD

  //bk   template for building 3HG
  AA[p].atom[17].ta[0] = 5; //   CG
  AA[p].atom[17].ta[1] = 4; //   CB
  AA[p].atom[17].ta[2] = 6; //   CD

  //bk   template for building 2HD
  AA[p].atom[18].ta[0] = 6; //   CD
  AA[p].atom[18].ta[1] = 5; //   CG
  AA[p].atom[18].ta[2] = 7; //   CE

  //bk   template for building 3HD
  AA[p].atom[19].ta[0] = 6; //   CD
  AA[p].atom[19].ta[1] = 5; //   CG
  AA[p].atom[19].ta[2] = 7; //   CE

  //bk   template for building 2HE
  AA[p].atom[20].ta[0] = 7; //   CE
  AA[p].atom[20].ta[1] = 6; //   CD
  AA[p].atom[20].ta[2] = 8; //   NZ

  //bk   template for building 3HE
  AA[p].atom[21].ta[0] = 7; //   CE
  AA[p].atom[21].ta[1] = 6; //   CD
  AA[p].atom[21].ta[2] = 8; //   NZ


  //bk   chi angles required to build atoms LYS

  //bk   chi angles needed for building  N
  //bk   chi angles needed for building  CA
  //bk   chi angles needed for building  C
  //bk   chi angles needed for building  O
  //bk   chi angles needed for building  CB
  //bk   chi angles needed for building  CG
  AA[p].chi_required[0] [5] = true;


  //bk   chi angles needed for building  CD
  AA[p].chi_required[0] [6] = true;
  AA[p].chi_required[1] [6] = true;

  //bk   chi angles needed for building  CE
  AA[p].chi_required[0] [7] = true;
  AA[p].chi_required[1] [7] = true;
  AA[p].chi_required[2] [7] = true;

  //bk   chi angles needed for building  NZ
  AA[p].chi_required[0] [8] = true;
  AA[p].chi_required[1] [8] = true;
  AA[p].chi_required[2] [8] = true;
  AA[p].chi_required[3] [8] = true;

  //bk   chi angles needed for building  H

  //bk   chi angles needed for building 1HZ
  AA[p].chi_required[0] [10] = true;
  AA[p].chi_required[1] [10] = true;
  AA[p].chi_required[2] [10] = true;
  AA[p].chi_required[3] [10] = true;

  //bk   chi angles needed for building 2HZ
  AA[p].chi_required[0] [11] = true;
  AA[p].chi_required[1] [11] = true;
  AA[p].chi_required[2] [11] = true;
  AA[p].chi_required[3] [11] = true;

  //bk   chi angles needed for building 3HZ
  AA[p].chi_required[0] [12] = true;
  AA[p].chi_required[1] [12] = true;
  AA[p].chi_required[2] [12] = true;
  AA[p].chi_required[3] [12] = true;

  //bk   chi angles needed for building  HA

  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [14] = true;

  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [15] = true;

  //bk   chi angles needed for building 2HG
  AA[p].chi_required[0] [16] = true;
  AA[p].chi_required[1] [16] = true;


  //bk   chi angles needed for building 3HG
  AA[p].chi_required[0] [17] = true;
  AA[p].chi_required[1] [17] = true;


  //bk   chi angles needed for building 2HD
  AA[p].chi_required[0] [18] = true;
  AA[p].chi_required[1] [18] = true;
  AA[p].chi_required[2] [18] = true;


  //bk   chi angles needed for building 3HD
  AA[p].chi_required[0] [19] = true;
  AA[p].chi_required[1] [19] = true;
  AA[p].chi_required[2] [19] = true;


  //bk   chi angles needed for building 2HE
  AA[p].chi_required[0] [20] = true;
  AA[p].chi_required[1] [20] = true;
  AA[p].chi_required[2] [20] = true;
  AA[p].chi_required[3] [20] = true;

  //bk   chi angles needed for building 3HE
  AA[p].chi_required[0] [21] = true;
  AA[p].chi_required[1] [21] = true;
  AA[p].chi_required[2] [21] = true;
  AA[p].chi_required[3] [21] = true;


  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   CG
  //bk   four atoms that define chi angle  2
  AA[p].chi_atoms[1] [0] = 1; //   CA
  AA[p].chi_atoms[1] [1] = 4; //   CB
  AA[p].chi_atoms[1] [2] = 5; //   CG
  AA[p].chi_atoms[1] [3] = 6; //   CD
  //bk   four atoms that define chi angle  3
  AA[p].chi_atoms[2] [0] = 4; //   CB
  AA[p].chi_atoms[2] [1] = 5; //   CG
  AA[p].chi_atoms[2] [2] = 6; //   CD
  AA[p].chi_atoms[2] [3] = 7; //   CE
  //bk   four atoms that define chi angle  4
  AA[p].chi_atoms[3] [0] = 5; //   CG
  AA[p].chi_atoms[3] [1] = 6; //   CD
  AA[p].chi_atoms[3] [2] = 7; //   CE
  AA[p].chi_atoms[3] [3] = 8; //   NZ


  AA[p].atom[0].icoor[0] = 0.0000; //   N
  AA[p].atom[0].icoor[1] = 13.9730; //   N
  AA[p].atom[0].icoor[2] = 29.3820; //   N
  AA[p].atom[1].icoor[0] = 0.0000; //   CA
  AA[p].atom[1].icoor[1] = 15.0130; //   CA
  AA[p].atom[1].icoor[2] = 30.4040; //   CA
  AA[p].atom[2].icoor[0] = 0.0000; //   C
  AA[p].atom[2].icoor[1] = 14.4110; //   C
  AA[p].atom[2].icoor[2] = 31.8030; //   C
  AA[p].atom[3].icoor[0] = 0.0030; //   O
  AA[p].atom[3].icoor[1] = 13.1890; //   O
  AA[p].atom[3].icoor[2] = 31.9540; //   O
  AA[p].atom[4].icoor[0] = -1.2070; //   CB
  AA[p].atom[4].icoor[1] = 15.9360; //   CB
  AA[p].atom[4].icoor[2] = 30.2290; //   CB
  AA[p].atom[5].icoor[0] = -1.2940; //   CG
  AA[p].atom[5].icoor[1] = 17.0590; //   CG
  AA[p].atom[5].icoor[2] = 31.2540; //   CG
  AA[p].atom[6].icoor[0] = -2.5170; //   CD
  AA[p].atom[6].icoor[1] = 17.9320; //   CD
  AA[p].atom[6].icoor[2] = 31.0160; //   CD
  AA[p].atom[7].icoor[0] = -2.6050; //   CE
  AA[p].atom[7].icoor[1] = 19.0540; //   CE
  AA[p].atom[7].icoor[2] = 32.0400; //   CE
  AA[p].atom[8].icoor[0] = -3.7970; //   NZ
  AA[p].atom[8].icoor[1] = 19.9170; //   NZ
  AA[p].atom[8].icoor[2] = 31.8190; //   NZ
  AA[p].atom[9].icoor[0] = -0.0400; //   H
  AA[p].atom[9].icoor[1] = 13.0090; //   H
  AA[p].atom[9].icoor[2] = 29.6800; //   H
  AA[p].atom[10].icoor[0] = -3.8180; //  1HZ
  AA[p].atom[10].icoor[1] = 20.6470; //  1HZ
  AA[p].atom[10].icoor[2] = 32.5170; //  1HZ
  AA[p].atom[11].icoor[0] = -3.7470; //  2HZ
  AA[p].atom[11].icoor[1] = 20.3320; //  2HZ
  AA[p].atom[11].icoor[2] = 30.8990; //  2HZ
  AA[p].atom[12].icoor[0] = -4.6350; //  3HZ
  AA[p].atom[12].icoor[1] = 19.3580; //  3HZ
  AA[p].atom[12].icoor[2] = 31.8910; //  3HZ
  AA[p].atom[13].icoor[0] = 0.9090; //   HA
  AA[p].atom[13].icoor[1] = 15.6100; //   HA
  AA[p].atom[13].icoor[2] = 30.3230; //   HA
  AA[p].atom[14].icoor[0] = -1.1400; //  2HB
  AA[p].atom[14].icoor[1] = 16.3630; //  2HB
  AA[p].atom[14].icoor[2] = 29.2280; //  2HB
  AA[p].atom[15].icoor[0] = -2.0990; //  3HB
  AA[p].atom[15].icoor[1] = 15.3130; //  3HB
  AA[p].atom[15].icoor[2] = 30.2980; //  3HB
  AA[p].atom[16].icoor[0] = -1.3510; //  2HG
  AA[p].atom[16].icoor[1] = 16.6160; //  2HG
  AA[p].atom[16].icoor[2] = 32.2490; //  2HG
  AA[p].atom[17].icoor[0] = -0.3930; //  3HG
  AA[p].atom[17].icoor[1] = 17.6670; //  3HG
  AA[p].atom[17].icoor[2] = 31.1790; //  3HG
  AA[p].atom[18].icoor[0] = -2.4500; //  2HD
  AA[p].atom[18].icoor[1] = 18.3580; //  2HD
  AA[p].atom[18].icoor[2] = 30.0140; //  2HD
  AA[p].atom[19].icoor[0] = -3.4080; //  3HD
  AA[p].atom[19].icoor[1] = 17.3080; //  3HD
  AA[p].atom[19].icoor[2] = 31.0850; //  3HD
  AA[p].atom[20].icoor[0] = -2.6590; //  2HE
  AA[p].atom[20].icoor[1] = 18.6080; //  2HE
  AA[p].atom[20].icoor[2] = 33.0320; //  2HE
  AA[p].atom[21].icoor[0] = -1.7010; //  3HE
  AA[p].atom[21].icoor[1] = 19.6580; //  3HE
  AA[p].atom[21].icoor[2] = 31.9630; //  3HE


  // atom number for backbone HN
  AA[p].HNpos=9;
  // atom number for backbone HA
  AA[p].HApos=13;

  // number of polar hydrogens
  AA[p].nH_polar_complete=4;
  // atom numbers for polar H
  AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
  AA[p].Hpos_polar_complete[0]=9;
  AA[p].Hpos_polar_complete[1]=10;
  AA[p].Hpos_polar_complete[2]=11;
  AA[p].Hpos_polar_complete[3]=12;

  // number of aromatic hydrogens
  AA[p].nH_aromatic_complete=0;

  // number of apolar hydrogens
  AA[p].nH_apolar_complete=9;
  //atom number for apolar hydrogens
  AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
  AA[p].Hpos_apolar_complete[0]=13;
  AA[p].Hpos_apolar_complete[1]=14;
  AA[p].Hpos_apolar_complete[2]=15;
  AA[p].Hpos_apolar_complete[3]=16;
  AA[p].Hpos_apolar_complete[4]=17;
  AA[p].Hpos_apolar_complete[5]=18;
  AA[p].Hpos_apolar_complete[6]=19;
  AA[p].Hpos_apolar_complete[7]=20;
  AA[p].Hpos_apolar_complete[8]=21;


  // number of acceptors
  AA[p].nacceptors=1;
  //acceptor information
  AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
  (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;

  //Number of Hydrogens connection
  AA[p].nH_hydrogen_connexions=13;
  //atoms hydrogens are connected too
  AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
  AA[p].Hbase[0][0]=9; AA[p].Hbase[0][1]=0;
  AA[p].Hbase[1][0]=10; AA[p].Hbase[1][1]=8;
  AA[p].Hbase[2][0]=11; AA[p].Hbase[2][1]=8;
  AA[p].Hbase[3][0]=12; AA[p].Hbase[3][1]=8;
  AA[p].Hbase[4][0]=13; AA[p].Hbase[4][1]=1;
  AA[p].Hbase[5][0]=14; AA[p].Hbase[5][1]=4;
  AA[p].Hbase[6][0]=15; AA[p].Hbase[6][1]=4;
  AA[p].Hbase[7][0]=16; AA[p].Hbase[7][1]=5;
  AA[p].Hbase[8][0]=17; AA[p].Hbase[8][1]=5;
  AA[p].Hbase[9][0]=18; AA[p].Hbase[9][1]=6;
  AA[p].Hbase[10][0]=19; AA[p].Hbase[10][1]=6;
  AA[p].Hbase[11][0]=20; AA[p].Hbase[11][1]=7;
  AA[p].Hbase[12][0]=21; AA[p].Hbase[12][1]=7;


  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=9;
  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=13;
  AA[p].atom[4].numHydrogens_atm=2;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=14;
  AA[p].atom[4].hydrogens_atm[1]=15;
  AA[p].atom[5].numHydrogens_atm=2;
  AA[p].atom[5].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[5].numHydrogens_atm);
  AA[p].atom[5].hydrogens_atm[0]=16;
  AA[p].atom[5].hydrogens_atm[1]=17;
  AA[p].atom[6].numHydrogens_atm=2;
  AA[p].atom[6].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[6].numHydrogens_atm);
  AA[p].atom[6].hydrogens_atm[0]=18;
  AA[p].atom[6].hydrogens_atm[1]=19;
  AA[p].atom[7].numHydrogens_atm=2;
  AA[p].atom[7].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[7].numHydrogens_atm);
  AA[p].atom[7].hydrogens_atm[0]=20;
  AA[p].atom[7].hydrogens_atm[1]=21;
  AA[p].atom[8].numHydrogens_atm=3;
  AA[p].atom[8].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[8].numHydrogens_atm);
  AA[p].atom[8].hydrogens_atm[0]=10;
  AA[p].atom[8].hydrogens_atm[1]=11;
  AA[p].atom[8].hydrogens_atm[2]=12;

}

//FIN DE AMINOACIDO  LYSINA

/**
* Initializes Leucine
*
* @param opt: Rosseta or ICM
*/
void init_LEU(Convention opt)
{
  //AMINOACIDO LEUCINA
  int k, p =LEU;


  strcpy( AA[p].aa_name3, "LEU" );
  AA[p].aa_name1 = 'L';
  AA[p].mass = 113.161;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = false;
  AA[p].aa_is_nonpolar = true;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = false;

  AA[p].natoms = 19;
  AA[p].nheavyatoms = 8;
  AA[p].nchi = 2;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );
  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;

  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

  AA[p].chi_types = ( int * ) new int[AA[p].nchi];
  AA[p].chi_types[0] = 2;
  AA[p].chi_types[1] = 2;
  AA[p].natoms_EEF1 = 8;


  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " CG " );
  strcpy( AA[p].atom[6].atom_name, " CD1" );
  strcpy( AA[p].atom[7].atom_name, " CD2" );
  strcpy( AA[p].atom[8].atom_name, " H  " );
  strcpy( AA[p].atom[9].atom_name, " HA " );
  strcpy( AA[p].atom[10].atom_name, "2HB " );
  strcpy( AA[p].atom[11].atom_name, "3HB " );
  strcpy( AA[p].atom[12].atom_name, " HG " );
  strcpy( AA[p].atom[13].atom_name, "1HD1" );
  strcpy( AA[p].atom[14].atom_name, "2HD1" );
  strcpy( AA[p].atom[15].atom_name, "3HD1" );
  strcpy( AA[p].atom[16].atom_name, "1HD2" );
  strcpy( AA[p].atom[17].atom_name, "2HD2" );
  strcpy( AA[p].atom[18].atom_name, "3HD2" );

  strcpy( AA[p].atom[19].atom_name," OXT" );
  strcpy( AA[p].atom[20].atom_name,"2H  " );
  strcpy( AA[p].atom[21].atom_name,"3H  " );

  switch(opt)
  {
    case Rosseta:
    case Sybil:
    AA[p].atom[0].charge = -0.47; //    N
      AA[p].atom[1].charge = 0.07; //    CA
      AA[p].atom[2].charge = 0.51; //    C
      AA[p].atom[3].charge = -0.51; //    O
      AA[p].atom[4].charge = -0.18; //    CB
      AA[p].atom[5].charge = -0.09; //    CG
      AA[p].atom[6].charge = -0.27; //    CD1
      AA[p].atom[7].charge = -0.27; //    CD2
      AA[p].atom[8].charge = 0.31; //    H
      AA[p].atom[9].charge = 0.09; //    HA
      AA[p].atom[10].charge = 0.09; //    2HB
      AA[p].atom[11].charge = 0.09; //    3HB
      AA[p].atom[12].charge = 0.09; //    HG
      AA[p].atom[13].charge = 0.09; //    1HD1
      AA[p].atom[14].charge = 0.09; //    2HD1
      AA[p].atom[15].charge = 0.09; //    3HD1
      AA[p].atom[16].charge = 0.09; //    1HD2
      AA[p].atom[17].charge = 0.09; //    2HD2
      AA[p].atom[18].charge = 0.09; //    3HD2
      AA[p].atom[19].charge = -0.67; //   0XT
      AA[p].atom[20].charge = 0.33; //  1H
      AA[p].atom[21].charge = 0.33; //  2H
      break;
    case ICM:
      AA[p].atom[0].charge = -0.356; //    N
      AA[p].atom[1].charge = 0.064; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = -0.03; //    CB
      AA[p].atom[5].charge = -0.011; //    CG
      AA[p].atom[6].charge = -0.074; //    CD1
      AA[p].atom[7].charge = -0.074; //    CD2
      AA[p].atom[8].charge = 0.176; //    H
      AA[p].atom[9].charge = 0.02; //    HA
      AA[p].atom[10].charge = 0.022; //    2HB
      AA[p].atom[11].charge = 0.022; //    3HB
      AA[p].atom[12].charge = 0.025; //    HG
      AA[p].atom[13].charge = 0.025; //    1HD1
      AA[p].atom[14].charge = 0.025; //    2HD1
      AA[p].atom[15].charge = 0.025; //    3HD1
      AA[p].atom[16].charge = 0.025; //    1HD2
      AA[p].atom[17].charge = 0.025; //    2HD2
      AA[p].atom[18].charge = 0.025; //    3HD2

      AA[p].atom[19].charge = -0.67; //   0XT
      AA[p].atom[20].charge = 0.33; //  1H
      AA[p].atom[21].charge = 0.33; //  2H
      break;
     case EEF1:
      AA[p].atom[0].charge = -0.35; //    N
      AA[p].atom[1].charge =  0.10; //    CA
      AA[p].atom[2].charge =  0.55; //    C
      AA[p].atom[3].charge = -0.55; //    O
      AA[p].atom[4].charge =  0.00; //    CB
      AA[p].atom[5].charge =  0.00; //    CG
      AA[p].atom[6].charge =  0.00; //    CD1
      AA[p].atom[7].charge =  0.00; //    CD2
      AA[p].atom[8].charge =  0.25; //    H
      AA[p].atom[9].charge =  0.00; //    HA
      AA[p].atom[10].charge = 0.00; //    2HB
      AA[p].atom[11].charge = 0.00; //    3HB
      AA[p].atom[12].charge = 0.00; //    HG
      AA[p].atom[13].charge = 0.00; //    1HD1
      AA[p].atom[14].charge = 0.00; //    2HD1
      AA[p].atom[15].charge = 0.00; //    3HD1
      AA[p].atom[16].charge = 0.00; //    1HD2
      AA[p].atom[17].charge = 0.00; //    2HD2
      AA[p].atom[18].charge = 0.00; //    3HD2
      break;
  }
  /*
   1.0 LEU-N:2H    0.33000 HC
   1.0 LEU-N:3HD2  0.09000 HA
   1.0 LEU-N:HA    0.10000 HB
   1.0 LEU-N:3H    0.33000 HC
   1.0 LEU-N:C     0.51000 C
   1.0 LEU-N:CB   -0.18000 CT2
   1.0 LEU-N:O    -0.51000 O
   1.0 LEU-N:HB1   0.09000 HA
   1.0 LEU-N:1HB   0.09000 HA
   1.0 LEU-N:CG   -0.09000 CT1
   1.0 LEU-N:HG    0.09000 HA
   1.0 LEU-N:CD1  -0.27000 CT3
   1.0 LEU-N:1HD1  0.09000 HA
   1.0 LEU-N:2HD1  0.09000 HA
   1.0 LEU-N:3HD1  0.09000 HA
   1.0 LEU-N:CD2  -0.27000 CT3
   1.0 LEU-N:N    -0.30000 NH3
   1.0 LEU-N:1HD2  0.09000 HA
   1.0 LEU-N:1H    0.33000 HC
   1.0 LEU-N:2HD2  0.09000 HA
   1.0 LEU-N:CA    0.21000 CT1

   1.0 LEU-C:3HD2  0.09000 HA
   1.0 LEU-C:HA    0.09000 HB
   1.0 LEU-C:C     0.34000 CC
   1.0 LEU-C:CB   -0.18000 CT2
   1.0 LEU-C:HB1   0.09000 HA
   1.0 LEU-C:1HB   0.09000 HA
   1.0 LEU-C:CG   -0.09000 CT1
   1.0 LEU-C:HG    0.09000 HA
   1.0 LEU-C:CD1  -0.27000 CT3
   1.0 LEU-C:1HD1  0.09000 HA
   1.0 LEU-C:O    -0.67000 OC
   1.0 LEU-C:2HD1  0.09000 HA
   1.0 LEU-C:OXT  -0.67000 OC
   1.0 LEU-C:3HD1  0.09000 HA
   1.0 LEU-C:CD2  -0.27000 CT3
   1.0 LEU-C:N    -0.47000 NH1
   1.0 LEU-C:1HD2  0.09000 HA
   1.0 LEU-C:HN    0.31000 H
   1.0 LEU-C:2HD2  0.09000 HA
   1.0 LEU-C:CA    0.07000 CT1
*/

  switch (opt) {
  	case EEF1:
	  AA[p].atom[0].fullatom_type =  12; // NH1    N
	  AA[p].atom[1].fullatom_type =   5; // CH1E   CA
	  AA[p].atom[2].fullatom_type =   4; // C      C
	  AA[p].atom[3].fullatom_type =  16; // O      O
	  AA[p].atom[4].fullatom_type =   6; // CH2E   CB
	  AA[p].atom[5].fullatom_type =   5; // CH1E   CG
	  AA[p].atom[6].fullatom_type =   7; // CH3E   CD1
	  AA[p].atom[7].fullatom_type =   7; // CH3E   CD2
	  AA[p].atom[8].fullatom_type =   0; // H      H
	  AA[p].atom[9].fullatom_type =  -1; // Hapo   HA
	  AA[p].atom[10].fullatom_type = -1; // Hapo   2HB
	  AA[p].atom[11].fullatom_type = -1; // Hapo   3HB
	  AA[p].atom[12].fullatom_type = -1; // Hapo   HG
	  AA[p].atom[13].fullatom_type = -1; // Hapo   1HD1
	  AA[p].atom[14].fullatom_type = -1; // Hapo   2HD1
	  AA[p].atom[15].fullatom_type = -1; // Hapo   3HD1
	  AA[p].atom[16].fullatom_type = -1; // Hapo   1HD2
	  AA[p].atom[17].fullatom_type = -1; // Hapo   2HD2
	  AA[p].atom[18].fullatom_type = -1; // Hapo   3HD2

	  AA[p].atom[19].fullatom_type = 17; // OC     OXT
	  AA[p].atom[20].fullatom_type =  1; // HC     2H
	  AA[p].atom[21].fullatom_type =  1; // HC     3H
	  break;
  	case Sybil:
  		  AA[p].atom[0].fullatom_type = 40; // Nbb    N
  		  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
  		  AA[p].atom[2].fullatom_type = 10; // CObb   C
  		  AA[p].atom[3].fullatom_type = 56; // OCbb   O
  		  AA[p].atom[4].fullatom_type = 8; // CH2    CB
  		  AA[p].atom[5].fullatom_type = 9; // CH1    CG
  		  AA[p].atom[6].fullatom_type = 6; // CH3    CD1
  		  AA[p].atom[7].fullatom_type = 6; // CH3    CD2
  		  AA[p].atom[8].fullatom_type = 66; // HNbb   H
  		  AA[p].atom[9].fullatom_type = 66; // Hapo   HA
  		  AA[p].atom[10].fullatom_type = 66; // Hapo  2HB
  		  AA[p].atom[11].fullatom_type = 66; // Hapo  3HB
  		  AA[p].atom[12].fullatom_type = 66; // Hapo   HG
  		  AA[p].atom[13].fullatom_type = 66; // Hapo  1HD1
  		  AA[p].atom[14].fullatom_type = 66; // Hapo  2HD1
  		  AA[p].atom[15].fullatom_type = 66; // Hapo  3HD1
  		  AA[p].atom[16].fullatom_type = 66; // Hapo  1HD2
  		  AA[p].atom[17].fullatom_type = 66; // Hapo  2HD2
  		  AA[p].atom[18].fullatom_type = 66; // Hapo  3HD2

  		  AA[p].atom[19].fullatom_type = 57; // OXT
  		  AA[p].atom[20].fullatom_type = 66; // Hpol
  		  AA[p].atom[21].fullatom_type = 66; // Hpol
  		  break;
  default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 3; // CH1    CG
	  AA[p].atom[6].fullatom_type = 5; // CH3    CD1
	  AA[p].atom[7].fullatom_type = 5; // CH3    CD2
	  AA[p].atom[8].fullatom_type = 25; // HNbb   H
	  AA[p].atom[9].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[10].fullatom_type = 23; // Hapo  2HB
	  AA[p].atom[11].fullatom_type = 23; // Hapo  3HB
	  AA[p].atom[12].fullatom_type = 23; // Hapo   HG
	  AA[p].atom[13].fullatom_type = 23; // Hapo  1HD1
	  AA[p].atom[14].fullatom_type = 23; // Hapo  2HD1
	  AA[p].atom[15].fullatom_type = 23; // Hapo  3HD1
	  AA[p].atom[16].fullatom_type = 23; // Hapo  1HD2
	  AA[p].atom[17].fullatom_type = 23; // Hapo  2HD2
	  AA[p].atom[18].fullatom_type = 23; // Hapo  3HD2

	  AA[p].atom[19].fullatom_type = 20; // OXT
	  AA[p].atom[20].fullatom_type = 22; // Hpol
	  AA[p].atom[21].fullatom_type = 22; // Hpol
	  break;
  }

  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 8; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 9; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--CG
  AA[p].atom[4].bonded_neighbor[2] = 10; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 11; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 4; // CG
  AA[p].atom[5].bonded_neighbor[0] = 4; // CG--CB
  AA[p].atom[5].bonded_neighbor[1] = 6; // CG--CD1
  AA[p].atom[5].bonded_neighbor[2] = 7; // CG--CD2
  AA[p].atom[5].bonded_neighbor[3] = 12; // CG--HG
  AA[p].atom[6].nbonded_neighbors = 4; // CD1
  AA[p].atom[6].bonded_neighbor[0] = 5; // CD1--CG
  AA[p].atom[6].bonded_neighbor[1] = 13; // CD1--1HD1
  AA[p].atom[6].bonded_neighbor[2] = 14; // CD1--2HD1
  AA[p].atom[6].bonded_neighbor[3] = 15; // CD1--3HD1
  AA[p].atom[7].nbonded_neighbors = 4; // CD2
  AA[p].atom[7].bonded_neighbor[0] = 5; // CD2--CG
  AA[p].atom[7].bonded_neighbor[1] = 16; // CD2--1HD2
  AA[p].atom[7].bonded_neighbor[2] = 17; // CD2--2HD2
  AA[p].atom[7].bonded_neighbor[3] = 18; // CD2--3HD2
  AA[p].atom[8].nbonded_neighbors = 1; // H
  AA[p].atom[8].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[9].nbonded_neighbors = 1; // HA
  AA[p].atom[9].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[10].nbonded_neighbors = 1; //2HB
  AA[p].atom[10].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[11].nbonded_neighbors = 1; //3HB
  AA[p].atom[11].bonded_neighbor[0] = 4; //3HB--CB
  AA[p].atom[12].nbonded_neighbors = 1; // HG
  AA[p].atom[12].bonded_neighbor[0] = 5; // HG--CG
  AA[p].atom[13].nbonded_neighbors = 1; //1HD1
  AA[p].atom[13].bonded_neighbor[0] = 6; //1HD1--CD1
  AA[p].atom[14].nbonded_neighbors = 1; //2HD1
  AA[p].atom[14].bonded_neighbor[0] = 6; //2HD1--CD1
  AA[p].atom[15].nbonded_neighbors = 1; //3HD1
  AA[p].atom[15].bonded_neighbor[0] = 6; //3HD1--CD1
  AA[p].atom[16].nbonded_neighbors = 1; //1HD2
  AA[p].atom[16].bonded_neighbor[0] = 7; //1HD2--CD2
  AA[p].atom[17].nbonded_neighbors = 1; //2HD2
  AA[p].atom[17].bonded_neighbor[0] = 7; //2HD2--CD2
  AA[p].atom[18].nbonded_neighbors = 1; //3HD2
  AA[p].atom[18].bonded_neighbor[0] = 7; //3HD2--CD2

  AA[p].atom[19].nbonded_neighbors = 1; // OXT
  AA[p].atom[19].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[20].nbonded_neighbors = 1; // 1H
  AA[p].atom[20].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[21].nbonded_neighbors = 1; // 2H
  AA[p].atom[21].bonded_neighbor[0] = 0; // H-N


  for ( k = 0; k < 9; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 9; k < 19; k++ )
    AA[p].atom[k].hastemplate = true;




  //bk   template for building  HA
  AA[p].atom[9].ta[0] = 1; //   CA
  AA[p].atom[9].ta[1] = 0; //   N
  AA[p].atom[9].ta[2] = 2; //   C

  //bk   template for building 2HB
  AA[p].atom[10].ta[0] = 4; //   CB
  AA[p].atom[10].ta[1] = 1; //   CA
  AA[p].atom[10].ta[2] = 5; //   CG

  //bk   template for building 3HB
  AA[p].atom[11].ta[0] = 4; //   CB
  AA[p].atom[11].ta[1] = 1; //   CA
  AA[p].atom[11].ta[2] = 5; //   CG

  //bk   template for building HG
  AA[p].atom[12].ta[0] = 5; //   CG
  AA[p].atom[12].ta[1] = 4; //   CB
  AA[p].atom[12].ta[2] = 6; //   CD1

  //bk   template for building 1HD1
  AA[p].atom[13].ta[0] = 6; //   CD1
  AA[p].atom[13].ta[1] = 5; //   CG
  AA[p].atom[13].ta[2] = 4; //   CB

  //bk   template for building 2HD1
  AA[p].atom[14].ta[0] = 6; //   CD1
  AA[p].atom[14].ta[1] = 5; //   CG
  AA[p].atom[14].ta[2] = 4; //   CB

  //bk   template for building 3HD1
  AA[p].atom[15].ta[0] = 6; //   CD1
  AA[p].atom[15].ta[1] = 5; //   CG
  AA[p].atom[15].ta[2] = 4; //   CB

  //bk   template for building 1HD2
  AA[p].atom[16].ta[0] = 7; //   CD2
  AA[p].atom[16].ta[1] = 5; //   CG
  AA[p].atom[16].ta[2] = 4; //   CB

  //bk   template for building 2HD2
  AA[p].atom[17].ta[0] = 7; //   CD2
  AA[p].atom[17].ta[1] = 5; //   CG
  AA[p].atom[17].ta[2] = 4; //   CB

  //bk   template for building 3HD2
  AA[p].atom[18].ta[0] = 7; //   CD2
  AA[p].atom[18].ta[1] = 5; //   CG
  AA[p].atom[18].ta[2] = 4; //   CB



  //bk   chi angles required to build atoms LEU



  //bk   chi angles needed for building  CG
  AA[p].chi_required[0] [5] = true;

  //bk   chi angles needed for building  CD1
  AA[p].chi_required[0] [6] = true;
  AA[p].chi_required[1] [6] = true;

  //bk   chi angles needed for building  CD2
  AA[p].chi_required[0] [7] = true;
  AA[p].chi_required[1] [7] = true;

  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [10] = true;

  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [11] = true;

  //bk   chi angles needed for building  HG
  AA[p].chi_required[0] [12] = true;
  AA[p].chi_required[1] [12] = true;

  //bk   chi angles needed for building 1HD1
  AA[p].chi_required[0] [13] = true;
  AA[p].chi_required[1] [13] = true;

  //bk   chi angles needed for building 2HD1
  AA[p].chi_required[0] [14] = true;
  AA[p].chi_required[1] [14] = true;

  //bk   chi angles needed for building 3HD1
  AA[p].chi_required[0] [15] = true;
  AA[p].chi_required[1] [15] = true;

  //bk   chi angles needed for building 1HD2
  AA[p].chi_required[0] [16] = true;
  AA[p].chi_required[1] [16] = true;

  //bk   chi angles needed for building 2HD2
  AA[p].chi_required[0] [17] = true;
  AA[p].chi_required[1] [17] = true;

  //bk   chi angles needed for building 3HD2
  AA[p].chi_required[0] [18] = true;
  AA[p].chi_required[1] [18] = true;




  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   CG

  //bk   four atoms that define chi angle  2
  AA[p].chi_atoms[1] [0] = 1; //   CA
  AA[p].chi_atoms[1] [1] = 4; //   CB
  AA[p].chi_atoms[1] [2] = 5; //   CG
  AA[p].chi_atoms[1] [3] = 6; //   CD1



  AA[p].atom[0].icoor[0] = 1.460; //   N
  AA[p].atom[0].icoor[1] = 0.000; //   N
  AA[p].atom[0].icoor[2] = 0.000; //   N
  AA[p].atom[1].icoor[0] = 0.000; //   CA
  AA[p].atom[1].icoor[1] = 0.000; //   CA
  AA[p].atom[1].icoor[2] = 0.000; //   CA
  AA[p].atom[2].icoor[0] = -0.559; //   C
  AA[p].atom[2].icoor[1] = 1.422; //   C
  AA[p].atom[2].icoor[2] = -0.000; //   C
  AA[p].atom[3].icoor[0] = -1.545; //   O
  AA[p].atom[3].icoor[1] = 1.703; //   O
  AA[p].atom[3].icoor[2] = -0.678; //   O
  AA[p].atom[4].icoor[0] = -0.530; //   CB
  AA[p].atom[4].icoor[1] = -0.778; //   CB
  AA[p].atom[4].icoor[2] = 1.211; //   CB
  AA[p].atom[5].icoor[0] = -0.286; //   CG
  AA[p].atom[5].icoor[1] = -2.292; //   CG
  AA[p].atom[5].icoor[2] = 1.172; //   CG
  AA[p].atom[6].icoor[0] = -0.695; //   CD1
  AA[p].atom[6].icoor[1] = -2.921; //   CD1
  AA[p].atom[6].icoor[2] = 2.497; //   CD1
  AA[p].atom[7].icoor[0] = -1.069; //   CD2
  AA[p].atom[7].icoor[1] = -2.904; //   CD2
  AA[p].atom[7].icoor[2] = 0.020; //   CD2
  AA[p].atom[8].icoor[0] = 2.424; //   H
  AA[p].atom[8].icoor[1] = -0.301; //   H
  AA[p].atom[8].icoor[2] = 0.000; //   H
  AA[p].atom[9].icoor[0] = -0.365; //   HA
  AA[p].atom[9].icoor[1] = -0.473; //   HA
  AA[p].atom[9].icoor[2] = -0.911; //   HA
  AA[p].atom[10].icoor[0] = 0.075; //  1HB
  AA[p].atom[10].icoor[1] = -0.322; //  1HB
  AA[p].atom[10].icoor[2] = 1.993; //  1HB
  AA[p].atom[11].icoor[0] = -1.583; //  2HB
  AA[p].atom[11].icoor[1] = -0.567; //  2HB
  AA[p].atom[11].icoor[2] = 1.397; //  2HB
  AA[p].atom[12].icoor[0] = 0.775; //   HG
  AA[p].atom[12].icoor[1] = -2.444; //   HG
  AA[p].atom[12].icoor[2] = 0.972; //   HG
  AA[p].atom[13].icoor[0] = -0.518; //  1HD1
  AA[p].atom[13].icoor[1] = -3.996; //  1HD1
  AA[p].atom[13].icoor[2] = 2.460; //  1HD1
  AA[p].atom[14].icoor[0] = -0.106; //  2HD1
  AA[p].atom[14].icoor[1] = -2.485; //  2HD1
  AA[p].atom[14].icoor[2] = 3.304; //  2HD1
  AA[p].atom[15].icoor[0] = -1.753; //  3HD1
  AA[p].atom[15].icoor[1] = -2.734; //  3HD1
  AA[p].atom[15].icoor[2] = 2.677; //  3HD1
  AA[p].atom[16].icoor[0] = -0.894; //  1HD2
  AA[p].atom[16].icoor[1] = -3.980; //  1HD2
  AA[p].atom[16].icoor[2] = -0.007; //  1HD2
  AA[p].atom[17].icoor[0] = -2.133; //  2HD2
  AA[p].atom[17].icoor[1] = -2.712; //  2HD2
  AA[p].atom[17].icoor[2] = 0.160; //  2HD2
  AA[p].atom[18].icoor[0] = -0.741; //  3HD2
  AA[p].atom[18].icoor[1] = -2.460; //  3HD2
  AA[p].atom[18].icoor[2] = -0.920; //  3HD2

  // atom number for backbone HN
  AA[p].HNpos=8;
  // atom number for backbone HA
  AA[p].HApos=9;

  // number of polar hydrogens
  AA[p].nH_polar_complete=1;
  // atom numbers for polar H
  AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
  AA[p].Hpos_polar_complete[0]=8;

  // number of aromatic hydrogens
  AA[p].nH_aromatic_complete=0;

  // number of apolar hydrogens
  AA[p].nH_apolar_complete=10;
  //atom number for apolar hydrogens
  AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
  AA[p].Hpos_apolar_complete[0]=9;
  AA[p].Hpos_apolar_complete[1]=10;
  AA[p].Hpos_apolar_complete[2]=11;
  AA[p].Hpos_apolar_complete[3]=12;
  AA[p].Hpos_apolar_complete[4]=13;
  AA[p].Hpos_apolar_complete[5]=14;
  AA[p].Hpos_apolar_complete[6]=15;
  AA[p].Hpos_apolar_complete[7]=16;
  AA[p].Hpos_apolar_complete[8]=17;
  AA[p].Hpos_apolar_complete[9]=18;

  // number of acceptors
  AA[p].nacceptors=1;
  //acceptor information
  AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
  (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;

  //Number of Hydrogens connection
  AA[p].nH_hydrogen_connexions=11;
  //atoms hydrogens are connected too
  AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
  AA[p].Hbase[0][0]=8; AA[p].Hbase[0][1]=0;
  AA[p].Hbase[1][0]=9; AA[p].Hbase[1][1]=1;
  AA[p].Hbase[2][0]=10; AA[p].Hbase[2][1]=4;
  AA[p].Hbase[3][0]=11; AA[p].Hbase[3][1]=4;
  AA[p].Hbase[4][0]=12; AA[p].Hbase[4][1]=5;
  AA[p].Hbase[5][0]=13; AA[p].Hbase[5][1]=6;
  AA[p].Hbase[6][0]=14; AA[p].Hbase[6][1]=6;
  AA[p].Hbase[7][0]=15; AA[p].Hbase[7][1]=6;
  AA[p].Hbase[8][0]=16; AA[p].Hbase[8][1]=7;
  AA[p].Hbase[9][0]=17; AA[p].Hbase[9][1]=7;
  AA[p].Hbase[10][0]=18; AA[p].Hbase[10][1]=7;


  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=8;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=9;

  AA[p].atom[4].numHydrogens_atm=2;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=10;
  AA[p].atom[4].hydrogens_atm[1]=11;

  AA[p].atom[5].numHydrogens_atm=1;
  AA[p].atom[5].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[5].numHydrogens_atm);
  AA[p].atom[5].hydrogens_atm[0]=12;

  AA[p].atom[6].numHydrogens_atm=3;
  AA[p].atom[6].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[6].numHydrogens_atm);
  AA[p].atom[6].hydrogens_atm[0]=13;
  AA[p].atom[6].hydrogens_atm[1]=14;
  AA[p].atom[6].hydrogens_atm[2]=15;

  AA[p].atom[7].numHydrogens_atm=3;
  AA[p].atom[7].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[7].numHydrogens_atm);
  AA[p].atom[7].hydrogens_atm[0]=16;
  AA[p].atom[7].hydrogens_atm[1]=17;
  AA[p].atom[7].hydrogens_atm[2]=18;


  //FIN DE AMINOACIDO LEUCINA
}

/**
* Initializes Methionine
*
* @param opt: Rosseta or ICM
*/
void init_MET(Convention opt)
{
  //AMINOACIDO METHIONINE
  int k, p =MET;


  strcpy( AA[p].aa_name3, "MET" );
  AA[p].aa_name1 = 'M';

  AA[p].mass = 131.2;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = false;
  AA[p].aa_is_nonpolar = true;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = false;

  AA[p].natoms = 17;
  AA[p].nheavyatoms = 8;
  AA[p].nchi = 3;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );
  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;

  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

  AA[p].chi_types = ( int * ) new int[AA[p].nchi];
  AA[p].chi_types[0] = 2;
  AA[p].chi_types[1] = 2;
  AA[p].chi_types[2] = 5;
  AA[p].natoms_EEF1 = 8;


  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " CG " );
  strcpy( AA[p].atom[6].atom_name, " SD " );
  strcpy( AA[p].atom[7].atom_name, " CE " );
  strcpy( AA[p].atom[8].atom_name, " H  " );
  strcpy( AA[p].atom[9].atom_name, " HA " );
  strcpy( AA[p].atom[10].atom_name, "2HB " );
  strcpy( AA[p].atom[11].atom_name, "3HB " );
  strcpy( AA[p].atom[12].atom_name, "2HG " );
  strcpy( AA[p].atom[13].atom_name, "3HG " );
  strcpy( AA[p].atom[14].atom_name, "1HE " );
  strcpy( AA[p].atom[15].atom_name, "2HE " );
  strcpy( AA[p].atom[16].atom_name, "3HE " );

  strcpy( AA[p].atom[17].atom_name," OXT" );
  strcpy( AA[p].atom[18].atom_name,"2H  " );
  strcpy( AA[p].atom[19].atom_name,"3H  " );


  switch(opt)
  {
    case Rosseta:
    case Sybil:
     AA[p].atom[0].charge = -0.47; //    N
      AA[p].atom[1].charge = 0.07; //    CA
      AA[p].atom[2].charge = 0.51; //    C
      AA[p].atom[3].charge = -0.51; //    O
      AA[p].atom[4].charge = -0.18; //    CB
      AA[p].atom[5].charge = -0.14; //    CG
      AA[p].atom[6].charge = -0.09; //    SD
      AA[p].atom[7].charge = -0.22; //    CE
      AA[p].atom[8].charge = 0.31; //    H
      AA[p].atom[9].charge = 0.09; //    HA
      AA[p].atom[10].charge = 0.09; //   2HB
      AA[p].atom[11].charge = 0.09; //   3HB
      AA[p].atom[12].charge = 0.09; //   2HG
      AA[p].atom[13].charge = 0.09; //   3HG
      AA[p].atom[14].charge = 0.09; //   1HE
      AA[p].atom[15].charge = 0.09; //   2HE
      AA[p].atom[16].charge = 0.09; //   3HE

      AA[p].atom[17].charge = -0.67; //   0XT
      AA[p].atom[18].charge = 0.33; //  1H
      AA[p].atom[19].charge = 0.33; //  2H
    break;

    case ICM:
      AA[p].atom[0].charge = -0.356; //    N
      AA[p].atom[1].charge = 0.064; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = -0.01; //    CB
      AA[p].atom[5].charge = -0.12; //    CG
      AA[p].atom[6].charge = -0.035; //    SD
      AA[p].atom[7].charge = -0.19; //    CE
      AA[p].atom[8].charge = 0.176; //    H
      AA[p].atom[9].charge = 0.02; //    HA
      AA[p].atom[10].charge = 0.03; //   2HB
      AA[p].atom[11].charge = 0.03; //   3HB
      AA[p].atom[12].charge = 0.045; //   2HG
      AA[p].atom[13].charge = 0.045; //   3HG
      AA[p].atom[14].charge = 0.055; //   1HE
      AA[p].atom[15].charge = 0.055; //   2HE
      AA[p].atom[16].charge = 0.055; //   3HE

      AA[p].atom[17].charge = -0.67; //   0XT
      AA[p].atom[18].charge = 0.33; //  1H
      AA[p].atom[19].charge = 0.33; //  2H
      break;
   case EEF1:
      AA[p].atom[0].charge = -0.35; //    N
      AA[p].atom[1].charge =  0.10; //    CA
      AA[p].atom[2].charge =  0.55; //    C
      AA[p].atom[3].charge = -0.55; //    O
      AA[p].atom[4].charge =  0.00; //    CB
      AA[p].atom[5].charge =  0.06; //    CG
      AA[p].atom[6].charge = -0.12; //    SD
      AA[p].atom[7].charge =  0.06; //    CE
      AA[p].atom[8].charge =  0.25; //    H
      AA[p].atom[9].charge =  0.00; //    HA
      AA[p].atom[10].charge = 0.00; //   2HB
      AA[p].atom[11].charge = 0.00; //   3HB
      AA[p].atom[12].charge = 0.00; //   2HG
      AA[p].atom[13].charge = 0.00; //   3HG
      AA[p].atom[14].charge = 0.00; //   1HE
      AA[p].atom[15].charge = 0.00; //   2HE
      AA[p].atom[16].charge = 0.00; //   3HE
      break;
  }
  /*
   1.0 MET-N:2H    0.33000 HC
   1.0 MET-N:O    -0.51000 O
   1.0 MET-N:HA    0.10000 HB
   1.0 MET-N:3H    0.33000 HC
   1.0 MET-N:CB   -0.18000 CT2
   1.0 MET-N:HB1   0.09000 HA
   1.0 MET-N:1HB   0.09000 HA
   1.0 MET-N:CG   -0.14000 CT2
   1.0 MET-N:HG1   0.09000 HA
   1.0 MET-N:1HG   0.09000 HA
   1.0 MET-N:SD   -0.09000 S
   1.0 MET-N:CE   -0.22000 CT3
   1.0 MET-N:1HE   0.09000 HA
   1.0 MET-N:2HE   0.09000 HA
   1.0 MET-N:N    -0.30000 NH3
   1.0 MET-N:3HE   0.09000 HA
   1.0 MET-N:1H    0.33000 HC
   1.0 MET-N:C     0.51000 C
   1.0 MET-N:CA    0.21000 CT1

   1.0 MET-C:HA    0.09000 HB
   1.0 MET-C:CB   -0.18000 CT2
   1.0 MET-C:HB1   0.09000 HA
   1.0 MET-C:1HB   0.09000 HA
   1.0 MET-C:CG   -0.14000 CT2
   1.0 MET-C:HG1   0.09000 HA
   1.0 MET-C:1HG   0.09000 HA
   1.0 MET-C:SD   -0.09000 S
   1.0 MET-C:O    -0.67000 OC
   1.0 MET-C:CE   -0.22000 CT3
   1.0 MET-C:OXT  -0.67000 OC
   1.0 MET-C:1HE   0.09000 HA
   1.0 MET-C:2HE   0.09000 HA
   1.0 MET-C:N    -0.47000 NH1
   1.0 MET-C:3HE   0.09000 HA
   1.0 MET-C:HN    0.31000 H
   1.0 MET-C:C     0.34000 CC
   1.0 MET-C:CA    0.07000 CT1
*/

  switch (opt) {
  	case EEF1:
	  AA[p].atom[0].fullatom_type =  12; // NH1    N
	  AA[p].atom[1].fullatom_type =   5; // CH1E   CA
	  AA[p].atom[2].fullatom_type =   4; // C      C
	  AA[p].atom[3].fullatom_type =  16; // O      O
	  AA[p].atom[4].fullatom_type =   6; // CH2E   CB
	  AA[p].atom[5].fullatom_type =   6; // CH2E   CG
	  AA[p].atom[6].fullatom_type =  20; // S      SD
	  AA[p].atom[7].fullatom_type =   7; // CH3E   CE
	  AA[p].atom[8].fullatom_type =   0; // H      H
	  AA[p].atom[9].fullatom_type =  -1; // Hapo   HA
	  AA[p].atom[10].fullatom_type = -1; // Hapo   2HB
	  AA[p].atom[11].fullatom_type = -1; // Hapo   3HB
	  AA[p].atom[12].fullatom_type = -1; // Hapo   2HG
	  AA[p].atom[13].fullatom_type = -1; // Hapo   3HG
	  AA[p].atom[14].fullatom_type = -1; // Hapo   1HE
	  AA[p].atom[15].fullatom_type = -1; // Hapo   2HE
	  AA[p].atom[16].fullatom_type = -1; // Hapo   3HE

	  AA[p].atom[17].fullatom_type = 17; // OC     OXT
	  AA[p].atom[18].fullatom_type =  1; // HC     Hpol
	  AA[p].atom[19].fullatom_type =  1; // HC     Hpol
	  break;
  	case Sybil:
 	  AA[p].atom[0].fullatom_type = 40; // Nbb    N
 	  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
 	  AA[p].atom[2].fullatom_type = 10; // CObb   C
 	  AA[p].atom[3].fullatom_type = 56; // OCbb   O
 	  AA[p].atom[4].fullatom_type = 8; // CH2    CB
 	  AA[p].atom[5].fullatom_type = 8; // CH2    CG
 	  AA[p].atom[6].fullatom_type = 63; // S      SD
 	  AA[p].atom[7].fullatom_type = 6; // CH3    CE
 	  AA[p].atom[8].fullatom_type = 66; // HNbb   H
 	  AA[p].atom[9].fullatom_type = 66; // Hapo   HA
 	  AA[p].atom[10].fullatom_type = 66; // Hapo  2HB
 	  AA[p].atom[11].fullatom_type = 66; // Hapo  3HB
 	  AA[p].atom[12].fullatom_type = 66; // Hapo  2HG
 	  AA[p].atom[13].fullatom_type = 66; // Hapo  3HG
 	  AA[p].atom[14].fullatom_type = 66; // Hapo  1HE
 	  AA[p].atom[15].fullatom_type = 66; // Hapo  2HE
 	  AA[p].atom[16].fullatom_type = 66; // Hapo  3HE

 	  AA[p].atom[17].fullatom_type = 57; // OXT
 	  AA[p].atom[18].fullatom_type = 66; // Hpol
 	  AA[p].atom[19].fullatom_type = 66; // Hpol
 	  break;
  	default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 4; // CH2    CG
	  AA[p].atom[6].fullatom_type = 16; // S      SD
	  AA[p].atom[7].fullatom_type = 5; // CH3    CE
	  AA[p].atom[8].fullatom_type = 25; // HNbb   H
	  AA[p].atom[9].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[10].fullatom_type = 23; // Hapo  2HB
	  AA[p].atom[11].fullatom_type = 23; // Hapo  3HB
	  AA[p].atom[12].fullatom_type = 23; // Hapo  2HG
	  AA[p].atom[13].fullatom_type = 23; // Hapo  3HG
	  AA[p].atom[14].fullatom_type = 23; // Hapo  1HE
	  AA[p].atom[15].fullatom_type = 23; // Hapo  2HE
	  AA[p].atom[16].fullatom_type = 23; // Hapo  3HE

	  AA[p].atom[17].fullatom_type = 20; // OXT
	  AA[p].atom[18].fullatom_type = 22; // Hpol
	  AA[p].atom[19].fullatom_type = 22; // Hpol
	  break;
  }

  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 8; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 9; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--CG
  AA[p].atom[4].bonded_neighbor[2] = 10; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 11; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 4; // CG
  AA[p].atom[5].bonded_neighbor[0] = 4; // CG--CB
  AA[p].atom[5].bonded_neighbor[1] = 6; // CG--CD
  AA[p].atom[5].bonded_neighbor[2] = 12; // CG--2HG
  AA[p].atom[5].bonded_neighbor[3] = 13; // CG--3HG
  AA[p].atom[6].nbonded_neighbors = 2; // SD
  AA[p].atom[6].bonded_neighbor[0] = 5; // SD--CG
  AA[p].atom[6].bonded_neighbor[1] = 7; // SD--CE
  AA[p].atom[7].nbonded_neighbors = 4; // CE
  AA[p].atom[7].bonded_neighbor[0] = 6; // CE--SD
  AA[p].atom[7].bonded_neighbor[1] = 14; // CE--1HE
  AA[p].atom[7].bonded_neighbor[2] = 15; // CE--2HE
  AA[p].atom[7].bonded_neighbor[3] = 16; // CE--3HE
  AA[p].atom[8].nbonded_neighbors = 1; // H
  AA[p].atom[8].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[9].nbonded_neighbors = 1; //HA
  AA[p].atom[9].bonded_neighbor[0] = 1; //HA--CA
  AA[p].atom[10].nbonded_neighbors = 1; //2HB
  AA[p].atom[10].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[11].nbonded_neighbors = 1; //3HB
  AA[p].atom[11].bonded_neighbor[0] = 4; //3HB--CB
  AA[p].atom[12].nbonded_neighbors = 1; // 2HG
  AA[p].atom[12].bonded_neighbor[0] = 5; //2HG--CG
  AA[p].atom[13].nbonded_neighbors = 1; //3HG
  AA[p].atom[13].bonded_neighbor[0] = 5; //3HG--CG
  AA[p].atom[14].nbonded_neighbors = 1; //1HE
  AA[p].atom[14].bonded_neighbor[0] = 7; //1HE--CE
  AA[p].atom[15].nbonded_neighbors = 1; //2HE
  AA[p].atom[15].bonded_neighbor[0] = 7; //2HE--CE
  AA[p].atom[16].nbonded_neighbors = 1; //3HE
  AA[p].atom[16].bonded_neighbor[0] = 7; //3HE--CG

  AA[p].atom[17].nbonded_neighbors = 1; // OXT
  AA[p].atom[17].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[18].nbonded_neighbors = 1; // 1H
  AA[p].atom[18].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[19].nbonded_neighbors = 1; // 2H
  AA[p].atom[19].bonded_neighbor[0] = 0; // H-N


  for ( k = 0; k < 9; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 9; k < 17; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template for building  HA
  AA[p].atom[9].ta[0] = 1; //   CA
  AA[p].atom[9].ta[1] = 0; //   N
  AA[p].atom[9].ta[2] = 2; //   C

  //bk   template for building   2HB
  AA[p].atom[10].ta[0] = 4; //   CB
  AA[p].atom[10].ta[1] = 1; //   CA
  AA[p].atom[10].ta[2] = 5; //   CG

  //bk   template for building  3HB
  AA[p].atom[11].ta[0] = 4; //   CB
  AA[p].atom[11].ta[1] = 1; //   CA
  AA[p].atom[11].ta[2] = 5; //   CG

  //bk   template for building  2HG
  AA[p].atom[12].ta[0] = 5; //   CG
  AA[p].atom[12].ta[1] = 4; //   CB
  AA[p].atom[12].ta[2] = 6; //   SD

  //bk   template for building  3HG
  AA[p].atom[13].ta[0] = 5; //   CG
  AA[p].atom[13].ta[1] = 4; //   CB
  AA[p].atom[13].ta[2] = 6; //   SD

  //bk   template for building  1HE
  AA[p].atom[14].ta[0] = 7; //   CE
  AA[p].atom[14].ta[1] = 6; //   SD
  AA[p].atom[14].ta[2] = 5; //   CG

  //bk   template for building  2HE
  AA[p].atom[15].ta[0] = 7; //   CE
  AA[p].atom[15].ta[1] = 6; //   SD
  AA[p].atom[15].ta[2] = 5; //   CG

  //bk   template for building  3HE
  AA[p].atom[16].ta[0] = 7; //   CE
  AA[p].atom[16].ta[1] = 6; //   SD
  AA[p].atom[16].ta[2] = 5; //   CG
  //bk   chi angles required to build atoms MET

  //bk   chi angles needed for building  CG
  AA[p].chi_required[0] [5] = true;

  //bk   chi angles needed for building  SD
  AA[p].chi_required[0] [6] = true;
  AA[p].chi_required[1] [6] = true;


  //bk   chi angles needed for building  CE
  AA[p].chi_required[0] [7] = true;
  AA[p].chi_required[1] [7] = true;
  AA[p].chi_required[2] [7] = true;

  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [10] = true;

  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [11] = true;


  //bk   chi angles needed for building 2HG
  AA[p].chi_required[0] [12] = true;
  AA[p].chi_required[1] [12] = true;


  //bk   chi angles needed for building 3HG
  AA[p].chi_required[0] [13] = true;
  AA[p].chi_required[1] [13] = true;


  //bk   chi angles needed for building 1HE
  AA[p].chi_required[0] [14] = true;
  AA[p].chi_required[1] [14] = true;
  AA[p].chi_required[2] [14] = true;


  //bk   chi angles needed for building 2HE
  AA[p].chi_required[0] [15] = true;
  AA[p].chi_required[1] [15] = true;
  AA[p].chi_required[2] [15] = true;


  //bk   chi angles needed for building 3HE
  AA[p].chi_required[0] [16] = true;
  AA[p].chi_required[1] [16] = true;
  AA[p].chi_required[2] [16] = true;





  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   CG
  //bk   four atoms that define chi angle  2
  AA[p].chi_atoms[1] [0] = 1; //   CA
  AA[p].chi_atoms[1] [1] = 4; //   CB
  AA[p].chi_atoms[1] [2] = 5; //   CG
  AA[p].chi_atoms[1] [3] = 6; //   SD
  //bk   four atoms that define chi angle  3
  AA[p].chi_atoms[2] [0] = 4; //   CB
  AA[p].chi_atoms[2] [1] = 5; //   CG
  AA[p].chi_atoms[2] [2] = 6; //   SD
  AA[p].chi_atoms[2] [3] = 7; //   CE



  AA[p].atom[0].icoor[0] = 0.0000; //   N
  AA[p].atom[0].icoor[1] = 15.5330; //   N
  AA[p].atom[0].icoor[2] = 36.4900; //   N
  AA[p].atom[1].icoor[0] = 0.0000; //   CA
  AA[p].atom[1].icoor[1] = 16.4650; //   CA
  AA[p].atom[1].icoor[2] = 37.6220; //   CA
  AA[p].atom[2].icoor[0] = 0.0000; //   C
  AA[p].atom[2].icoor[1] = 15.5850; //   C
  AA[p].atom[2].icoor[2] = 38.8610; //   C
  AA[p].atom[3].icoor[0] = 0.0630; //   O
  AA[p].atom[3].icoor[1] = 14.3610; //   O
  AA[p].atom[3].icoor[2] = 38.7570; //   O
  AA[p].atom[4].icoor[0] = -1.2080; //   CB
  AA[p].atom[4].icoor[1] = 17.3960; //   CB
  AA[p].atom[4].icoor[2] = 37.5390; //   CB
  AA[p].atom[5].icoor[0] = -1.3030; //   CG
  AA[p].atom[5].icoor[1] = 18.4070; //   CG
  AA[p].atom[5].icoor[2] = 38.6730; //   CG
  AA[p].atom[6].icoor[0] = -2.7450; //   SD
  AA[p].atom[6].icoor[1] = 19.4810; //   SD
  AA[p].atom[6].icoor[2] = 38.5280; //   SD
  AA[p].atom[7].icoor[0] = -2.5470; //   CE
  AA[p].atom[7].icoor[1] = 20.5110; //   CE
  AA[p].atom[7].icoor[2] = 39.9790; //   CE
  AA[p].atom[8].icoor[0] = -0.0400; //   H
  AA[p].atom[8].icoor[1] = 14.5430; //   H
  AA[p].atom[8].icoor[2] = 36.6850; //   H
  AA[p].atom[9].icoor[0] = 0.9070; //   HA
  AA[p].atom[9].icoor[1] = 17.0690; //   HA
  AA[p].atom[9].icoor[2] = 37.6060; //   HA
  AA[p].atom[10].icoor[0] = -1.1380; //   2HB
  AA[p].atom[10].icoor[1] = 17.9220; //   2HB
  AA[p].atom[10].icoor[2] = 36.5880; //   2HB
  AA[p].atom[11].icoor[0] = -2.0970; //   3HB
  AA[p].atom[11].icoor[1] = 16.7640; //   3HB
  AA[p].atom[11].icoor[2] = 37.5410; //   3HB
  AA[p].atom[12].icoor[0] = -1.3580; //   2HG
  AA[p].atom[12].icoor[1] = 17.8590; //   2HG
  AA[p].atom[12].icoor[2] = 39.6130; //   2HG
  AA[p].atom[13].icoor[0] = -0.4000; //   3HG
  AA[p].atom[13].icoor[1] = 19.0170; //   3HG
  AA[p].atom[13].icoor[2] = 38.6600; //   3HG
  AA[p].atom[14].icoor[0] = -3.3650; //   1HE
  AA[p].atom[14].icoor[1] = 21.2300; //   1HE
  AA[p].atom[14].icoor[2] = 40.0300; //   1HE
  AA[p].atom[15].icoor[0] = -2.5560; //   2HE
  AA[p].atom[15].icoor[1] = 19.8860; //   2HE
  AA[p].atom[15].icoor[2] = 40.8730; //   2HE
  AA[p].atom[16].icoor[0] = -1.5980; //   3HE
  AA[p].atom[16].icoor[1] = 21.0450; //   3HE
  AA[p].atom[16].icoor[2] = 39.9190; //   3HE

  // atom number for backbone HN
  AA[p].HNpos=8;
  // atom number for backbone HA
  AA[p].HApos=9;

  // number of polar hydrogens
  AA[p].nH_polar_complete=1;
  // atom numbers for polar H
  AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
  AA[p].Hpos_polar_complete[0]=8;

  // number of aromatic hydrogens
  AA[p].nH_aromatic_complete=0;

  // number of apolar hydrogens
  AA[p].nH_apolar_complete=8;
  //atom number for apolar hydrogens
  AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
  AA[p].Hpos_apolar_complete[0]=9;
  AA[p].Hpos_apolar_complete[1]=10;
  AA[p].Hpos_apolar_complete[2]=11;
  AA[p].Hpos_apolar_complete[3]=12;
  AA[p].Hpos_apolar_complete[4]=13;
  AA[p].Hpos_apolar_complete[5]=14;
  AA[p].Hpos_apolar_complete[6]=15;
  AA[p].Hpos_apolar_complete[7]=16;

  // number of acceptors
  AA[p].nacceptors=1;
  //acceptor information
  AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
  (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;

  //Number of Hydrogens connection
  AA[p].nH_hydrogen_connexions=9;
  //atoms hydrogens are connected too
  AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
  AA[p].Hbase[0][0]=8; AA[p].Hbase[0][1]=0;
  AA[p].Hbase[1][0]=9; AA[p].Hbase[1][1]=1;
  AA[p].Hbase[2][0]=10; AA[p].Hbase[2][1]=4;
  AA[p].Hbase[3][0]=11; AA[p].Hbase[3][1]=4;
  AA[p].Hbase[4][0]=12; AA[p].Hbase[4][1]=5;
  AA[p].Hbase[5][0]=13; AA[p].Hbase[5][1]=5;
  AA[p].Hbase[6][0]=14; AA[p].Hbase[6][1]=7;
  AA[p].Hbase[7][0]=15; AA[p].Hbase[7][1]=7;
  AA[p].Hbase[8][0]=16; AA[p].Hbase[8][1]=7;

  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=8;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=9;

  AA[p].atom[4].numHydrogens_atm=2;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=10;
  AA[p].atom[4].hydrogens_atm[1]=11;

  AA[p].atom[7].numHydrogens_atm=3;
  AA[p].atom[7].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[7].numHydrogens_atm);
  AA[p].atom[7].hydrogens_atm[0]=14;
  AA[p].atom[7].hydrogens_atm[1]=15;
  AA[p].atom[7].hydrogens_atm[2]=16;


  //FIN DE AMINOACIDO METHIONINE
}

/**
* Initializes Aspargine
*
* @param opt: Rosseta or ICM
*/
void init_ASN(Convention opt)
{
  //AMINOACIDO ASPARAGINA
  int k, p = ASN;

  strcpy( AA[p].aa_name3, "ASN" );
  AA[p].aa_name1 = 'N';
  AA[p].mass = 114.106;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = true;
  AA[p].aa_is_nonpolar = false;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = false;

  AA[p].natoms = 14;
  AA[p].nheavyatoms = 8;
  AA[p].nchi = 2;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );
  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;

  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

  AA[p].chi_types = ( int * ) new int[AA[p].nchi];
  AA[p].chi_types[0] = 2;
  AA[p].chi_types[1] = 4;
  AA[p].natoms_EEF1 = 10;


  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " CG " );
  strcpy( AA[p].atom[6].atom_name, " OD1" );
  strcpy( AA[p].atom[7].atom_name, " ND2" );
  strcpy( AA[p].atom[8].atom_name, " H  " );
  strcpy( AA[p].atom[9].atom_name, "1HD2" );
  strcpy( AA[p].atom[10].atom_name, "2HD2" );
  strcpy( AA[p].atom[11].atom_name, " HA " );
  strcpy( AA[p].atom[12].atom_name, "2HB " );
  strcpy( AA[p].atom[13].atom_name, "3HB " );

  strcpy( AA[p].atom[14].atom_name," OXT" );
  strcpy( AA[p].atom[15].atom_name,"2H  " );
  strcpy( AA[p].atom[16].atom_name,"3H  " );

  switch(opt)
  {
    case Rosseta:
    case Sybil:
    AA[p].atom[0].charge = -0.47; //    N
      AA[p].atom[1].charge = 0.07; //    CA
      AA[p].atom[2].charge = 0.51; //    C
      AA[p].atom[3].charge = -0.51; //    O
      AA[p].atom[4].charge = -0.18; //    CB
      AA[p].atom[5].charge = 0.55; //    CG
      AA[p].atom[6].charge = -0.55; //    OD1
      AA[p].atom[7].charge = -0.62; //    ND2
      AA[p].atom[8].charge = 0.31; //    H
      AA[p].atom[9].charge = 0.32; //    1HD2
      AA[p].atom[10].charge = 0.30; //    2HD2
      AA[p].atom[11].charge = 0.09; //    HA
      AA[p].atom[12].charge = 0.09; //    2HB
      AA[p].atom[13].charge = 0.09; //    3HB

      AA[p].atom[14].charge = -0.67; //   0XT
      AA[p].atom[15].charge = 0.33; //  1H
      AA[p].atom[16].charge = 0.33; //  2H
    break;

    case ICM:
      AA[p].atom[0].charge = -0.356; //    N
      AA[p].atom[1].charge = 0.064; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = -0.12; //    CB
      AA[p].atom[5].charge = 0.46; //    CG
      AA[p].atom[6].charge = -0.375; //    OD1
      AA[p].atom[7].charge = -0.445; //    ND2
      AA[p].atom[8].charge = 0.176; //    H
      AA[p].atom[9].charge = 0.2; //    1HD2
      AA[p].atom[10].charge = 0.2; //    2HD2
      AA[p].atom[11].charge = 0.02; //    HA
      AA[p].atom[12].charge = 0.055; //    2HB
      AA[p].atom[13].charge = 0.055; //    3HB

      AA[p].atom[14].charge = -0.67; //   0XT
      AA[p].atom[15].charge = 0.33; //  1H
      AA[p].atom[16].charge = 0.33; //  2H
     break;
   case EEF1:
      AA[p].atom[0].charge = -0.35; //    N
      AA[p].atom[1].charge =  0.10; //    CA
      AA[p].atom[2].charge =  0.55; //    C
      AA[p].atom[3].charge = -0.55; //    O
      AA[p].atom[4].charge =  0.00; //    CB
      AA[p].atom[5].charge =  0.55; //    CG
      AA[p].atom[6].charge = -0.55; //    OD1
      AA[p].atom[7].charge = -0.60; //    ND2
      AA[p].atom[8].charge =  0.25; //    H
      AA[p].atom[9].charge =  0.30; //    1HD2
      AA[p].atom[10].charge = 0.30; //    2HD2
      AA[p].atom[11].charge = 0.00; //    HA
      AA[p].atom[12].charge = 0.00; //    2HB
      AA[p].atom[13].charge = 0.00; //    3HB
      break;
  }
  /*
   1.0 ASN-N:2H    0.33000 HC
   1.0 ASN-N:HA    0.10000 HB
   1.0 ASN-N:3H    0.33000 HC
   1.0 ASN-N:CB   -0.18000 CT2
   1.0 ASN-N:HB1   0.09000 HA
   1.0 ASN-N:1HB   0.09000 HA
   1.0 ASN-N:CG    0.55000 CC
   1.0 ASN-N:OD1  -0.55000 O
   1.0 ASN-N:ND2  -0.62000 NH2
   1.0 ASN-N:1HD2  0.32000 H
   1.0 ASN-N:2HD2  0.30000 H
   1.0 ASN-N:C     0.51000 C
   1.0 ASN-N:O    -0.51000 O
   1.0 ASN-N:N    -0.30000 NH3
   1.0 ASN-N:1H    0.33000 HC
   1.0 ASN-N:CA    0.21000 CT1

   1.0 ASN-C:HA    0.09000 HB
   1.0 ASN-C:CB   -0.18000 CT2
   1.0 ASN-C:HB1   0.09000 HA
   1.0 ASN-C:1HB   0.09000 HA
   1.0 ASN-C:CG    0.55000 CC
   1.0 ASN-C:OD1  -0.55000 O
   1.0 ASN-C:ND2  -0.62000 NH2
   1.0 ASN-C:1HD2  0.32000 H
   1.0 ASN-C:O    -0.67000 OC
   1.0 ASN-C:2HD2  0.30000 H
   1.0 ASN-C:OXT  -0.67000 OC
   1.0 ASN-C:C     0.34000 CC
   1.0 ASN-C:N    -0.47000 NH1
   1.0 ASN-C:HN    0.31000 H
   1.0 ASN-C:CA    0.07000 CT1
*/

  switch (opt) {
  	case EEF1:
	  AA[p].atom[0].fullatom_type =  12; // NH1    N
	  AA[p].atom[1].fullatom_type =   5; // CH1E   CA
	  AA[p].atom[2].fullatom_type =   4; // C      C
	  AA[p].atom[3].fullatom_type =  16; // O      O
	  AA[p].atom[4].fullatom_type =   6; // CH2E   CB
	  AA[p].atom[5].fullatom_type =   4; // C      CG
	  AA[p].atom[6].fullatom_type =  16; // O      OD1
	  AA[p].atom[7].fullatom_type =  13; // NH2    ND2
	  AA[p].atom[8].fullatom_type =   0; // H      H
	  AA[p].atom[9].fullatom_type =   0; // H      1HD2
	  AA[p].atom[10].fullatom_type =  0; // H      2HD2
	  AA[p].atom[11].fullatom_type = -1; // Hapo   HA
	  AA[p].atom[12].fullatom_type = -1; // Hapo   2HB
	  AA[p].atom[13].fullatom_type = -1; // Hapo   3HB

	  AA[p].atom[14].fullatom_type = 17; // OC     OXT
	  AA[p].atom[15].fullatom_type =  1; // HC     Hpol
	  AA[p].atom[16].fullatom_type =  1; // HC     Hpol
	  break;
  	case Sybil:
	  AA[p].atom[0].fullatom_type = 40; // Nbb    N
	  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 10; // CObb   C
	  AA[p].atom[3].fullatom_type = 56; // OCbb   O
	  AA[p].atom[4].fullatom_type = 8; // CH2    CB
	  AA[p].atom[5].fullatom_type = 10; // CNH2   CG
	  AA[p].atom[6].fullatom_type = 56; // ONH2   OD1
	  AA[p].atom[7].fullatom_type = 39; // NH2O   ND2
	  AA[p].atom[8].fullatom_type = 66; // HNbb   H
	  AA[p].atom[9].fullatom_type = 66; // Hpol  1HD2
	  AA[p].atom[10].fullatom_type = 66; // Hpol  2HD2
	  AA[p].atom[11].fullatom_type = 66; // Hapo   HA
	  AA[p].atom[12].fullatom_type = 66; // Hapo  2HB
	  AA[p].atom[13].fullatom_type = 66; // Hapo  3HB

	  AA[p].atom[14].fullatom_type = 57; // OXT
	  AA[p].atom[15].fullatom_type = 66; // Hpol
	  AA[p].atom[16].fullatom_type = 66; // Hpol
	  break;
  	default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 1; // CNH2   CG
	  AA[p].atom[6].fullatom_type = 14; // ONH2   OD1
	  AA[p].atom[7].fullatom_type = 9; // NH2O   ND2
	  AA[p].atom[8].fullatom_type = 25; // HNbb   H
	  AA[p].atom[9].fullatom_type = 22; // Hpol  1HD2
	  AA[p].atom[10].fullatom_type = 22; // Hpol  2HD2
	  AA[p].atom[11].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[12].fullatom_type = 23; // Hapo  2HB
	  AA[p].atom[13].fullatom_type = 23; // Hapo  3HB

	  AA[p].atom[14].fullatom_type = 20; // OXT
	  AA[p].atom[15].fullatom_type = 22; // Hpol
	  AA[p].atom[16].fullatom_type = 22; // Hpol
	  break;
  }

  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 8; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 11; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--CG
  AA[p].atom[4].bonded_neighbor[2] = 12; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 13; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 3; // CG
  AA[p].atom[5].bonded_neighbor[0] = 4; // CG--CB
  AA[p].atom[5].bonded_neighbor[1] = 6; // CG--OD1
  AA[p].atom[5].bonded_neighbor[2] = 7; // CG--ND2
  AA[p].atom[6].nbonded_neighbors = 1; // OD1
  AA[p].atom[6].bonded_neighbor[0] = 5; // OD1--CG
  AA[p].atom[7].nbonded_neighbors = 3; // ND2
  AA[p].atom[7].bonded_neighbor[0] = 5; // ND2--CG
  AA[p].atom[7].bonded_neighbor[0] = 9; // ND2--1HD2
  AA[p].atom[7].bonded_neighbor[0] = 10; // ND2--2HD2
  AA[p].atom[8].nbonded_neighbors = 1; // H
  AA[p].atom[8].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[9].nbonded_neighbors = 1; //1HD2
  AA[p].atom[9].bonded_neighbor[0] = 7; //1HD2--ND2
  AA[p].atom[10].nbonded_neighbors = 1; //2HD2
  AA[p].atom[10].bonded_neighbor[0] = 7; //2HD2--ND2
  AA[p].atom[11].nbonded_neighbors = 1; // HA
  AA[p].atom[11].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[12].nbonded_neighbors = 1; //2HB
  AA[p].atom[12].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[13].nbonded_neighbors = 1; //3HB
  AA[p].atom[13].bonded_neighbor[0] = 4; //3HB--CB

  AA[p].atom[14].nbonded_neighbors = 1; // OXT
  AA[p].atom[14].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[15].nbonded_neighbors = 1; // 1H
  AA[p].atom[15].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[16].nbonded_neighbors = 1; // 2H
  AA[p].atom[16].bonded_neighbor[0] = 0; // H-N

  for ( k = 0; k < 9; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 9; k < 14; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template atoms used for placing atoms ASN

  //bk   template for building 1HD2
  AA[p].atom[9].ta[0] = 7; //   ND2
  AA[p].atom[9].ta[1] = 5; //   CG
  AA[p].atom[9].ta[2] = 4; //   CB

  //bk   template for building 2HD2
  AA[p].atom[10].ta[0] = 7; //   ND2
  AA[p].atom[10].ta[1] = 5; //   CG
  AA[p].atom[10].ta[2] = 4; //   CB

  //bk   template for building  HA
  AA[p].atom[11].ta[0] = 1; //   CA
  AA[p].atom[11].ta[1] = 0; //   N
  AA[p].atom[11].ta[2] = 2; //   C

  //bk   template for building 2HB
  AA[p].atom[12].ta[0] = 4; //   CB
  AA[p].atom[12].ta[1] = 1; //   CA
  AA[p].atom[12].ta[2] = 5; //   CG

  //bk   template for building 3HB
  AA[p].atom[13].ta[0] = 4; //   CB
  AA[p].atom[13].ta[1] = 1; //   CA
  AA[p].atom[13].ta[2] = 5; //   CG


  //bk   chi angles required to build atoms ASN


  //bk   chi angles needed for building  CG
  AA[p].chi_required[0] [5] = true;


  //bk   chi angles needed for building  OD1
  AA[p].chi_required[0] [6] = true;
  AA[p].chi_required[1] [6] = true;


  //bk   chi angles needed for building  ND2
  AA[p].chi_required[0] [7] = true;
  AA[p].chi_required[1] [7] = true;


  //   chi angles needed for building 1HD2
  AA[p].chi_required[0] [9] = true;
  AA[p].chi_required[1] [9] = true;


  //bk   chi angles needed for building 2HD2
  AA[p].chi_required[0] [10] = true;
  AA[p].chi_required[1] [10] = true;


  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [12] = true;


  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [13] = true;



  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   CG

  //bk   four atoms that define chi angle  2
  AA[p].chi_atoms[1] [0] = 1; //   CA
  AA[p].chi_atoms[1] [1] = 4; //   CB
  AA[p].chi_atoms[1] [2] = 5; //   CG
  AA[p].chi_atoms[1] [3] = 6; //   OD1


  //bk   template coordinates for the amino acid
  //jjh new template from Kevin Karplus

  AA[p].atom[0].icoor[0] = 1.463; //   N
  AA[p].atom[0].icoor[1] = 0.000; //   N
  AA[p].atom[0].icoor[2] = 0.000; //   N
  AA[p].atom[1].icoor[0] = 0.000; //   CA
  AA[p].atom[1].icoor[1] = 0.000; //   CA
  AA[p].atom[1].icoor[2] = 0.000; //   CA
  AA[p].atom[2].icoor[0] = -0.565; //   C
  AA[p].atom[2].icoor[1] = 1.411; //   C
  AA[p].atom[2].icoor[2] = 0.000; //   C
  AA[p].atom[3].icoor[0] = -1.619; //   O
  AA[p].atom[3].icoor[1] = 1.666; //   O
  AA[p].atom[3].icoor[2] = -0.603; //   O
  AA[p].atom[4].icoor[0] = -0.555; //   CB
  AA[p].atom[4].icoor[1] = -0.778; //   CB
  AA[p].atom[4].icoor[2] = 1.179; //   CB
  AA[p].atom[5].icoor[0] = -0.282; //   CG
  AA[p].atom[5].icoor[1] = -2.255; //   CG
  AA[p].atom[5].icoor[2] = 1.112; //   CG
  AA[p].atom[6].icoor[0] = 0.050; //   OD1
  AA[p].atom[6].icoor[1] = -2.798; //   OD1
  AA[p].atom[6].icoor[2] = 0.052; //   OD1
  AA[p].atom[7].icoor[0] = -0.499; //   ND2
  AA[p].atom[7].icoor[1] = -2.920; //   ND2
  AA[p].atom[7].icoor[2] = 2.218; //   ND2
  AA[p].atom[8].icoor[0] = 2.425; //   H
  AA[p].atom[8].icoor[1] = -0.306; //   H
  AA[p].atom[8].icoor[2] = 0.000; //   H
  AA[p].atom[9].icoor[0] = -0.338; //   1HD2
  AA[p].atom[9].icoor[1] = -3.907; //   1HD2
  AA[p].atom[9].icoor[2] = 2.247; //   1HD2
  AA[p].atom[10].icoor[0] = -0.825; //   2HD2
  AA[p].atom[10].icoor[1] = -2.442; //   2HD2
  AA[p].atom[10].icoor[2] = 3.033; //   2HD2
  AA[p].atom[11].icoor[0] = -0.370; //   HA
  AA[p].atom[11].icoor[1] = -0.475; //   HA
  AA[p].atom[11].icoor[2] = -0.910; //   HA
  AA[p].atom[12].icoor[0] = -0.397; //   2HB
  AA[p].atom[12].icoor[1] = -0.431; //   2HB
  AA[p].atom[12].icoor[2] = 2.201; //   2HB
  AA[p].atom[13].icoor[0] = -1.591; //   3HB
  AA[p].atom[13].icoor[1] = -0.599; //   3HB
  AA[p].atom[13].icoor[2] = 0.891; //   3HB

  // atom number for backbone HN
  AA[p].HNpos=8;
  // atom number for backbone HA
  AA[p].HApos=11;

  // number of polar hydrogens
  AA[p].nH_polar_complete=3;
  // atom numbers for polar H
  AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
  AA[p].Hpos_polar_complete[0]=8;
  AA[p].Hpos_polar_complete[1]=9;
  AA[p].Hpos_polar_complete[2]=10;

 // number of aromatic hydrogens
 AA[p].nH_aromatic_complete=0;

 // number of apolar hydrogens
 AA[p].nH_apolar_complete=3;
 //atom number for apolar hydrogens
 AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
 AA[p].Hpos_apolar_complete[0]=11;
 AA[p].Hpos_apolar_complete[1]=12;
 AA[p].Hpos_apolar_complete[2]=13;

 // number of acceptors
 AA[p].nacceptors=2;
 //acceptor information
 AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
 AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
 (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
 (AA[p].aBase[1])[0]=6;(AA[p].aBase[1])[1]=5;
 (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;
 (AA[p].accpt_pos[1])[0]=1;(AA[p].accpt_pos[1])[1]=6;

 //Number of Hydrogens connection
 AA[p].nH_hydrogen_connexions=6;
 //atoms hydrogens are connected too
 AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
 AA[p].Hbase[0][0]=8; AA[p].Hbase[0][1]=0;
 AA[p].Hbase[1][0]=9; AA[p].Hbase[1][1]=7;
 AA[p].Hbase[2][0]=10; AA[p].Hbase[2][1]=7;
 AA[p].Hbase[3][0]=11; AA[p].Hbase[3][1]=1;
 AA[p].Hbase[4][0]=12; AA[p].Hbase[4][1]=4;
 AA[p].Hbase[5][0]=13; AA[p].Hbase[5][1]=4;

  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=8;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=11;

  AA[p].atom[4].numHydrogens_atm=2;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=12;
  AA[p].atom[4].hydrogens_atm[1]=13;

  AA[p].atom[7].numHydrogens_atm=2;
  AA[p].atom[7].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[7].numHydrogens_atm);
  AA[p].atom[7].hydrogens_atm[0]=9;
  AA[p].atom[7].hydrogens_atm[1]=10;


}

//FIN DE AMINOACIDO ASPARAGINA

/**
* Initializes Proline
*
* @param opt: Rosseta or ICM
*/
void init_PRO(Convention opt)
{
  //AMINOACIDO PROLINE
  int k, p =PRO;

  strcpy( AA[p].aa_name3, "PRO" );
  AA[p].aa_name1 = 'P';
  AA[p].mass = 97.118;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = false;
  AA[p].aa_is_nonpolar = true;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = false;

  AA[p].natoms = 14;
  AA[p].nheavyatoms = 7;
  //AA[p].nchi = 1;
  AA[p].nchi = 2;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );
  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;

  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

  AA[p].chi_types = ( int * ) new int[AA[p].nchi];
  AA[p].chi_types[0] = 2;
  AA[p].chi_types[1] = 2;
  AA[p].natoms_EEF1 = 6;

  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " CG " );
  strcpy( AA[p].atom[6].atom_name, " CD " );
  strcpy( AA[p].atom[7].atom_name, "2HD " );
  strcpy( AA[p].atom[8].atom_name, "3HD " );
  strcpy( AA[p].atom[9].atom_name, "2HG " );
  strcpy( AA[p].atom[10].atom_name, "3HG " );
  strcpy( AA[p].atom[11].atom_name, "2HB " );
  strcpy( AA[p].atom[12].atom_name, "3HB " );
  strcpy( AA[p].atom[13].atom_name, " HA " );

  strcpy( AA[p].atom[14].atom_name," OXT" );
  strcpy( AA[p].atom[15].atom_name,"2H  " );
  strcpy( AA[p].atom[16].atom_name,"3H  " );

  switch(opt)
  {
    case Rosseta:
    case Sybil:
    AA[p].atom[0].charge = -0.29; //    N
      AA[p].atom[1].charge = 0.02; //    CA
      AA[p].atom[2].charge = 0.51; //    C
      AA[p].atom[3].charge = -0.51; //    O
      AA[p].atom[4].charge = -0.18; //    CB
      AA[p].atom[5].charge = -0.18; //    CG
      AA[p].atom[6].charge = 0.00; //    CD
      AA[p].atom[7].charge = 0.09; //   2HD
      AA[p].atom[8].charge = 0.09; //   3HD
      AA[p].atom[9].charge = 0.09; //   2HG
      AA[p].atom[10].charge = 0.09; //   3HG
      AA[p].atom[11].charge = 0.09; //   2HB
      AA[p].atom[12].charge = 0.09; //   3HB
      AA[p].atom[13].charge = 0.09; //    HA

      AA[p].atom[14].charge = -0.67; //   0XT
      AA[p].atom[15].charge = 0.33; //  1H
      AA[p].atom[16].charge = 0.33; //  2H
      break;

    case ICM:
      AA[p].atom[0].charge = -0.285; //    N
      AA[p].atom[1].charge = 0.05; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.385; //    O
      AA[p].atom[4].charge = -0.025; //    CB
      AA[p].atom[5].charge = -0.05; //    CG
      AA[p].atom[6].charge = 0.1; //    CD
      AA[p].atom[7].charge = 0.01; //   2HD
      AA[p].atom[8].charge = 0.01; //   3HD
      AA[p].atom[9].charge = 0.025; //   2HG
      AA[p].atom[10].charge = 0.025; //   3HG
      AA[p].atom[11].charge = 0.015; //   2HB
      AA[p].atom[12].charge = 0.015; //   3HB
      AA[p].atom[13].charge = 0.04; //    HA

      AA[p].atom[14].charge = -0.67; //   0XT
      AA[p].atom[15].charge = 0.33; //  1H
      AA[p].atom[16].charge = 0.33; //  2H
      break;
   case EEF1:
      AA[p].atom[0].charge = -0.20; //    N
      AA[p].atom[1].charge =  0.10; //    CA
      AA[p].atom[2].charge =  0.55; //    C
      AA[p].atom[3].charge = -0.55; //    O
      AA[p].atom[4].charge =  0.00; //    CB
      AA[p].atom[5].charge =  0.00; //    CG
      AA[p].atom[6].charge =  0.10; //    CD
      AA[p].atom[7].charge =  0.00; //   2HD
      AA[p].atom[8].charge =  0.00; //   3HD
      AA[p].atom[9].charge =  0.00; //   2HG
      AA[p].atom[10].charge = 0.00; //   3HG
      AA[p].atom[11].charge = 0.00; //   2HB
      AA[p].atom[12].charge = 0.00; //   3HB
      AA[p].atom[13].charge = 0.00; //    HA
      break;
  }
  /*
   1.0 PRO-C:1HD   0.09000 HA
   1.0 PRO-C:CA    0.02000 CP1
   1.0 PRO-C:HA    0.09000 HB
   1.0 PRO-C:CB   -0.18000 CP2
   1.0 PRO-C:HB1   0.09000 HA
   1.0 PRO-C:1HB   0.09000 HA
   1.0 PRO-C:CG   -0.18000 CP2
   1.0 PRO-C:HG1   0.09000 HA
   1.0 PRO-C:O    -0.67000 OC
   1.0 PRO-C:1HG   0.09000 HA
   1.0 PRO-C:OXT  -0.67000 OC
   1.0 PRO-C:C     0.34000 CC
   1.0 PRO-C:N    -0.29000 N
   1.0 PRO-C:CD    0.00000 CP3
   1.0 PRO-C:HD1   0.09000 HA
*/

  switch (opt) {
  	case EEF1:
	  AA[p].atom[0].fullatom_type =   9; // N      N
	  AA[p].atom[1].fullatom_type =   5; // CH1E   CA
	  AA[p].atom[2].fullatom_type =   4; // C      C
	  AA[p].atom[3].fullatom_type =  16; // O      O
	  AA[p].atom[4].fullatom_type =   6; // CH2E   CB
	  AA[p].atom[5].fullatom_type =   6; // CH2E   CG
	  AA[p].atom[6].fullatom_type =   6; // CH2E   CD
	  AA[p].atom[7].fullatom_type =  -1; // Hapo   2HD
	  AA[p].atom[8].fullatom_type =  -1; // Hapo   3HD
	  AA[p].atom[9].fullatom_type =  -1; // Hapo   2HG
	  AA[p].atom[10].fullatom_type = -1; // Hapo   3HG
	  AA[p].atom[11].fullatom_type = -1; // Hapo   2HB
	  AA[p].atom[12].fullatom_type = -1; // Hapo   3HB
	  AA[p].atom[13].fullatom_type = -1; // Hapo   HA

	  AA[p].atom[14].fullatom_type = 17; // OC     OXT
	  AA[p].atom[15].fullatom_type =  1; // HC     2H
	  AA[p].atom[16].fullatom_type =  1; // HC     3H
	  break;
  	case Sybil:
	  AA[p].atom[0].fullatom_type = 41; // Nbb    N
	  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 10; // CObb   C
	  AA[p].atom[3].fullatom_type = 56; // OCbb   O
	  AA[p].atom[4].fullatom_type = 8; // CH2    CB
	  AA[p].atom[5].fullatom_type = 8; // CH2   CG
	  AA[p].atom[6].fullatom_type = 5; // CH2    CD
	  AA[p].atom[7].fullatom_type = 66; // Hapo  2HD
	  AA[p].atom[8].fullatom_type = 66; // Hapo  3HD
	  AA[p].atom[9].fullatom_type = 66; // Hapo  2HG
	  AA[p].atom[10].fullatom_type = 66; // Hapo  3HG
	  AA[p].atom[11].fullatom_type = 66; // Hapo  2HB
	  AA[p].atom[12].fullatom_type = 66; // Hapo  3HB
	  AA[p].atom[13].fullatom_type = 66; // Hapo  HA

	  AA[p].atom[14].fullatom_type = 57; // OXT
	  AA[p].atom[15].fullatom_type = 66; // Hpol
	  AA[p].atom[16].fullatom_type = 66; // Hpol
	  break;
  	default:
	  AA[p].atom[0].fullatom_type = 12; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 4; // CH2   CG
	  AA[p].atom[6].fullatom_type = 4; // CH2    CD
	  AA[p].atom[7].fullatom_type = 23; // Hapo  2HD
	  AA[p].atom[8].fullatom_type = 23; // Hapo  3HD
	  AA[p].atom[9].fullatom_type = 23; // Hapo  2HG
	  AA[p].atom[10].fullatom_type = 23; // Hapo  3HG
	  AA[p].atom[11].fullatom_type = 23; // Hapo  2HB
	  AA[p].atom[12].fullatom_type = 23; // Hapo  3HB
	  AA[p].atom[13].fullatom_type = 23; // Hapo  HA

	  AA[p].atom[14].fullatom_type = 20; // OXT
	  AA[p].atom[15].fullatom_type = 22; // Hpol
	  AA[p].atom[16].fullatom_type = 22; // Hpol
	  break;
  }

  //jjh intra residue bonding

  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 6; // N--CD
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 13; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--CG
  AA[p].atom[4].bonded_neighbor[2] = 11; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 12; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 4; // CG
  AA[p].atom[5].bonded_neighbor[0] = 4; // CG--CB
  AA[p].atom[5].bonded_neighbor[1] = 6; // CG--CD
  AA[p].atom[5].bonded_neighbor[2] = 9; // CG--2HG
  AA[p].atom[5].bonded_neighbor[3] = 10; // CG--3HG
  AA[p].atom[6].nbonded_neighbors = 4; // CD
  AA[p].atom[6].bonded_neighbor[0] = 0; // CD--N
  AA[p].atom[6].bonded_neighbor[1] = 5; // CD--CG
  AA[p].atom[6].bonded_neighbor[2] = 7; // CD--2HD
  AA[p].atom[6].bonded_neighbor[3] = 8; // CD--3HD
  AA[p].atom[7].nbonded_neighbors = 1; //2HD
  AA[p].atom[7].bonded_neighbor[0] = 6; //2HD--CD
  AA[p].atom[8].nbonded_neighbors = 1; //3HD
  AA[p].atom[8].bonded_neighbor[0] = 6; //3HD--CD
  AA[p].atom[9].nbonded_neighbors = 1; //2HG
  AA[p].atom[9].bonded_neighbor[0] = 5; //2HG--CG
  AA[p].atom[10].nbonded_neighbors = 1; //3HG
  AA[p].atom[10].bonded_neighbor[0] = 5; //3HG--CG
  AA[p].atom[11].nbonded_neighbors = 1; //2HB
  AA[p].atom[11].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[12].nbonded_neighbors = 1; //3HB
  AA[p].atom[12].bonded_neighbor[0] = 4; //3HB--CB
  AA[p].atom[13].nbonded_neighbors = 1; // HA
  AA[p].atom[13].bonded_neighbor[0] = 1; // HA--CA

  AA[p].atom[14].nbonded_neighbors = 1; // OXT
  AA[p].atom[14].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[15].nbonded_neighbors = 1; // 1H
  AA[p].atom[15].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[16].nbonded_neighbors = 1; // 2H
  AA[p].atom[16].bonded_neighbor[0] = 0; // H-N

  for ( k = 0; k < 7; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 7; k < 14; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template atoms used for placing atoms PRO

  //bk   template for building 2HD
  AA[p].atom[7].ta[0] = 6; //   CD
  AA[p].atom[7].ta[1] = 0; //   N
  AA[p].atom[7].ta[2] = 5; //   CG

  //bk   template for building 3HD
  AA[p].atom[8].ta[0] = 6; //   CD
  AA[p].atom[8].ta[1] = 0; //   N
  AA[p].atom[8].ta[2] = 5; //   CG

  //bk   template for building  2HG
  AA[p].atom[9].ta[0] = 5; //   CG
  AA[p].atom[9].ta[1] = 4; //   CB
  AA[p].atom[9].ta[2] = 6; //   CD

  //bk   template for building 3HG
  AA[p].atom[10].ta[0] = 5; //   CG
  AA[p].atom[10].ta[1] = 4; //   CB
  AA[p].atom[10].ta[2] = 6; //   CD

  //bk   template for building 2HB
  AA[p].atom[11].ta[0] = 4; //   CB
  AA[p].atom[11].ta[1] = 1; //   CA
  AA[p].atom[11].ta[2] = 5; //   CG

  //bk   template for building 3HB
  AA[p].atom[12].ta[0] = 4; //   CB
  AA[p].atom[12].ta[1] = 1; //   CA
  AA[p].atom[12].ta[2] = 5; //   CG

  //bk   template for building HA
  AA[p].atom[13].ta[0] = 1; //   CA
  AA[p].atom[13].ta[1] = 0; //   N
  AA[p].atom[13].ta[2] = 2; //   C

  //bk   chi angles required to build atoms PRO
  //bk   chi angles needed for building  CG
  AA[p].chi_required[0] [5] = true;

  //bk   chi angles needed for building 2HG
  AA[p].chi_required[0] [9] = true;
  AA[p].chi_required[1] [9] = true;

  //bk   chi angles needed for building 3HG
  AA[p].chi_required[0] [10] = true;
  AA[p].chi_required[1] [10] = true;

  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [11] = true;


  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [12] = true;

  //bk   chi angles needed for building CD
  AA[p].chi_required[0] [6] = true;
  AA[p].chi_required[1] [6] = true;

  //bk   chi angles needed for building 1HD
  AA[p].chi_required[0] [7] = true;
  AA[p].chi_required[1] [7] = true;

  //bk   chi angles needed for building 2HD
  AA[p].chi_required[0] [8] = true;
  AA[p].chi_required[1] [8] = true;

  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   CG

  //bk   four atoms that define chi angle  2
  AA[p].chi_atoms[1] [0] = 1; //   CA
  AA[p].chi_atoms[1] [1] = 4; //   CB
  AA[p].chi_atoms[1] [2] = 5; //   CG
  AA[p].chi_atoms[1] [3] = 6; //   CD

  //bk   template coordinates for the amino acid
  //jjh new template from Kevin Karplus

  AA[p].atom[0].icoor[0] = 13.895; //   N
  AA[p].atom[0].icoor[1] = 1.827; //   N
  AA[p].atom[0].icoor[2] = 9.451; //   N
  AA[p].atom[1].icoor[0] = 15.259; //   CA
  AA[p].atom[1].icoor[1] = 2.242; //   CA
  AA[p].atom[1].icoor[2] = 9.192; //   CA
  AA[p].atom[2].icoor[0] = 15.959; //   C
  AA[p].atom[2].icoor[1] = 1.321; //   C
  AA[p].atom[2].icoor[2] = 8.222; //   C
  AA[p].atom[3].icoor[0] = 15.702; //   O
  AA[p].atom[3].icoor[1] = 0.132; //   O
  AA[p].atom[3].icoor[2] = 8.143; //   O
  AA[p].atom[4].icoor[0] = 15.927; //   CB
  AA[p].atom[4].icoor[1] = 2.223; //   CB
  AA[p].atom[4].icoor[2] = 10.577; //   CB
  AA[p].atom[5].icoor[0] = 15.069; //   CG
  AA[p].atom[5].icoor[1] = 1.310; //   CG
  AA[p].atom[5].icoor[2] = 11.409; //   CG
  AA[p].atom[6].icoor[0] = 13.672; //   CD
  AA[p].atom[6].icoor[1] = 1.471; //   CD
  AA[p].atom[6].icoor[2] = 10.862; //   CD
  AA[p].atom[7].icoor[0] = 13.189; //   2HD
  AA[p].atom[7].icoor[1] = 0.650; //   2HD
  AA[p].atom[7].icoor[2] = 10.821; //   2HD
  AA[p].atom[8].icoor[0] = 13.219; //   3HD
  AA[p].atom[8].icoor[1] = 2.235; //   3HD
  AA[p].atom[8].icoor[2] = 11.206; //   3HD
  AA[p].atom[9].icoor[0] = 15.229; //   2HG
  AA[p].atom[9].icoor[1] = 0.510; //   2HG
  AA[p].atom[9].icoor[2] = 11.182; //   2HG
  AA[p].atom[10].icoor[0] = 14.957; //   3HG
  AA[p].atom[10].icoor[1] = 1.715; //   3HG
  AA[p].atom[10].icoor[2] = 12.247; //   3HG
  AA[p].atom[11].icoor[0] = 16.695; //   2HB
  AA[p].atom[11].icoor[1] = 1.763; //   2HB
  AA[p].atom[11].icoor[2] = 10.389; //   2HB
  AA[p].atom[12].icoor[0] = 15.782; //   3HB
  AA[p].atom[12].icoor[1] = 3.034; //   3HB
  AA[p].atom[12].icoor[2] = 10.843; //   3HB
  AA[p].atom[13].icoor[0] = 15.270; //   HA
  AA[p].atom[13].icoor[1] = 3.160; //   HA
  AA[p].atom[13].icoor[2] = 8.805; //   HA


  // atom number for backbone HN
  AA[p].HNpos=-1;
  // atom number for backbone HA
  AA[p].HApos=13;

  // number of polar hydrogens
  AA[p].nH_polar_complete=0;

  // number of aromatic hydrogens
  AA[p].nH_aromatic_complete=0;

  // number of apolar hydrogens
  AA[p].nH_apolar_complete=7;
  //atom number for apolar hydrogens
  AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
  AA[p].Hpos_apolar_complete[0]=7;
  AA[p].Hpos_apolar_complete[1]=8;
  AA[p].Hpos_apolar_complete[2]=9;
  AA[p].Hpos_apolar_complete[3]=10;
  AA[p].Hpos_apolar_complete[4]=11;
  AA[p].Hpos_apolar_complete[5]=12;
  AA[p].Hpos_apolar_complete[6]=13;

  // number of acceptors
  AA[p].nacceptors=1;
  //acceptor information
  AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
  (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;

  //Number of Hydrogens connection
  AA[p].nH_hydrogen_connexions=7;
  //atoms hydrogens are connected too
  AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)* AA[p].nH_hydrogen_connexions);

  AA[p].Hbase[0][0]=7; AA[p].Hbase[0][1]=6;
  AA[p].Hbase[1][0]=8; AA[p].Hbase[1][1]=6;
  AA[p].Hbase[2][0]=9; AA[p].Hbase[2][1]=5;
  AA[p].Hbase[3][0]=10; AA[p].Hbase[3][1]=5;
  AA[p].Hbase[4][0]=11; AA[p].Hbase[4][1]=4;
  AA[p].Hbase[5][0]=12; AA[p].Hbase[5][1]=4;
  AA[p].Hbase[6][0]=13;
  AA[p].Hbase[6][1]=1;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=13;

  AA[p].atom[4].numHydrogens_atm=2;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=11;
  AA[p].atom[4].hydrogens_atm[1]=12;

  AA[p].atom[5].numHydrogens_atm=2;
  AA[p].atom[5].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[5].numHydrogens_atm);
  AA[p].atom[5].hydrogens_atm[0]=9;
  AA[p].atom[5].hydrogens_atm[1]=10;

  AA[p].atom[6].numHydrogens_atm=2;
  AA[p].atom[6].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[6].numHydrogens_atm);
  AA[p].atom[6].hydrogens_atm[0]=7;
  AA[p].atom[6].hydrogens_atm[1]=8;



  //FIN DE AMINOACIDO PROLINE
}
/**
* Initializes Glutamine
*
* @param opt: Rosseta or ICM
*/
void init_GLN(Convention opt)
{
  //AMINOACIDO GLUTAMINA
  int k, p = GLN;

  strcpy( AA[p].aa_name3, "GLN" );
  AA[p].aa_name1 = 'Q';
  AA[p].mass = 128.133;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = true;
  AA[p].aa_is_nonpolar = false;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = false;

  AA[p].natoms = 17;
  AA[p].nheavyatoms = 9;
  AA[p].nchi = 3;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );
  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;

  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

  AA[p].chi_types = ( int * ) new int[AA[p].nchi];
  AA[p].chi_types[0] = 2;
  AA[p].chi_types[1] = 2;
  AA[p].chi_types[2] = 4;
  AA[p].natoms_EEF1 = 11;

  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " CG " );
  strcpy( AA[p].atom[6].atom_name, " CD " );
  strcpy( AA[p].atom[7].atom_name, " OE1" );
  strcpy( AA[p].atom[8].atom_name, " NE2" );
  strcpy( AA[p].atom[9].atom_name, " H  " );
  strcpy( AA[p].atom[10].atom_name, "1HE2" );
  strcpy( AA[p].atom[11].atom_name, "2HE2" );
  strcpy( AA[p].atom[12].atom_name, " HA " );
  strcpy( AA[p].atom[13].atom_name, "2HB " );
  strcpy( AA[p].atom[14].atom_name, "3HB " );
  strcpy( AA[p].atom[15].atom_name, "2HG " );
  strcpy( AA[p].atom[16].atom_name, "3HG " );

  strcpy( AA[p].atom[17].atom_name," OXT" );
  strcpy( AA[p].atom[18].atom_name,"2H  " );
  strcpy( AA[p].atom[19].atom_name,"3H  " );

  switch(opt)
  {
    case Rosseta:
    case Sybil:
    AA[p].atom[0].charge = -0.47; //    N
      AA[p].atom[1].charge = 0.07; //    CA
      AA[p].atom[2].charge = 0.51; //    C
      AA[p].atom[3].charge = -0.51; //    O
      AA[p].atom[4].charge = -0.18; //    CB
      AA[p].atom[5].charge = -0.18; //    CG
      AA[p].atom[6].charge = 0.55; //    CD
      AA[p].atom[7].charge = -0.55; //    OE1
      AA[p].atom[8].charge = -0.62; //    NE2
      AA[p].atom[9].charge = 0.31; //    H
      AA[p].atom[10].charge = 0.32; //   1HE2
      AA[p].atom[11].charge = 0.30; //   2HE2
      AA[p].atom[12].charge = 0.09; //    HA
      AA[p].atom[13].charge = 0.09; //   2HB
      AA[p].atom[14].charge = 0.09; //   3HB
      AA[p].atom[15].charge = 0.09; //   2HG
      AA[p].atom[16].charge = 0.09; //   3HG

      AA[p].atom[17].charge = -0.67; //   0XT
      AA[p].atom[18].charge = 0.33; //  1H
      AA[p].atom[19].charge = 0.33; //  2H
      break;
     case ICM:
      AA[p].atom[0].charge = -0.356; //    N
      AA[p].atom[1].charge = 0.064; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = -0.03; //    CB
      AA[p].atom[5].charge = -0.11; //    CG
      AA[p].atom[6].charge = 0.47; //    CD
      AA[p].atom[7].charge = -0.38; //    OE1
      AA[p].atom[8].charge = -0.43; //    NE2
      AA[p].atom[9].charge = 0.176; //    H
      AA[p].atom[10].charge = 0.185; //   1HE2
      AA[p].atom[11].charge = 0.185; //   2HE2
      AA[p].atom[12].charge = 0.02; //    HA
      AA[p].atom[13].charge = 0.02; //   2HB
      AA[p].atom[14].charge = 0.02; //   3HB
      AA[p].atom[15].charge = 0.05; //   2HG
      AA[p].atom[16].charge = 0.05; //   3HG

      AA[p].atom[17].charge = -0.67; //   0XT
      AA[p].atom[18].charge = 0.33; //  1H
      AA[p].atom[19].charge = 0.33; //  2H
      break;
  case EEF1:
      AA[p].atom[0].charge = -0.35; //    N
      AA[p].atom[1].charge =  0.10; //    CA
      AA[p].atom[2].charge =  0.55; //    C
      AA[p].atom[3].charge = -0.55; //    O
      AA[p].atom[4].charge =  0.00; //    CB
      AA[p].atom[5].charge =  0.00; //    CG
      AA[p].atom[6].charge =  0.55; //    CD
      AA[p].atom[7].charge = -0.55; //    OE1
      AA[p].atom[8].charge = -0.60; //    NE2
      AA[p].atom[9].charge =  0.25; //    H
      AA[p].atom[10].charge = 0.30; //   1HE2
      AA[p].atom[11].charge = 0.30; //   2HE2
      AA[p].atom[12].charge = 0.00; //    HA
      AA[p].atom[13].charge = 0.00; //   2HB
      AA[p].atom[14].charge = 0.00; //   3HB
      AA[p].atom[15].charge = 0.00; //   2HG
      AA[p].atom[16].charge = 0.00; //   3HG
      break;
  }
  /*
   1.0 GLN-N:2H    0.33000 HC
   1.0 GLN-N:O    -0.51000 O
   1.0 GLN-N:HA    0.10000 HB
   1.0 GLN-N:3H    0.33000 HC
   1.0 GLN-N:CB   -0.18000 CT2
   1.0 GLN-N:HB1   0.09000 HA
   1.0 GLN-N:1HB   0.09000 HA
   1.0 GLN-N:CG   -0.18000 CT2
   1.0 GLN-N:HG1   0.09000 HA
   1.0 GLN-N:1HG   0.09000 HA
   1.0 GLN-N:CD    0.55000 CC
   1.0 GLN-N:OE1  -0.55000 O
   1.0 GLN-N:NE2  -0.62000 NH2
   1.0 GLN-N:1HE2  0.32000 H
   1.0 GLN-N:N    -0.30000 NH3
   1.0 GLN-N:2HE2  0.30000 H
   1.0 GLN-N:1H    0.33000 HC
   1.0 GLN-N:C     0.51000 C
   1.0 GLN-N:CA    0.21000 CT1

   1.0 GLN-C:HA    0.09000 HB
   1.0 GLN-C:CB   -0.18000 CT2
   1.0 GLN-C:HB1   0.09000 HA
   1.0 GLN-C:1HB   0.09000 HA
   1.0 GLN-C:CG   -0.18000 CT2
   1.0 GLN-C:HG1   0.09000 HA
   1.0 GLN-C:1HG   0.09000 HA
   1.0 GLN-C:CD    0.55000 CC
   1.0 GLN-C:O    -0.67000 OC
   1.0 GLN-C:OE1  -0.55000 O
   1.0 GLN-C:OXT  -0.67000 OC
   1.0 GLN-C:NE2  -0.62000 NH2
   1.0 GLN-C:1HE2  0.32000 H
   1.0 GLN-C:N    -0.47000 NH1
   1.0 GLN-C:2HE2  0.30000 H
   1.0 GLN-C:HN    0.31000 H
   1.0 GLN-C:C     0.34000 CC
   1.0 GLN-C:CA    0.07000 CT1
*/

  switch (opt) {
  	case EEF1:
	  AA[p].atom[0].fullatom_type =  12; // NH1    N
	  AA[p].atom[1].fullatom_type =   5; // CH1E   CA
	  AA[p].atom[2].fullatom_type =   4; // C      C
	  AA[p].atom[3].fullatom_type =  16; // O      O
	  AA[p].atom[4].fullatom_type =   6; // CH2E   CB
	  AA[p].atom[5].fullatom_type =   6; // CH2E   CG
	  AA[p].atom[6].fullatom_type =   4; // C      CD
	  AA[p].atom[7].fullatom_type =  16; // O      OE1
	  AA[p].atom[8].fullatom_type =  13; // NH2    NE2
	  AA[p].atom[9].fullatom_type =   0; // H      H
	  AA[p].atom[10].fullatom_type =  0; // H      1HE2
	  AA[p].atom[11].fullatom_type =  0; // H      2HE2
	  AA[p].atom[12].fullatom_type = -1; // Hapo   HA
	  AA[p].atom[13].fullatom_type = -1; // Hapo   2HB
	  AA[p].atom[14].fullatom_type = -1; // Hapo   3HB
	  AA[p].atom[15].fullatom_type = -1; // Hapo   2HG
	  AA[p].atom[16].fullatom_type = -1; // Hapo   3HG

	  AA[p].atom[17].fullatom_type = 17; // OC     OXT
	  AA[p].atom[18].fullatom_type =  1; // HC     2H
	  AA[p].atom[19].fullatom_type =  1; // HC     3H
	  break;
  	case Sybil:
	  AA[p].atom[0].fullatom_type = 40; // Nbb    N
	  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 10; // CObb   C
	  AA[p].atom[3].fullatom_type = 56; // OCbb   O
	  AA[p].atom[4].fullatom_type = 8; // CH2    CB
	  AA[p].atom[5].fullatom_type = 8; // CH2    CG
	  AA[p].atom[6].fullatom_type = 10; // CNH2   CD
	  AA[p].atom[7].fullatom_type = 56; // ONH2   OE1
	  AA[p].atom[8].fullatom_type = 39; // NH2O   NE2
	  AA[p].atom[9].fullatom_type = 66; // HNbb   H
	  AA[p].atom[10].fullatom_type = 66; // Hpol  1HE2
	  AA[p].atom[11].fullatom_type = 66; // Hpol  2HE2
	  AA[p].atom[12].fullatom_type = 66; // Hapo   HA
	  AA[p].atom[13].fullatom_type = 66; // Hapo  2HB
	  AA[p].atom[14].fullatom_type = 66; // Hapo  3HB
	  AA[p].atom[15].fullatom_type = 66; // Hapo  2HG
	  AA[p].atom[16].fullatom_type = 66; // Hapo  3HG

	  AA[p].atom[17].fullatom_type = 57; // OXT
	  AA[p].atom[18].fullatom_type = 66; // Hpol
	  AA[p].atom[19].fullatom_type = 66; // Hpol
	  break;
  	default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 4; // CH2    CG
	  AA[p].atom[6].fullatom_type = 1; // CNH2   CD
	  AA[p].atom[7].fullatom_type = 14; // ONH2   OE1
	  AA[p].atom[8].fullatom_type = 9; // NH2O   NE2
	  AA[p].atom[9].fullatom_type = 25; // HNbb   H
	  AA[p].atom[10].fullatom_type = 22; // Hpol  1HE2
	  AA[p].atom[11].fullatom_type = 22; // Hpol  2HE2
	  AA[p].atom[12].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[13].fullatom_type = 23; // Hapo  2HB
	  AA[p].atom[14].fullatom_type = 23; // Hapo  3HB
	  AA[p].atom[15].fullatom_type = 23; // Hapo  2HG
	  AA[p].atom[16].fullatom_type = 23; // Hapo  3HG

	  AA[p].atom[17].fullatom_type = 20; // OXT
	  AA[p].atom[18].fullatom_type = 22; // Hpol
	  AA[p].atom[19].fullatom_type = 22; // Hpol
	  break;
  }

  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 9; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 12; // CA--H1A
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--CG
  AA[p].atom[4].bonded_neighbor[2] = 13; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 14; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 4; // CG
  AA[p].atom[5].bonded_neighbor[0] = 4; // CG--CB
  AA[p].atom[5].bonded_neighbor[1] = 6; // CG--CD
  AA[p].atom[5].bonded_neighbor[2] = 15; // CG--2HG
  AA[p].atom[5].bonded_neighbor[3] = 16; // CG--3HG
  AA[p].atom[6].nbonded_neighbors = 3; // CD
  AA[p].atom[6].bonded_neighbor[0] = 5; // CD--CG
  AA[p].atom[6].bonded_neighbor[1] = 7; // CD--OE1
  AA[p].atom[6].bonded_neighbor[2] = 8; // CD--NE2
  AA[p].atom[7].nbonded_neighbors = 1; // OE1
  AA[p].atom[7].bonded_neighbor[0] = 6; // OE1--CD
  AA[p].atom[8].nbonded_neighbors = 3; // NE2
  AA[p].atom[8].bonded_neighbor[0] = 6; // NE2--CD
  AA[p].atom[8].bonded_neighbor[1] = 10; // NE2--1HE2
  AA[p].atom[8].bonded_neighbor[2] = 11; // NE2--2HE2
  AA[p].atom[9].nbonded_neighbors = 1; // H
  AA[p].atom[9].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[10].nbonded_neighbors = 1; //1HE2
  AA[p].atom[10].bonded_neighbor[0] = 8; //1HE2--NE2
  AA[p].atom[11].nbonded_neighbors = 1; //2HE2
  AA[p].atom[11].bonded_neighbor[0] = 8; //2HE2--NE2
  AA[p].atom[12].nbonded_neighbors = // HA
  AA[p].atom[12].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[13].nbonded_neighbors = 1; //2HB
  AA[p].atom[13].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[14].nbonded_neighbors = 1; //3HB
  AA[p].atom[14].bonded_neighbor[0] = 4; //3HB--CB
  AA[p].atom[15].nbonded_neighbors = 1; //2HG
  AA[p].atom[15].bonded_neighbor[0] = 5; //2HG--CG
  AA[p].atom[16].nbonded_neighbors = 1; //3HG
  AA[p].atom[16].bonded_neighbor[0] = 5; //3HG--CG

  AA[p].atom[17].nbonded_neighbors = 1; // OXT
  AA[p].atom[17].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[18].nbonded_neighbors = 1; // 1H
  AA[p].atom[18].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[19].nbonded_neighbors = 1; // 2H
  AA[p].atom[19].bonded_neighbor[0] = 0; // H-N


  for ( k = 0; k < 10; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 10; k < 17; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template atoms used for placing atoms GLN

  //bk   template for building 1HE2
  AA[p].atom[10].ta[0] = 8; //   NE2
  AA[p].atom[10].ta[1] = 6; //   CD
  AA[p].atom[10].ta[2] = 5; //   CG

  //bk   template for building 2HE2
  AA[p].atom[11].ta[0] = 8; //   NE2
  AA[p].atom[11].ta[1] = 6; //   CD
  AA[p].atom[11].ta[2] = 5; //   CG

  //bk   template for building  HA
  AA[p].atom[12].ta[0] = 1; //   CA
  AA[p].atom[12].ta[1] = 0; //   N
  AA[p].atom[12].ta[2] = 2; //   C

  //bk   template for building 2HB
  AA[p].atom[13].ta[0] = 4; //   CB
  AA[p].atom[13].ta[1] = 1; //   CA
  AA[p].atom[13].ta[2] = 5; //   CG

  //bk   template for building 3HB
  AA[p].atom[14].ta[0] = 4; //   CB
  AA[p].atom[14].ta[1] = 1; //   CA
  AA[p].atom[14].ta[2] = 5; //   CG

  //bk   template for building 2HG
  AA[p].atom[15].ta[0] = 5; //   CG
  AA[p].atom[15].ta[1] = 4; //   CB
  AA[p].atom[15].ta[2] = 6; //   CD

  //bk   template for building 3HG
  AA[p].atom[16].ta[0] = 5; //   CG
  AA[p].atom[16].ta[1] = 4; //   CB
  AA[p].atom[16].ta[2] = 6; //   CD


  //bk   chi angles needed for building  CG
  AA[p].chi_required[0] [5] = true;


  //bk   chi angles needed for building  CD
  AA[p].chi_required[0] [6] = true;
  AA[p].chi_required[1] [6] = true;


  //bk   chi angles needed for building  OE1
  AA[p].chi_required[0] [7] = true;
  AA[p].chi_required[1] [7] = true;
  AA[p].chi_required[2] [7] = true;


  //bk   chi angles needed for building  NE2
  AA[p].chi_required[0] [8] = true;
  AA[p].chi_required[1] [8] = true;
  AA[p].chi_required[2] [8] = true;


  //bk   chi angles needed for building 1HE2
  AA[p].chi_required[0] [10] = true;
  AA[p].chi_required[1] [10] = true;
  AA[p].chi_required[2] [10] = true;


  //bk   chi angles needed for building 2HE2
  AA[p].chi_required[0] [11] = true;
  AA[p].chi_required[1] [11] = true;
  AA[p].chi_required[2] [11] = true;

  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [13] = true;

  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [14] = true;

  //bk   chi angles needed for building 2HG
  AA[p].chi_required[0] [15] = true;
  AA[p].chi_required[1] [15] = true;

  //bk   chi angles needed for building 3HG
  AA[p].chi_required[0] [16] = true;
  AA[p].chi_required[1] [16] = true;


  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   CG

  //bk   four atoms that define chi angle  2
  AA[p].chi_atoms[1] [0] = 1; //   CA
  AA[p].chi_atoms[1] [1] = 4; //   CB
  AA[p].chi_atoms[1] [2] = 5; //   CG
  AA[p].chi_atoms[1] [3] = 6; //   CD

  //bk   four atoms that define chi angle  3
  AA[p].chi_atoms[2] [0] = 4; //   CB
  AA[p].chi_atoms[2] [1] = 5; //   CG
  AA[p].chi_atoms[2] [2] = 6; //   CD
  AA[p].chi_atoms[2] [3] = 7; //   OE1


  //bk   template coordinates for the amino acid
  //jjh new template from Kevin Karplus

  AA[p].atom[0].icoor[0] = 1.458; //   N
  AA[p].atom[0].icoor[1] = 0.000; //   N
  AA[p].atom[0].icoor[2] = 0.000; //   N
  AA[p].atom[1].icoor[0] = 0.000; //   CA
  AA[p].atom[1].icoor[1] = 0.000; //   CA
  AA[p].atom[1].icoor[2] = 0.000; //   CA
  AA[p].atom[2].icoor[0] = -0.553; //   C
  AA[p].atom[2].icoor[1] = 1.420; //   C
  AA[p].atom[2].icoor[2] = 0.000; //   C
  AA[p].atom[3].icoor[0] = -1.607; //   O
  AA[p].atom[3].icoor[1] = 1.682; //   O
  AA[p].atom[3].icoor[2] = -0.580; //   O
  AA[p].atom[4].icoor[0] = -0.536; //   CB
  AA[p].atom[4].icoor[1] = -0.762; //   CB
  AA[p].atom[4].icoor[2] = 1.215; //   CB
  AA[p].atom[5].icoor[0] = -0.317; //   CG
  AA[p].atom[5].icoor[1] = -2.264; //   CG
  AA[p].atom[5].icoor[2] = 1.154; //   CG
  AA[p].atom[6].icoor[0] = -0.856; //   CD
  AA[p].atom[6].icoor[1] = -2.876; //   CD
  AA[p].atom[6].icoor[2] = -0.125; //   CD
  AA[p].atom[7].icoor[0] = -2.024; //   OE1
  AA[p].atom[7].icoor[1] = -2.687; //   OE1
  AA[p].atom[7].icoor[2] = -0.476; //   OE1
  AA[p].atom[8].icoor[0] = -0.008; //   NE2
  AA[p].atom[8].icoor[1] = -3.618; //   NE2
  AA[p].atom[8].icoor[2] = -0.828; //   NE2
  AA[p].atom[9].icoor[0] = 2.424; //   H
  AA[p].atom[9].icoor[1] = -0.294; //   H
  AA[p].atom[9].icoor[2] = 0.000; //   H
  AA[p].atom[10].icoor[0] = -0.307; //   1HE2
  AA[p].atom[10].icoor[1] = -4.048; //   1HE2
  AA[p].atom[10].icoor[2] = -1.681; //   1HE2
  AA[p].atom[11].icoor[0] = 0.930; //   2HE2
  AA[p].atom[11].icoor[1] = -3.747; //   2HE2
  AA[p].atom[11].icoor[2] = -0.506; //   2HE2
  AA[p].atom[12].icoor[0] = -0.365; //   HA
  AA[p].atom[12].icoor[1] = -0.474; //   HA
  AA[p].atom[12].icoor[2] = -0.911; //   HA
  AA[p].atom[13].icoor[0] = -0.034; //   2HB
  AA[p].atom[13].icoor[1] = -0.351; //   2HB
  AA[p].atom[13].icoor[2] = 2.091; //   2HB
  AA[p].atom[14].icoor[0] = -1.603; //   3HB
  AA[p].atom[14].icoor[1] = -0.548; //   3HB
  AA[p].atom[14].icoor[2] = 1.275; //   3HB
  AA[p].atom[15].icoor[0] = 0.645; //   2HG
  AA[p].atom[15].icoor[1] = -2.729; //   2HG
  AA[p].atom[15].icoor[2] = 1.368; //   2HG
  AA[p].atom[16].icoor[0] = -1.014; //   3HG
  AA[p].atom[16].icoor[1] = -2.486; //   3HG
  AA[p].atom[16].icoor[2] = 1.963; //   3HG

  // atom number for backbone HN
  AA[p].HNpos=9;
  // atom number for backbone HA
  AA[p].HApos=12;

  // number of polar hydrogens
  AA[p].nH_polar_complete=3;
  //atom number for apolar hydrogens
  AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
  AA[p].Hpos_polar_complete[0]=9;
  AA[p].Hpos_polar_complete[1]=10;
  AA[p].Hpos_polar_complete[2]=11;

  // number of aromatic hydrogens
  AA[p].nH_aromatic_complete=0;

  // number of apolar hydrogens
  AA[p].nH_apolar_complete=5;
  //atom number for apolar hydrogens
  AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
  AA[p].Hpos_apolar_complete[0]=12;
  AA[p].Hpos_apolar_complete[1]=13;
  AA[p].Hpos_apolar_complete[2]=14;
  AA[p].Hpos_apolar_complete[3]=15;
  AA[p].Hpos_apolar_complete[4]=16;

  // number of acceptors
  AA[p].nacceptors=2;
  //acceptor information
  AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
  (AA[p].aBase[1])[0]=7;(AA[p].aBase[1])[1]=6;
  (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;
  (AA[p].accpt_pos[1])[0]=1;(AA[p].accpt_pos[1])[1]=7;

  //Number of Hydrogens connection
  AA[p].nH_hydrogen_connexions=8;
  //atoms hydrogens are connected too
  AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
  AA[p].Hbase[0][0]=9; AA[p].Hbase[0][1]=0;
  AA[p].Hbase[1][0]=10; AA[p].Hbase[1][1]=8;
  AA[p].Hbase[2][0]=11; AA[p].Hbase[2][1]=8;
  AA[p].Hbase[3][0]=12; AA[p].Hbase[3][1]=1;
  AA[p].Hbase[4][0]=13; AA[p].Hbase[4][1]=4;
  AA[p].Hbase[5][0]=14; AA[p].Hbase[5][1]=4;
  AA[p].Hbase[6][0]=15; AA[p].Hbase[6][1]=5;
  AA[p].Hbase[7][0]=16; AA[p].Hbase[7][1]=5;

  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=9;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=12;

  AA[p].atom[4].numHydrogens_atm=2;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=13;
  AA[p].atom[4].hydrogens_atm[1]=14;

  AA[p].atom[5].numHydrogens_atm=2;
  AA[p].atom[5].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[5].numHydrogens_atm);
  AA[p].atom[5].hydrogens_atm[0]=15;
  AA[p].atom[5].hydrogens_atm[1]=16;

  AA[p].atom[8].numHydrogens_atm=2;
  AA[p].atom[8].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[8].numHydrogens_atm);
  AA[p].atom[8].hydrogens_atm[0]=10;
  AA[p].atom[8].hydrogens_atm[1]=11;


}

//FIN DE AMINOACIDO GLUTAMINA
/**
* Initializes Arginine
*
* @param opt: Rosseta or ICM
*/
void init_ARG(Convention opt)
{
  //AMINOACIDO ARGININA
  int k, p = ARG;

  strcpy( AA[p].aa_name3, "ARG" );
  AA[p].aa_name1 = 'R';
  AA[p].mass = 157.198;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = true;
  AA[p].aa_is_nonpolar = false;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = true;

  AA[p].natoms = 24;
  AA[p].nheavyatoms = 11;
  AA[p].nchi = 4;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );
  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;


  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

  AA[p].chi_types = ( int * ) new int[AA[p].nchi];
  AA[p].chi_types[0] = 2;
  AA[p].chi_types[1] = 2;
  AA[p].chi_types[2] = 2;
  AA[p].chi_types[3] = 3;
  AA[p].natoms_EEF1 = 16;


  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " CG " );
  strcpy( AA[p].atom[6].atom_name, " CD " );
  strcpy( AA[p].atom[7].atom_name, " NE " );
  strcpy( AA[p].atom[8].atom_name, " CZ " );
  strcpy( AA[p].atom[9].atom_name, " NH1" );
  strcpy( AA[p].atom[10].atom_name, " NH2" );
  strcpy( AA[p].atom[11].atom_name, " H  " );
  strcpy( AA[p].atom[12].atom_name, "1HH1" );
  strcpy( AA[p].atom[13].atom_name, "2HH1" );
  strcpy( AA[p].atom[14].atom_name, "1HH2" );
  strcpy( AA[p].atom[15].atom_name, "2HH2" );
  strcpy( AA[p].atom[16].atom_name, " HE " );
  strcpy( AA[p].atom[17].atom_name, " HA " );
  strcpy( AA[p].atom[18].atom_name, "2HB " );
  strcpy( AA[p].atom[19].atom_name, "3HB " );
  strcpy( AA[p].atom[20].atom_name, "2HG " );
  strcpy( AA[p].atom[21].atom_name, "3HG " );
  strcpy( AA[p].atom[22].atom_name, "2HD " );
  strcpy( AA[p].atom[23].atom_name, "3HD " );

  strcpy( AA[p].atom[24].atom_name," OXT" );
  strcpy( AA[p].atom[25].atom_name,"2H  " );
  strcpy( AA[p].atom[26].atom_name,"3H  " );

  switch(opt)
  {
    case Rosseta:
    case Sybil:
      AA[p].atom[0].charge = -0.47; //    N
      AA[p].atom[1].charge = 0.07; //    CA
      AA[p].atom[2].charge = 0.51; //    C
      AA[p].atom[3].charge = -0.51; //    O
      AA[p].atom[4].charge = -0.18; //    CB
      AA[p].atom[5].charge = -0.18; //    CG
      AA[p].atom[6].charge = 0.20; //    CD
      AA[p].atom[7].charge = -0.70; //    NE
      AA[p].atom[8].charge = 0.64; //    CZ
      AA[p].atom[9].charge = -0.80; //    NH1
      AA[p].atom[10].charge = -0.80; //    NH2
      AA[p].atom[11].charge = 0.31; //    H
      AA[p].atom[12].charge = 0.46; //   1HH1
      AA[p].atom[13].charge = 0.46; //   2HH1
      AA[p].atom[14].charge = 0.46; //   1HH2
      AA[p].atom[15].charge = 0.46; //   2HH2
      AA[p].atom[16].charge = 0.44; //    HE
      AA[p].atom[17].charge = 0.09; //    HA
      AA[p].atom[18].charge = 0.09; //   2HB
      AA[p].atom[19].charge = 0.09; //   3HB
      AA[p].atom[20].charge = 0.09; //   2HG
      AA[p].atom[21].charge = 0.09; //   3HG
      AA[p].atom[22].charge = 0.09; //   2HD
      AA[p].atom[23].charge = 0.09; //   3HD

      AA[p].atom[24].charge = -0.67; //   0XT
      AA[p].atom[25].charge = 0.33; //  1H
      AA[p].atom[26].charge = 0.33; //  2H
      break;

    case ICM:
      AA[p].atom[0].charge = -0.356; //    N
      AA[p].atom[1].charge = 0.064; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = -0.03; //    CB
      AA[p].atom[5].charge = -0.03; //    CG
      AA[p].atom[6].charge = 0.12; //    CD
      AA[p].atom[7].charge = -0.3; //    NE
      AA[p].atom[8].charge = 0.58; //    CZ
      AA[p].atom[9].charge = -0.39; //    NH1
      AA[p].atom[10].charge = -0.39; //    NH2
      AA[p].atom[11].charge = 0.176; //    H
      AA[p].atom[12].charge = 0.28; //   1HH1
      AA[p].atom[13].charge = 0.28; //   2HH1
      AA[p].atom[14].charge = 0.28; //   1HH2
      AA[p].atom[15].charge = 0.28; //   2HH2
      AA[p].atom[16].charge = 0.23; //    HE
      AA[p].atom[17].charge = 0.020; //    HA
      AA[p].atom[18].charge = 0.015; //   2HB
      AA[p].atom[19].charge = 0.015; //   3HB
      AA[p].atom[20].charge = 0.03; //   2HG
      AA[p].atom[21].charge = 0.03; //   3HG
      AA[p].atom[22].charge = 0.01; //   2HD
      AA[p].atom[23].charge = 0.01; //   3HD

      AA[p].atom[24].charge = -0.67; //   0XT
      AA[p].atom[25].charge = 0.33; //  1H
      AA[p].atom[26].charge = 0.33; //  2H
      break;
    case EEF1:
      AA[p].atom[0].charge = -0.35; //    N
      AA[p].atom[1].charge =  0.10; //    CA
      AA[p].atom[2].charge =  0.55; //    C
      AA[p].atom[3].charge = -0.55; //    O
      AA[p].atom[4].charge =  0.00; //    CB
      AA[p].atom[5].charge =  0.00; //    CG
      AA[p].atom[6].charge =  0.10; //    CD
      AA[p].atom[7].charge = -0.40; //    NE
      AA[p].atom[8].charge =  0.25; //    CZ
      AA[p].atom[9].charge = -0.85; //    NH1
      AA[p].atom[10].charge =-0.85; //    NH2
      AA[p].atom[11].charge = 0.25; //    H
      AA[p].atom[12].charge = 0.40; //   1HH1
      AA[p].atom[13].charge = 0.40; //   2HH1
      AA[p].atom[14].charge = 0.40; //   1HH2
      AA[p].atom[15].charge = 0.40; //   2HH2
      AA[p].atom[16].charge = 0.15; //    HE
      AA[p].atom[17].charge = 0.00; //    HA
      AA[p].atom[18].charge = 0.00; //   2HB
      AA[p].atom[19].charge = 0.00; //   3HB
      AA[p].atom[20].charge = 0.00; //   2HG
      AA[p].atom[21].charge = 0.00; //   3HG
      AA[p].atom[22].charge = 0.00; //   2HD
      AA[p].atom[23].charge = 0.00; //   3HD
      break;
  }
  /*
   1.0 ARG-N:2H    0.33000 HC
   1.0 ARG-N:NH1  -0.80000 NC2
   1.0 ARG-N:HA    0.10000 HB
   1.0 ARG-N:3H    0.33000 HC
   1.0 ARG-N:1HH1  0.46000 HC
   1.0 ARG-N:CB   -0.18000 CT2
   1.0 ARG-N:2HH1  0.46000 HC
   1.0 ARG-N:HB1   0.09000 HA
   1.0 ARG-N:1HB   0.09000 HA
   1.0 ARG-N:CG   -0.18000 CT2
   1.0 ARG-N:HG1   0.09000 HA
   1.0 ARG-N:1HG   0.09000 HA
   1.0 ARG-N:NH2  -0.80000 NC2
   1.0 ARG-N:CD    0.20000 CT2
   1.0 ARG-N:1HH2  0.46000 HC
   1.0 ARG-N:HD1   0.09000 HA
   1.0 ARG-N:2HH2  0.46000 HC
   1.0 ARG-N:1HD   0.09000 HA
   1.0 ARG-N:C     0.51000 C
   1.0 ARG-N:NE   -0.70000 NC2
   1.0 ARG-N:N    -0.30000 NH3
   1.0 ARG-N:O    -0.51000 O
   1.0 ARG-N:HE    0.44000 HC
   1.0 ARG-N:1H    0.33000 HC
   1.0 ARG-N:CZ    0.64000 C
   1.0 ARG-N:CA    0.21000 CT1

   1.0 ARG-C:NH1  -0.80000 NC2
   1.0 ARG-C:HA    0.09000 HB
   1.0 ARG-C:1HH1  0.46000 HC
   1.0 ARG-C:CB   -0.18000 CT2
   1.0 ARG-C:2HH1  0.46000 HC
   1.0 ARG-C:HB1   0.09000 HA
   1.0 ARG-C:1HB   0.09000 HA
   1.0 ARG-C:CG   -0.18000 CT2
   1.0 ARG-C:HG1   0.09000 HA
   1.0 ARG-C:1HG   0.09000 HA
   1.0 ARG-C:NH2  -0.80000 NC2
   1.0 ARG-C:CD    0.20000 CT2
   1.0 ARG-C:O    -0.67000 OC
   1.0 ARG-C:1HH2  0.46000 HC
   1.0 ARG-C:HD1   0.09000 HA
   1.0 ARG-C:OXT  -0.67000 OC
   1.0 ARG-C:2HH2  0.46000 HC
   1.0 ARG-C:1HD   0.09000 HA
   1.0 ARG-C:C     0.34000 CC
   1.0 ARG-C:NE   -0.70000 NC2
   1.0 ARG-C:N    -0.47000 NH1
   1.0 ARG-C:HE    0.44000 HC
   1.0 ARG-C:HN    0.31000 H
   1.0 ARG-C:CZ    0.64000 C
   1.0 ARG-C:CA    0.07000 CT1
*/

  switch (opt) {
  	case EEF1:
	  AA[p].atom[0].fullatom_type =  12; // NH1    N
	  AA[p].atom[1].fullatom_type =   5; // CH1E   CA
	  AA[p].atom[2].fullatom_type =   4; // C      C
	  AA[p].atom[3].fullatom_type =  16; // O      O
	  AA[p].atom[4].fullatom_type =   6; // CH2E   CB
	  AA[p].atom[5].fullatom_type =   6; // CH2E   CG
	  AA[p].atom[6].fullatom_type =   6; // CH2E   CD
	  AA[p].atom[7].fullatom_type =  12; // NH1    NE
	  AA[p].atom[8].fullatom_type =  24; // CR     CZ
	  AA[p].atom[9].fullatom_type =  15; // NC2    NH1
	  AA[p].atom[10].fullatom_type = 15; // NC2    NH2
	  AA[p].atom[11].fullatom_type =  0; // H      H
	  AA[p].atom[12].fullatom_type =  1; // HC     1HH1
	  AA[p].atom[13].fullatom_type =  1; // HC     2HH1
	  AA[p].atom[14].fullatom_type =  1; // HC     1HH2
	  AA[p].atom[15].fullatom_type =  1; // HC     2HH2
	  AA[p].atom[16].fullatom_type =  0; // H      HE
	  AA[p].atom[17].fullatom_type = -1; // Hapo   HA
	  AA[p].atom[18].fullatom_type = -1; // Hapo   2HB
	  AA[p].atom[19].fullatom_type = -1; // Hapo   3HB
	  AA[p].atom[20].fullatom_type = -1; // Hapo   2HG
	  AA[p].atom[21].fullatom_type = -1; // Hapo   3HG
	  AA[p].atom[22].fullatom_type = -1; // Hapo   2HD
	  AA[p].atom[23].fullatom_type = -1; // Hapo   3HD

	  AA[p].atom[24].fullatom_type = 17; // OC     OXT
	  AA[p].atom[25].fullatom_type =  1; // HC     Hpol
	  AA[p].atom[26].fullatom_type =  1; // HC     Hpol
	  break;
  	case Sybil:
	  AA[p].atom[0].fullatom_type = 40; // Nbb    N
	  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 10; // CObb   C
	  AA[p].atom[3].fullatom_type = 56; // OCbb   O
	  AA[p].atom[4].fullatom_type = 8; // CH2    CB
	  AA[p].atom[5].fullatom_type = 8; // CH2    CG
	  AA[p].atom[6].fullatom_type = 10; // CH2    CD
	  AA[p].atom[7].fullatom_type = 47; // Narg   NE
	  AA[p].atom[8].fullatom_type = 16; // aroC   CZ
	  AA[p].atom[9].fullatom_type = 47; // Narg   NH1
	  AA[p].atom[10].fullatom_type = 47; // Narg   NH2
	  AA[p].atom[11].fullatom_type = 66; // HNbb   H
	  AA[p].atom[12].fullatom_type = 66; // Hpol  1HH1
	  AA[p].atom[13].fullatom_type = 66; // Hpol  2HH1
	  AA[p].atom[14].fullatom_type = 66; // Hpol  1HH2
	  AA[p].atom[15].fullatom_type = 66; // Hpol  2HH2
	  AA[p].atom[16].fullatom_type = 66; // Hpol   HE
	  AA[p].atom[17].fullatom_type = 66; // Hapo   HA
	  AA[p].atom[18].fullatom_type = 66; // Hapo  2HB
	  AA[p].atom[19].fullatom_type = 66; // Hapo  3HB
	  AA[p].atom[20].fullatom_type = 66; // Hapo  2HG
	  AA[p].atom[21].fullatom_type = 66; // Hapo  3HG
	  AA[p].atom[22].fullatom_type = 66; // Hapo  2HD
	  AA[p].atom[23].fullatom_type = 66; // Hapo  3HD

	  AA[p].atom[24].fullatom_type = 57; // OXT
	  AA[p].atom[25].fullatom_type = 66; // Hpol
	  AA[p].atom[26].fullatom_type = 66; // Hpol
	  break;
  	default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 4; // CH2    CG
	  AA[p].atom[6].fullatom_type = 4; // CH2    CD
	  AA[p].atom[7].fullatom_type = 11; // Narg   NE
	  AA[p].atom[8].fullatom_type = 6; // aroC   CZ
	  AA[p].atom[9].fullatom_type = 11; // Narg   NH1
	  AA[p].atom[10].fullatom_type = 11; // Narg   NH2
	  AA[p].atom[11].fullatom_type = 25; // HNbb   H
	  AA[p].atom[12].fullatom_type = 22; // Hpol  1HH1
	  AA[p].atom[13].fullatom_type = 22; // Hpol  2HH1
	  AA[p].atom[14].fullatom_type = 22; // Hpol  1HH2
	  AA[p].atom[15].fullatom_type = 22; // Hpol  2HH2
	  AA[p].atom[16].fullatom_type = 22; // Hpol   HE
	  AA[p].atom[17].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[18].fullatom_type = 23; // Hapo  2HB
	  AA[p].atom[19].fullatom_type = 23; // Hapo  3HB
	  AA[p].atom[20].fullatom_type = 23; // Hapo  2HG
	  AA[p].atom[21].fullatom_type = 23; // Hapo  3HG
	  AA[p].atom[22].fullatom_type = 23; // Hapo  2HD
	  AA[p].atom[23].fullatom_type = 23; // Hapo  3HD

	  AA[p].atom[24].fullatom_type = 20; // OXT
	  AA[p].atom[25].fullatom_type = 22; // Hpol
	  AA[p].atom[26].fullatom_type = 22; // Hpol
	  break;
  }

  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 11; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 17; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--CG
  AA[p].atom[4].bonded_neighbor[2] = 18; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 19; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 4; // CG
  AA[p].atom[5].bonded_neighbor[0] = 4; // CG--CB
  AA[p].atom[5].bonded_neighbor[1] = 6; // CG--CD
  AA[p].atom[5].bonded_neighbor[2] = 20; // CG--2HG
  AA[p].atom[5].bonded_neighbor[3] = 21; // CG--3HG
  AA[p].atom[6].nbonded_neighbors = 4; // CD
  AA[p].atom[6].bonded_neighbor[0] = 5; // CD--CG
  AA[p].atom[6].bonded_neighbor[1] = 7; // CD--NE
  AA[p].atom[6].bonded_neighbor[2] = 22; // CD--2HD
  AA[p].atom[6].bonded_neighbor[3] = 23; // CD--3HD
  AA[p].atom[7].nbonded_neighbors = 3; // NE
  AA[p].atom[7].bonded_neighbor[0] = 6; // NE--CD
  AA[p].atom[7].bonded_neighbor[1] = 8; // NE--CZ
  AA[p].atom[7].bonded_neighbor[2] = 16; // NE--HE
  AA[p].atom[8].nbonded_neighbors = 3; // CZ
  AA[p].atom[8].bonded_neighbor[0] = 7; // CZ--NE
  AA[p].atom[8].bonded_neighbor[1] = 9; // CZ--NH1
  AA[p].atom[8].bonded_neighbor[2] = 10; // CZ--NH2
  AA[p].atom[9].nbonded_neighbors = 3; // NH1
  AA[p].atom[9].bonded_neighbor[0] = 8; // NH1--CZ
  AA[p].atom[9].bonded_neighbor[1] = 12; // NH1--1HH1
  AA[p].atom[9].bonded_neighbor[2] = 13; // NH1--2HH1
  AA[p].atom[10].nbonded_neighbors = 3; // NH2
  AA[p].atom[10].bonded_neighbor[0] = 8; // NH2--CZ
  AA[p].atom[10].bonded_neighbor[1] = 14; // NH2--1HH2
  AA[p].atom[10].bonded_neighbor[2] = 15; // NH2--2HH2
  AA[p].atom[11].nbonded_neighbors = 1; // H
  AA[p].atom[11].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[12].nbonded_neighbors = 1; //1HH1
  AA[p].atom[12].bonded_neighbor[0] = 9; //1HH1--NH1
  AA[p].atom[13].nbonded_neighbors = 1; //2HH1
  AA[p].atom[13].bonded_neighbor[0] = 9; //2HH1--NH1
  AA[p].atom[14].nbonded_neighbors = 1; //1HH2
  AA[p].atom[14].bonded_neighbor[0] = 10; //1HH2--NH2
  AA[p].atom[15].nbonded_neighbors = 1; //2HH2
  AA[p].atom[15].bonded_neighbor[0] = 10; //2HH2--NH2
  AA[p].atom[16].nbonded_neighbors = 1; // HE
  AA[p].atom[16].bonded_neighbor[0] = 7; // HE--NE
  AA[p].atom[17].nbonded_neighbors = 1; // HA
  AA[p].atom[17].bonded_neighbor[0] = 1; //HA--C
  AA[p].atom[18].nbonded_neighbors = 1; //2HB
  AA[p].atom[18].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[19].nbonded_neighbors = 1; //3HB
  AA[p].atom[19].bonded_neighbor[0] = 4; //3HB--CB
  AA[p].atom[20].nbonded_neighbors = 1; //2HG
  AA[p].atom[20].bonded_neighbor[0] = 5; //2HG--CG
  AA[p].atom[21].nbonded_neighbors = 1; //3HG
  AA[p].atom[21].bonded_neighbor[0] = 5; //3HG--CG
  AA[p].atom[22].nbonded_neighbors = 1; //2HD
  AA[p].atom[22].bonded_neighbor[0] = 6; //2HD--CD
  AA[p].atom[23].nbonded_neighbors = 1; //3HD
  AA[p].atom[23].bonded_neighbor[0] = 6; //3HD--CD

  AA[p].atom[24].nbonded_neighbors = 1; // OXT
  AA[p].atom[24].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[25].nbonded_neighbors = 1; // 1H
  AA[p].atom[25].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[26].nbonded_neighbors = 1; // 2H
  AA[p].atom[26].bonded_neighbor[0] = 0; // H-N


  for ( k = 0; k < 12; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 12; k < 24; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template atoms used for placing atoms ARG

  //bk   template for building 1HH1
  AA[p].atom[12].ta[0] = 9; //   NH1
  AA[p].atom[12].ta[1] = 8; //   CZ
  AA[p].atom[12].ta[2] = 7; //   NE

  //bk   template for building 2HH1
  AA[p].atom[13].ta[0] = 9; //   NH1
  AA[p].atom[13].ta[1] = 8; //   CZ
  AA[p].atom[13].ta[2] = 7; //   NE

  //bk   template for building  1HH2
  AA[p].atom[14].ta[0] = 10; //   NH2
  AA[p].atom[14].ta[1] = 8; //   CZ
  AA[p].atom[14].ta[2] = 7; //   NE

  //bk   template for building 2HH2
  AA[p].atom[15].ta[0] = 10; //   NH2
  AA[p].atom[15].ta[1] = 8; //   CZ
  AA[p].atom[15].ta[2] = 7; //   NE

  //bk   template for building HE
  AA[p].atom[16].ta[0] = 7; //   NE
  AA[p].atom[16].ta[1] = 6; //   CD
  AA[p].atom[16].ta[2] = 8; //   CZ

  //bk   template for building HA
  AA[p].atom[17].ta[0] = 1; //   CA
  AA[p].atom[17].ta[1] = 0; //   N
  AA[p].atom[17].ta[2] = 2; //   C

  //bk   template for building 2HB
  AA[p].atom[18].ta[0] = 4; //   CB
  AA[p].atom[18].ta[1] = 1; //   CA
  AA[p].atom[18].ta[2] = 5; //   CG

  //bk   template for building  3HB
  AA[p].atom[19].ta[0] = 4; //   CB
  AA[p].atom[19].ta[1] = 1; //   CA
  AA[p].atom[19].ta[2] = 5; //   CG

  //bk   template for building 2HG
  AA[p].atom[20].ta[0] = 5; //   CG
  AA[p].atom[20].ta[1] = 4; //   CB
  AA[p].atom[20].ta[2] = 6; //   CD

  //bk   template for building 3HG
  AA[p].atom[21].ta[0] = 5; //   CG
  AA[p].atom[21].ta[1] = 4; //   CB
  AA[p].atom[21].ta[2] = 6; //   CD

  //bk   template for building 2HD
  AA[p].atom[22].ta[0] = 6; //   CD
  AA[p].atom[22].ta[1] = 5; //   CG
  AA[p].atom[22].ta[2] = 7; //   NE

  //bk   template for building 3HD
  AA[p].atom[23].ta[0] = 6; //   CD
  AA[p].atom[23].ta[1] = 5; //   CG
  AA[p].atom[23].ta[2] = 7; //   NE


  // ARG
  //bk   chi angles needed for building  N
  //bk   chi angles needed for building  CA
  //bk   chi angles needed for building  C
  //bk   chi angles needed for building  O
  //bk   chi angles needed for building  CB
  //bk   chi angles needed for building  CG
  AA[p].chi_required[0] [5] = true;

  //bk   chi angles needed for building  CD
  AA[p].chi_required[0] [6] = true;
  AA[p].chi_required[1] [6] = true;


  //bk   chi angles needed for building  NE
  AA[p].chi_required[0] [7] = true;
  AA[p].chi_required[1] [7] = true;
  AA[p].chi_required[2] [7] = true;

  //bk   chi angles needed for building  CZ
  AA[p].chi_required[0] [8] = true;
  AA[p].chi_required[1] [8] = true;
  AA[p].chi_required[2] [8] = true;
  AA[p].chi_required[3] [8] = true;

  //bk   chi angles needed for building  NH1
  AA[p].chi_required[0] [9] = true;
  AA[p].chi_required[1] [9] = true;
  AA[p].chi_required[2] [9] = true;
  AA[p].chi_required[3] [9] = true;

  //bk   chi angles needed for building  NH2
  AA[p].chi_required[0] [10] = true;
  AA[p].chi_required[1] [10] = true;
  AA[p].chi_required[2] [10] = true;
  AA[p].chi_required[3] [10] = true;

  //bk   chi angles needed for building  H

  //bk   chi angles needed for building 1HH1
  AA[p].chi_required[0] [12] = true;
  AA[p].chi_required[1] [12] = true;
  AA[p].chi_required[2] [12] = true;
  AA[p].chi_required[3] [12] = true;

  //bk   chi angles needed for building 2HH1
  AA[p].chi_required[0] [13] = true;
  AA[p].chi_required[1] [13] = true;
  AA[p].chi_required[2] [13] = true;
  AA[p].chi_required[3] [13] = true;

  //bk   chi angles needed for building 1HH2
  AA[p].chi_required[0] [14] = true;
  AA[p].chi_required[1] [14] = true;
  AA[p].chi_required[2] [14] = true;
  AA[p].chi_required[3] [14] = true;

  //bk   chi angles needed for building 2HH2
  AA[p].chi_required[0] [15] = true;
  AA[p].chi_required[1] [15] = true;
  AA[p].chi_required[2] [15] = true;
  AA[p].chi_required[3] [15] = true;

  //bk   chi angles needed for building  HE
  AA[p].chi_required[0] [16] = true;
  AA[p].chi_required[1] [16] = true;
  AA[p].chi_required[2] [16] = true;
  AA[p].chi_required[3] [16] = true;


  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [18] = true;

  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [19] = true;


  //bk   chi angles needed for building 2HG
  AA[p].chi_required[0] [20] = true;
  AA[p].chi_required[1] [20] = true;

  //bk   chi angles needed for building 3HG
  AA[p].chi_required[0] [21] = true;
  AA[p].chi_required[1] [21] = true;

  //bk   chi angles needed for building 2HD
  AA[p].chi_required[0] [22] = true;
  AA[p].chi_required[1] [22] = true;
  AA[p].chi_required[2] [22] = true;

  //bk   chi angles needed for building 3HD
  AA[p].chi_required[0] [23] = true;
  AA[p].chi_required[1] [23] = true;
  AA[p].chi_required[2] [23] = true;




  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   CG

  //bk   four atoms that define chi angle  2
  AA[p].chi_atoms[1] [0] = 1; //   CA
  AA[p].chi_atoms[1] [1] = 4; //   CB
  AA[p].chi_atoms[1] [2] = 5; //   CG
  AA[p].chi_atoms[1] [3] = 6; //   CD

  //bk   four atoms that define chi angle  3
  AA[p].chi_atoms[2] [0] = 4; //   CB
  AA[p].chi_atoms[2] [1] = 5; //   CG
  AA[p].chi_atoms[2] [2] = 6; //   CD
  AA[p].chi_atoms[2] [3] = 7; //   NE

  //bk   four atoms that define chi angle  4
  AA[p].chi_atoms[3] [0] = 5; //   CG
  AA[p].chi_atoms[3] [1] = 6; //   CD
  AA[p].chi_atoms[3] [2] = 7; //   NE
  AA[p].chi_atoms[3] [3] = 8; //   CZ


  //bk   template coordinates for the amino acid
  //jjh new template from Kevin Karplus

  AA[p].atom[0].icoor[0] = 1.455; //  N
  AA[p].atom[0].icoor[1] = 0.000; //  N
  AA[p].atom[0].icoor[2] = 0.000; //  N
  AA[p].atom[1].icoor[0] = 0.000; //  CA
  AA[p].atom[1].icoor[1] = 0.000; //  CA
  AA[p].atom[1].icoor[2] = 0.000; //  CA
  AA[p].atom[2].icoor[0] = -0.569; //  C
  AA[p].atom[2].icoor[1] = 1.399; //  C
  AA[p].atom[2].icoor[2] = -0.000; //  C
  AA[p].atom[3].icoor[0] = -1.568; //  O
  AA[p].atom[3].icoor[1] = 1.674; //  O
  AA[p].atom[3].icoor[2] = -0.695; //  O
  AA[p].atom[4].icoor[0] = -0.574; //  CB
  AA[p].atom[4].icoor[1] = -0.820; //  CB
  AA[p].atom[4].icoor[2] = 1.146; //  CB
  AA[p].atom[5].icoor[0] = -0.335; //  CG
  AA[p].atom[5].icoor[1] = -2.318; //  CG
  AA[p].atom[5].icoor[2] = 1.043; //  CG
  AA[p].atom[6].icoor[0] = -1.198; //  CD
  AA[p].atom[6].icoor[1] = -3.009; //  CD
  AA[p].atom[6].icoor[2] = 0.051; //  CD
  AA[p].atom[7].icoor[0] = -2.624; //  NE
  AA[p].atom[7].icoor[1] = -2.901; //  NE
  AA[p].atom[7].icoor[2] = 0.314; //  NE
  AA[p].atom[8].icoor[0] = -3.279; //  CZ
  AA[p].atom[8].icoor[1] = -3.584; //  CZ
  AA[p].atom[8].icoor[2] = 1.273; //  CZ
  AA[p].atom[9].icoor[0] = -2.651; //  NH1
  AA[p].atom[9].icoor[1] = -4.451; //  NH1
  AA[p].atom[9].icoor[2] = 2.036; //  NH1
  AA[p].atom[10].icoor[0] = -4.577; //  NH2
  AA[p].atom[10].icoor[1] = -3.377; //  NH2
  AA[p].atom[10].icoor[2] = 1.411; //  NH2
  AA[p].atom[11].icoor[0] = 2.418; //  H
  AA[p].atom[11].icoor[1] = -0.306; //  H
  AA[p].atom[11].icoor[2] = 0.000; //  H
  AA[p].atom[12].icoor[0] = -1.662; //  1HH1
  AA[p].atom[12].icoor[1] = -4.611; //  1HH1
  AA[p].atom[12].icoor[2] = 1.906; //  1HH1
  AA[p].atom[13].icoor[0] = -3.160; //  2HH1
  AA[p].atom[13].icoor[1] = -4.953; //  2HH1
  AA[p].atom[13].icoor[2] = 2.749; //  2HH1
  AA[p].atom[14].icoor[0] = -5.048; //  1HH2
  AA[p].atom[14].icoor[1] = -2.720; //  1HH2
  AA[p].atom[14].icoor[2] = 0.804; //  1HH2
  AA[p].atom[15].icoor[0] = -5.092; //  2HH2
  AA[p].atom[15].icoor[1] = -3.876; //  2HH2
  AA[p].atom[15].icoor[2] = 2.121; //  2HH2
  AA[p].atom[16].icoor[0] = -3.323; //  HE
  AA[p].atom[16].icoor[1] = -2.333; //  HE
  AA[p].atom[16].icoor[2] = -0.146; //  HE
  AA[p].atom[17].icoor[0] = -0.369; //  HA
  AA[p].atom[17].icoor[1] = -0.476; //  HA
  AA[p].atom[17].icoor[2] = -0.909; //  HA
  AA[p].atom[18].icoor[0] = -0.121; //  2HB
  AA[p].atom[18].icoor[1] = -0.446; //  2HB
  AA[p].atom[18].icoor[2] = 2.063; //  2HB
  AA[p].atom[19].icoor[0] = -1.647; //  3HB
  AA[p].atom[19].icoor[1] = -0.629; //  3HB
  AA[p].atom[19].icoor[2] = 1.167; //  3HB
  AA[p].atom[20].icoor[0] = 0.704; //  2HG
  AA[p].atom[20].icoor[1] = -2.485; //  2HG
  AA[p].atom[20].icoor[2] = 0.758; //  2HG
  AA[p].atom[21].icoor[0] = -0.520; //  3HG
  AA[p].atom[21].icoor[1] = -2.766; //  3HG
  AA[p].atom[21].icoor[2] = 2.020; //  3HG
  AA[p].atom[22].icoor[0] = -1.015; //  2HD
  AA[p].atom[22].icoor[1] = -2.583; //  2HD
  AA[p].atom[22].icoor[2] = -0.935; //  2HD
  AA[p].atom[23].icoor[0] = -0.947; //  3HD
  AA[p].atom[23].icoor[1] = -4.069; //  3HD
  AA[p].atom[23].icoor[2] = 0.041; //  3HD

  // atom number for backbone HN
  AA[p].HNpos=11;
  // atom number for backbone HA
  AA[p].HApos=17;

  // number of polar hydrogens
  AA[p].nH_polar_complete=6;
  //atom number for apolar hydrogens
  AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
  AA[p].Hpos_polar_complete[0]=11;
  AA[p].Hpos_polar_complete[1]=12;
  AA[p].Hpos_polar_complete[2]=13;
  AA[p].Hpos_polar_complete[3]=14;
  AA[p].Hpos_polar_complete[4]=15;
  AA[p].Hpos_polar_complete[5]=16;

  // number of aromatic hydrogens
  AA[p].nH_aromatic_complete=0;

  // number of apolar hydrogens
  AA[p].nH_apolar_complete=7;
  //atom number for apolar hydrogens
  AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
  AA[p].Hpos_apolar_complete[0]=17;
  AA[p].Hpos_apolar_complete[1]=18;
  AA[p].Hpos_apolar_complete[2]=19;
  AA[p].Hpos_apolar_complete[3]=20;
  AA[p].Hpos_apolar_complete[4]=21;
  AA[p].Hpos_apolar_complete[5]=22;
  AA[p].Hpos_apolar_complete[6]=23;

  // number of acceptors
  AA[p].nacceptors=1;
  //acceptor information
  AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
  (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;

  //Number of Hydrogens connection
  AA[p].nH_hydrogen_connexions=13;
  //atoms hydrogens are connected too
  AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
  AA[p].Hbase[0][0]=11; AA[p].Hbase[0][1]=0;
  AA[p].Hbase[1][0]=12; AA[p].Hbase[1][1]=9;
  AA[p].Hbase[2][0]=13; AA[p].Hbase[2][1]=9;
  AA[p].Hbase[3][0]=14; AA[p].Hbase[3][1]=10;
  AA[p].Hbase[4][0]=15; AA[p].Hbase[4][1]=10;
  AA[p].Hbase[5][0]=16; AA[p].Hbase[5][1]=7;
  AA[p].Hbase[6][0]=17; AA[p].Hbase[6][1]=1;
  AA[p].Hbase[7][0]=18; AA[p].Hbase[7][1]=4;
  AA[p].Hbase[8][0]=19; AA[p].Hbase[8][1]=4;
  AA[p].Hbase[9][0]=20; AA[p].Hbase[9][1]=5;
  AA[p].Hbase[10][0]=21; AA[p].Hbase[10][1]=5;
  AA[p].Hbase[11][0]=22; AA[p].Hbase[11][1]=6;
  AA[p].Hbase[12][0]=23; AA[p].Hbase[12][1]=6;

  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=11;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=17;

  AA[p].atom[4].numHydrogens_atm=2;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=18;
  AA[p].atom[4].hydrogens_atm[1]=19;

  AA[p].atom[5].numHydrogens_atm=2;
  AA[p].atom[5].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[5].numHydrogens_atm);
  AA[p].atom[5].hydrogens_atm[0]=20;
  AA[p].atom[5].hydrogens_atm[1]=21;

  AA[p].atom[6].numHydrogens_atm=2;
  AA[p].atom[6].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[6].numHydrogens_atm);
  AA[p].atom[6].hydrogens_atm[0]=22;
  AA[p].atom[6].hydrogens_atm[1]=23;

  AA[p].atom[7].numHydrogens_atm=1;
  AA[p].atom[7].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[7].numHydrogens_atm);
  AA[p].atom[7].hydrogens_atm[0]=16;

  AA[p].atom[9].numHydrogens_atm=2;
  AA[p].atom[9].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[9].numHydrogens_atm);
  AA[p].atom[9].hydrogens_atm[0]=12;
  AA[p].atom[9].hydrogens_atm[1]=13;

  AA[p].atom[10].numHydrogens_atm=2;
  AA[p].atom[10].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[10].numHydrogens_atm);
  AA[p].atom[10].hydrogens_atm[0]=14;
  AA[p].atom[10].hydrogens_atm[1]=15;

  //FIN DE AMINOACIDO ARGININA
}

/**
* Initializes Serine
*
* @param opt: Rosseta or ICM
*/
void init_SER(Convention opt)
{
  //AMINOACIDO SERINA
  int k, p =SER;

  strcpy( AA[p].aa_name3, "SER" );
  AA[p].aa_name1 = 'S';
  AA[p].mass = 87.08;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = true;
  AA[p].aa_is_nonpolar = false;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = false;

  AA[p].natoms = 11;
  AA[p].nheavyatoms = 6;
  AA[p].nchi = 1;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );

  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;

  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

  AA[p].chi_types = ( int * ) new int[AA[p].nchi];
  AA[p].chi_types[0] = 2;
  AA[p].natoms_EEF1 = 7;

  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " OG " );
  strcpy( AA[p].atom[6].atom_name, " H  " );
  strcpy( AA[p].atom[7].atom_name, " HG " );
  strcpy( AA[p].atom[8].atom_name, " HA " );
  strcpy( AA[p].atom[9].atom_name,  "2HB " );
  strcpy( AA[p].atom[10].atom_name, "3HB " );

  strcpy( AA[p].atom[11].atom_name," OXT" );
  strcpy( AA[p].atom[12].atom_name,"2H  " );
  strcpy( AA[p].atom[13].atom_name,"3H  " );

  switch(opt)
  {
    case Rosseta:
    case Sybil:
     AA[p].atom[0].charge = -0.47; //    N
      AA[p].atom[1].charge = 0.07; //    CA
      AA[p].atom[2].charge = 0.51; //    C
      AA[p].atom[3].charge = -0.51; //    O
      AA[p].atom[4].charge = 0.05; //    CB
      AA[p].atom[5].charge = -0.66; //    OG
      AA[p].atom[6].charge = 0.31; //    H
      AA[p].atom[7].charge = 0.43; //    HG
      AA[p].atom[8].charge = 0.09; //    HA
      AA[p].atom[9].charge = 0.09; //   2HB
      AA[p].atom[10].charge = 0.09; //   3HB

      AA[p].atom[11].charge = -0.67; //   0XT
      AA[p].atom[12].charge = 0.33; //  1H
      AA[p].atom[13].charge = 0.33; //  2H
      break;
    case ICM:
      AA[p].atom[0].charge = -0.356; //    N
      AA[p].atom[1].charge = 0.064; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = 0.13; //    CB
      AA[p].atom[5].charge = -0.31; //    OG
      AA[p].atom[6].charge = 0.176; //    H
      AA[p].atom[7].charge = 0.17; //    HG
      AA[p].atom[8].charge = 0.02; //    HA
      AA[p].atom[9].charge = 0.02; //   2HB
      AA[p].atom[10].charge = 0.02; //   3HB

      AA[p].atom[11].charge = -0.67; //   0XT
      AA[p].atom[12].charge = 0.33; //  1H
      AA[p].atom[13].charge = 0.33; //  2H
     break;
    case EEF1:
      AA[p].atom[0].charge = -0.35; //    N
      AA[p].atom[1].charge =  0.10; //    CA
      AA[p].atom[2].charge =  0.55; //    C
      AA[p].atom[3].charge = -0.55; //    O
      AA[p].atom[4].charge =  0.25; //    CB
      AA[p].atom[5].charge = -0.65; //    OG
      AA[p].atom[6].charge =  0.25; //    H
      AA[p].atom[7].charge =  0.40; //    HG
      AA[p].atom[8].charge =  0.00; //    HA
      AA[p].atom[9].charge =  0.00; //   2HB
      AA[p].atom[10].charge = 0.00; //   3HB
      break;
  }
  /*
   1.0 SER-N:2H    0.33000 HC
   1.0 SER-N:HA    0.10000 HB
   1.0 SER-N:3H    0.33000 HC
   1.0 SER-N:CB    0.05000 CT2
   1.0 SER-N:HB1   0.09000 HA
   1.0 SER-N:1HB   0.09000 HA
   1.0 SER-N:OG   -0.66000 OH1
   1.0 SER-N:HG1   0.43000 H
   1.0 SER-N:C     0.51000 C
   1.0 SER-N:O    -0.51000 O
   1.0 SER-N:N    -0.30000 NH3
   1.0 SER-N:1H    0.33000 HC
   1.0 SER-N:CA    0.21000 CT1

   1.0 SER-C:HA    0.09000 HB
   1.0 SER-C:CB    0.05000 CT2
   1.0 SER-C:HB1   0.09000 HA
   1.0 SER-C:1HB   0.09000 HA
   1.0 SER-C:OG   -0.66000 OH1
   1.0 SER-C:HG1   0.43000 H
   1.0 SER-C:C     0.34000 CC
   1.0 SER-C:O    -0.67000 OC
   1.0 SER-C:OXT  -0.67000 OC
   1.0 SER-C:N    -0.47000 NH1
   1.0 SER-C:HN    0.31000 H
   1.0 SER-C:CA    0.07000 CT1
*/

  switch (opt) {
  	case EEF1:
	  AA[p].atom[0].fullatom_type =  12; // NH1    N
	  AA[p].atom[1].fullatom_type =   5; // CH1E   CA
	  AA[p].atom[2].fullatom_type =   4; // C      C
	  AA[p].atom[3].fullatom_type =  16; // O      O
	  AA[p].atom[4].fullatom_type =   6; // CH2E   CB
	  AA[p].atom[5].fullatom_type =  18; // OH1    OG
	  AA[p].atom[6].fullatom_type =   0; // H      H
	  AA[p].atom[7].fullatom_type =   0; // H      HG
	  AA[p].atom[8].fullatom_type =  -1; // Hapo   HA
	  AA[p].atom[9].fullatom_type =  -1; // Hapo   2HB
	  AA[p].atom[10].fullatom_type = -1; // Hapo   3HB

	  AA[p].atom[11].fullatom_type = 17; // OC     OXT
	  AA[p].atom[12].fullatom_type =  1; // HC     Hpol
	  AA[p].atom[13].fullatom_type =  1; // HC     Hpol
	  break;
  	case Sybil:
  		  AA[p].atom[0].fullatom_type = 40; // Nbb    N
  		  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
  		  AA[p].atom[2].fullatom_type = 10; // CObb   C
  		  AA[p].atom[3].fullatom_type = 56; // OCbb   O
  		  AA[p].atom[4].fullatom_type = 8; // CH2    CB
  		  AA[p].atom[5].fullatom_type = 55; // OH     OG
  		  AA[p].atom[6].fullatom_type = 66; // HNbb   H
  		  AA[p].atom[7].fullatom_type = 66; // Hpol   HG
  		  AA[p].atom[8].fullatom_type = 66; // Hapo   HA
  		  AA[p].atom[9].fullatom_type = 66; // Hapo   2HB
  		  AA[p].atom[10].fullatom_type = 66; // Hapo  3HB

  		  AA[p].atom[11].fullatom_type = 57; // OXT
  		  AA[p].atom[12].fullatom_type = 66; // Hpol
  		  AA[p].atom[13].fullatom_type = 66; // Hpol
  		  break;
  	default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 13; // OH     OG
	  AA[p].atom[6].fullatom_type = 25; // HNbb   H
	  AA[p].atom[7].fullatom_type = 22; // Hpol   HG
	  AA[p].atom[8].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[9].fullatom_type = 23; // Hapo   2HB
	  AA[p].atom[10].fullatom_type = 23; // Hapo  3HB

	  AA[p].atom[11].fullatom_type = 20; // OXT
	  AA[p].atom[12].fullatom_type = 22; // Hpol
	  AA[p].atom[13].fullatom_type = 22; // Hpol
	  break;
  }

  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 6; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 8; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--OG
  AA[p].atom[4].bonded_neighbor[2] = 9; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 10; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 2; // OG
  AA[p].atom[5].bonded_neighbor[0] = 4; // OG--CB
  AA[p].atom[5].bonded_neighbor[1] = 7; // OG--HG
  AA[p].atom[6].nbonded_neighbors = 1; // H
  AA[p].atom[6].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[7].nbonded_neighbors = 1; // HG
  AA[p].atom[7].bonded_neighbor[0] = 5; // HG--OG
  AA[p].atom[8].nbonded_neighbors = 1; // HA
  AA[p].atom[8].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[9].nbonded_neighbors = 1; //2HB
  AA[p].atom[9].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[10].nbonded_neighbors = 1; //3HB
  AA[p].atom[10].bonded_neighbor[0] = 4; //3HB--CB

  AA[p].atom[11].nbonded_neighbors = 1; // OXT
  AA[p].atom[11].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[12].nbonded_neighbors = 1; // 1H
  AA[p].atom[12].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[13].nbonded_neighbors = 1; // 2H
  AA[p].atom[13].bonded_neighbor[0] = 0; // H-N

  for ( k = 0; k < 7; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 7; k < 11; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template atoms used for placing atoms SER

  //bk   template for building HG
  AA[p].atom[7].ta[0] = 5; //   OG
  AA[p].atom[7].ta[1] = 4; //   CB
  AA[p].atom[7].ta[2] = 1; //   CA

  //bk   template for building HA
  AA[p].atom[8].ta[0] = 1; //   CA
  AA[p].atom[8].ta[1] = 0; //   N
  AA[p].atom[8].ta[2] = 2; //   C

  //bk   template for building  2HB
  AA[p].atom[9].ta[0] = 4; //   CB
  AA[p].atom[9].ta[1] = 1; //   CA
  AA[p].atom[9].ta[2] = 5; //   OG

  //bk   template for building 3HB
  AA[p].atom[10].ta[0] = 4; //   CB
  AA[p].atom[10].ta[1] = 1; //   CA
  AA[p].atom[10].ta[2] = 5; //   OG

  //bk   chi angles required to build atoms SER

  //bk   chi angles needed for building  OG
  AA[p].chi_required[0] [5] = true;

  //bk   chi angles needed for building  HG
  AA[p].chi_required[0] [7] = true;
  //        AA[p].chi_required[1][7] =  true;

  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [9] = true;

  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [10] = true;


  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   OG

  //bk   four atoms that define chi angle  2
  //AA[p].chi_atoms[1] [0] = 1; //   CA
  //AA[p].chi_atoms[1] [1] = 4; //   CB
  //AA[p].chi_atoms[1] [2] = 5; //   OG
  //AA[p].chi_atoms[1] [3] = 7; //   HG


  //bk   template coordinates for the amino acid
  //jjh new template from Kevin Karplus

  AA[p].atom[0].icoor[0] = 1.470; //  N
  AA[p].atom[0].icoor[1] = 0.000; //  N
  AA[p].atom[0].icoor[2] = 0.000; //  N
  AA[p].atom[1].icoor[0] = 0.000; //  CA
  AA[p].atom[1].icoor[1] = 0.000; //  CA
  AA[p].atom[1].icoor[2] = 0.000; //  CA
  AA[p].atom[2].icoor[0] = -0.572; //  C
  AA[p].atom[2].icoor[1] = 1.422; //  C
  AA[p].atom[2].icoor[2] = 0.000; //  C
  AA[p].atom[3].icoor[0] = -1.553; //  O
  AA[p].atom[3].icoor[1] = 1.698; //  O
  AA[p].atom[3].icoor[2] = -0.675; //  O
  AA[p].atom[4].icoor[0] = -0.522; //  CB
  AA[p].atom[4].icoor[1] = -0.769; //  CB
  AA[p].atom[4].icoor[2] = 1.198; //  CB
  AA[p].atom[5].icoor[0] = -0.231; //  OG
  AA[p].atom[5].icoor[1] = -2.137; //  OG
  AA[p].atom[5].icoor[2] = 1.113; //  OG
  AA[p].atom[6].icoor[0] = 2.435; //  H
  AA[p].atom[6].icoor[1] = -0.299; //  H
  AA[p].atom[6].icoor[2] = 0.000; //  H
  AA[p].atom[7].icoor[0] = -0.577; //  HG
  AA[p].atom[7].icoor[1] = -2.586; //  HG
  AA[p].atom[7].icoor[2] = 1.888; //  HG
  AA[p].atom[8].icoor[0] = -0.433; //  HA
  AA[p].atom[8].icoor[1] = -0.566; //  HA
  AA[p].atom[8].icoor[2] = -0.826; //  HA
  AA[p].atom[9].icoor[0] = -0.062; //  2HB
  AA[p].atom[9].icoor[1] = -0.364; //  2HB
  AA[p].atom[9].icoor[2] = 2.099; //  2HB
  AA[p].atom[10].icoor[0] = -1.602; //  3HB
  AA[p].atom[10].icoor[1] = -0.639; //  3HB
  AA[p].atom[10].icoor[2] = 1.253; //  3HB

  // atom number for backbone HN
  AA[p].HNpos=6;
  // atom number for backbone HA
  AA[p].HApos=8;

  // number of polar hydrogens
  AA[p].nH_polar_complete=2;
  //atom number for apolar hydrogens
  AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
  AA[p].Hpos_polar_complete[0]=6;
  AA[p].Hpos_polar_complete[1]=7;

  // number of aromatic hydrogens
  AA[p].nH_aromatic_complete=0;

  // number of apolar hydrogens
  AA[p].nH_apolar_complete=3;
  //atom number for apolar hydrogens
  AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
  AA[p].Hpos_apolar_complete[0]=8;
  AA[p].Hpos_apolar_complete[1]=9;
  AA[p].Hpos_apolar_complete[2]=10;

  // number of acceptors
  AA[p].nacceptors=2;
  //acceptor information
  AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
  (AA[p].aBase[1])[0]=5;(AA[p].aBase[1])[1]=7;
  (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;
  (AA[p].accpt_pos[1])[0]=1;(AA[p].accpt_pos[1])[1]=5;

  //Number of Hydrogens connection
  AA[p].nH_hydrogen_connexions=5;
  //atoms hydrogens are connected too
  AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
  AA[p].Hbase[0][0]=6; AA[p].Hbase[0][1]=0;
  AA[p].Hbase[1][0]=7; AA[p].Hbase[1][1]=5;
  AA[p].Hbase[2][0]=8; AA[p].Hbase[2][1]=1;
  AA[p].Hbase[3][0]=9; AA[p].Hbase[3][1]=4;
  AA[p].Hbase[4][0]=10; AA[p].Hbase[4][1]=4;

  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=6;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=8;

  AA[p].atom[4].numHydrogens_atm=2;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=9;
  AA[p].atom[4].hydrogens_atm[1]=10;

  AA[p].atom[5].numHydrogens_atm=1;
  AA[p].atom[5].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[5].numHydrogens_atm);
  AA[p].atom[5].hydrogens_atm[0]=7;

}

//FIN DE AMINOACIDO SERINA

/**
* Initializes Threonine
*
* @param opt: Rosseta or ICM
*/
void init_THR(Convention opt)
{
  //AMINOACIDO TREONINA
  int k, p =THR;

  strcpy( AA[p].aa_name3, "THR" );
  AA[p].aa_name1 = 'T';
  AA[p].mass = 101.107;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = true;
  AA[p].aa_is_nonpolar = false;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = false;

  AA[p].natoms = 14;
  AA[p].nheavyatoms = 7;
  AA[p].nchi = 1;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );
  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;

  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

  AA[p].chi_types = ( int * ) new int[AA[p].nchi];
  AA[p].chi_types[0] = 2;
  AA[p].natoms_EEF1 = 8;


  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " OG1" );
  strcpy( AA[p].atom[6].atom_name, " CG2" );
  strcpy( AA[p].atom[7].atom_name, " H  " );
  strcpy( AA[p].atom[8].atom_name, "1HG " );
  strcpy( AA[p].atom[9].atom_name, " HA " );
  strcpy( AA[p].atom[10].atom_name," HB " );
  strcpy( AA[p].atom[11].atom_name,"1HG2" );
  strcpy( AA[p].atom[12].atom_name,"2HG2" );
  strcpy( AA[p].atom[13].atom_name,"3HG2" );

  strcpy( AA[p].atom[14].atom_name," OXT" );
  strcpy( AA[p].atom[15].atom_name,"2H  " );
  strcpy( AA[p].atom[16].atom_name,"3H  " );

  switch(opt)
  {
    case Rosseta:
    case Sybil:
    AA[p].atom[0].charge = -0.47; //    N
      AA[p].atom[1].charge = 0.07; //    CA
      AA[p].atom[2].charge = 0.51; //    C
      AA[p].atom[3].charge = -0.51; //    O
      AA[p].atom[4].charge = 0.14; //    CB
      AA[p].atom[5].charge = -0.66; //    OG1
      AA[p].atom[6].charge = -0.27; //    CG2
      AA[p].atom[7].charge = 0.31; //    H
      AA[p].atom[8].charge = 0.43; //    HG1
      AA[p].atom[9].charge = 0.09; //    HA
      AA[p].atom[10].charge = 0.09; //    HB
      AA[p].atom[11].charge = 0.09; //    1HG2
      AA[p].atom[12].charge = 0.09; //    2HG2
      AA[p].atom[13].charge = 0.09; //    3HG2

      AA[p].atom[14].charge = -0.67; //   0XT
      AA[p].atom[15].charge = 0.33; //  1H
      AA[p].atom[16].charge = 0.33; //  2H
      break;

    case ICM:
      AA[p].atom[0].charge = -0.356; //    N
      AA[p].atom[1].charge = 0.064; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = 0.16; //    CB
      AA[p].atom[5].charge = -0.31; //    OG1
      AA[p].atom[6].charge = -0.095; //    CG2
      AA[p].atom[7].charge = 0.176; //    H
      AA[p].atom[8].charge = 0.17; //    HG1
      AA[p].atom[9].charge = 0.02; //    HA
      AA[p].atom[10].charge = 0.015; //    HB
      AA[p].atom[11].charge = 0.03; //    1HG2
      AA[p].atom[12].charge = 0.03; //    2HG2
      AA[p].atom[13].charge = 0.03; //    3HG2

      AA[p].atom[14].charge = -0.67; //   0XT
      AA[p].atom[15].charge = 0.33; //  1H
      AA[p].atom[16].charge = 0.33; //  2H
     break;
    case EEF1:
      AA[p].atom[0].charge = -0.35; //    N
      AA[p].atom[1].charge =  0.10; //    CA
      AA[p].atom[2].charge =  0.55; //    C
      AA[p].atom[3].charge = -0.55; //    O
      AA[p].atom[4].charge =  0.25; //    CB
      AA[p].atom[5].charge = -0.65; //    OG1
      AA[p].atom[6].charge =  0.00; //    CG2
      AA[p].atom[7].charge =  0.25; //    H
      AA[p].atom[8].charge =  0.40; //    HG1
      AA[p].atom[9].charge =  0.00; //    HA
      AA[p].atom[10].charge = 0.00; //    HB
      AA[p].atom[11].charge = 0.00; //    1HG2
      AA[p].atom[12].charge = 0.00; //    2HG2
      AA[p].atom[13].charge = 0.00; //    3HG2
     break;
  }
  /*
   1.0 THR-N:2H    0.33000 HC
   1.0 THR-N:HA    0.10000 HB
   1.0 THR-N:3H    0.33000 HC
   1.0 THR-N:CB    0.14000 CT1
   1.0 THR-N:HB    0.09000 HA
   1.0 THR-N:OG1  -0.66000 OH1
   1.0 THR-N:HG1   0.43000 H
   1.0 THR-N:CG2  -0.27000 CT3
   1.0 THR-N:1HG2  0.09000 HA
   1.0 THR-N:2HG2  0.09000 HA
   1.0 THR-N:3HG2  0.09000 HA
   1.0 THR-N:C     0.51000 C
   1.0 THR-N:O    -0.51000 O
   1.0 THR-N:N    -0.30000 NH3
   1.0 THR-N:1H    0.33000 HC
   1.0 THR-N:CA    0.21000 CT1

   1.0 THR-C:HA    0.09000 HB
   1.0 THR-C:CB    0.14000 CT1
   1.0 THR-C:HB    0.09000 HA
   1.0 THR-C:OG1  -0.66000 OH1
   1.0 THR-C:HG1   0.43000 H
   1.0 THR-C:CG2  -0.27000 CT3
   1.0 THR-C:1HG2  0.09000 HA
   1.0 THR-C:2HG2  0.09000 HA
   1.0 THR-C:O    -0.67000 OC
   1.0 THR-C:3HG2  0.09000 HA
   1.0 THR-C:OXT  -0.67000 OC
   1.0 THR-C:C     0.34000 CC
   1.0 THR-C:N    -0.47000 NH1
   1.0 THR-C:HN    0.31000 H
   1.0 THR-C:CA    0.07000 CT1
*/

  switch (opt) {
  	case EEF1:
	  AA[p].atom[0].fullatom_type =  12; // NH1    N
	  AA[p].atom[1].fullatom_type =   5; // CH1E   CA
	  AA[p].atom[2].fullatom_type =   4; // C      C
	  AA[p].atom[3].fullatom_type =  16; // O      O
	  AA[p].atom[4].fullatom_type =   5; // CH1E   CB
	  AA[p].atom[5].fullatom_type =  18; // OH1    OG1
	  AA[p].atom[6].fullatom_type =   7; // CH3E   CG2
	  AA[p].atom[7].fullatom_type =   0; // H      H
	  AA[p].atom[8].fullatom_type =   0; // H      HG1
	  AA[p].atom[9].fullatom_type =  -1; // Hapo   HA
	  AA[p].atom[10].fullatom_type = -1; // Hapo   HB
	  AA[p].atom[11].fullatom_type = -1; // Hapo   1HG2
	  AA[p].atom[12].fullatom_type = -1; // Hapo   2HG2
	  AA[p].atom[13].fullatom_type = -1; // Hapo   3HG2

	  AA[p].atom[14].fullatom_type = 17; // OC     OXT
	  AA[p].atom[15].fullatom_type =  1; // HC     2H
	  AA[p].atom[16].fullatom_type =  1; // HC     3H
	  break;
  	case Sybil:
	  AA[p].atom[0].fullatom_type = 40; // Nbb    N
	  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 10; // CObb   C
	  AA[p].atom[3].fullatom_type = 56; // OCbb   O
	  AA[p].atom[4].fullatom_type = 9; // CH1    CB
	  AA[p].atom[5].fullatom_type = 55; // OH     OG1
	  AA[p].atom[6].fullatom_type = 6; // CH3    CG2
	  AA[p].atom[7].fullatom_type = 66; // HNbb   H
	  AA[p].atom[8].fullatom_type = 66; // Hpol   HG1
	  AA[p].atom[9].fullatom_type = 66; // Hapo   HA
	  AA[p].atom[10].fullatom_type = 66; // Hapo   HB
	  AA[p].atom[11].fullatom_type = 66; // Hapo  1HG2
	  AA[p].atom[12].fullatom_type = 66; // Hapo  2HG2
	  AA[p].atom[13].fullatom_type = 66; // Hapo  3HG2

	  AA[p].atom[14].fullatom_type = 57; // OXT
	  AA[p].atom[15].fullatom_type = 66; // Hpol
	  AA[p].atom[16].fullatom_type = 66; // Hpol
	  break;
	default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 3; // CH1    CB
	  AA[p].atom[5].fullatom_type = 13; // OH     OG1
	  AA[p].atom[6].fullatom_type = 5; // CH3    CG2
	  AA[p].atom[7].fullatom_type = 25; // HNbb   H
	  AA[p].atom[8].fullatom_type = 22; // Hpol   HG1
	  AA[p].atom[9].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[10].fullatom_type = 23; // Hapo   HB
	  AA[p].atom[11].fullatom_type = 23; // Hapo  1HG2
	  AA[p].atom[12].fullatom_type = 23; // Hapo  2HG2
	  AA[p].atom[13].fullatom_type = 23; // Hapo  3HG2

	  AA[p].atom[14].fullatom_type = 20; // OXT
	  AA[p].atom[15].fullatom_type = 22; // Hpol
	  AA[p].atom[16].fullatom_type = 22; // Hpol
	  break;
  }

  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 7; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 9; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--OG1
  AA[p].atom[4].bonded_neighbor[2] = 6; // CB--CG2
  AA[p].atom[4].bonded_neighbor[3] = 10; // CB--HB
  AA[p].atom[5].nbonded_neighbors = 2; // OG1
  AA[p].atom[5].bonded_neighbor[0] = 4; // OG1--CB
  AA[p].atom[5].bonded_neighbor[1] = 8; // OG1--HG1
  AA[p].atom[6].nbonded_neighbors = 4; // CG2
  AA[p].atom[6].bonded_neighbor[0] = 4; // CG2--CB
  AA[p].atom[6].bonded_neighbor[1] = 11; // CG2--1HG2
  AA[p].atom[6].bonded_neighbor[2] = 12; // CG2--2HG2
  AA[p].atom[6].bonded_neighbor[3] = 13; // CG2--3HG2
  AA[p].atom[7].nbonded_neighbors = 1; // H
  AA[p].atom[7].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[8].nbonded_neighbors = 1; // HG1
  AA[p].atom[8].bonded_neighbor[0] = 5; // HG1--OG1
  AA[p].atom[9].nbonded_neighbors = 1; // HA
  AA[p].atom[9].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[10].nbonded_neighbors = 1; // HB
  AA[p].atom[10].bonded_neighbor[0] = 4; // HB--CB
  AA[p].atom[11].nbonded_neighbors = 1; //1HG2
  AA[p].atom[11].bonded_neighbor[0] = 6; //1HG2--CG2
  AA[p].atom[12].nbonded_neighbors = 1; //2HG2
  AA[p].atom[12].bonded_neighbor[0] = 6; //2HG2--CG2
  AA[p].atom[13].nbonded_neighbors = 1; //3HG2
  AA[p].atom[13].bonded_neighbor[0] = 6; //3HG2--CG2

  AA[p].atom[14].nbonded_neighbors = 1; // OXT
  AA[p].atom[14].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[15].nbonded_neighbors = 1; // 1H
  AA[p].atom[15].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[16].nbonded_neighbors = 1; // 2H
  AA[p].atom[16].bonded_neighbor[0] = 0; // H-N


  for ( k = 0; k < 8; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 8; k < 14; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template atoms used for placing atoms THR

  //bk   template for building HG1
  AA[p].atom[8].ta[0] = 5; //   OG1
  AA[p].atom[8].ta[1] = 4; //   CB
  AA[p].atom[8].ta[2] = 1; //   CA

  //bk   template for building HA
  AA[p].atom[9].ta[0] = 1; //   CA
  AA[p].atom[9].ta[1] = 0; //   N
  AA[p].atom[9].ta[2] = 2; //   C

  //bk   template for building  HB
  AA[p].atom[10].ta[0] = 4; //   CB
  AA[p].atom[10].ta[1] = 1; //   CA
  AA[p].atom[10].ta[2] = 5; //   OG1

  //bk   template for building 1HG2
  AA[p].atom[11].ta[0] = 6; //   CG2
  AA[p].atom[11].ta[1] = 4; //   CB
  AA[p].atom[11].ta[2] = 1; //   CA

  //bk   template for building 2HG2
  AA[p].atom[12].ta[0] = 6; //   CG2
  AA[p].atom[12].ta[1] = 4; //   CB
  AA[p].atom[12].ta[2] = 1; //   CA

  //bk   template for building 3HG2
  AA[p].atom[13].ta[0] = 6; //   CG2
  AA[p].atom[13].ta[1] = 4; //   CB
  AA[p].atom[13].ta[2] = 1; //   CA


  //bk   chi angles required to build atoms THR

  //bk   chi angles needed for building  OG1
  AA[p].chi_required[0] [5] = true;

  //bk   chi angles needed for building  CG2
  AA[p].chi_required[0] [6] = true;

  //bk   chi angles needed for building  HG1
  AA[p].chi_required[0] [8] = true;
  //        AA[p].chi_required[1][8] =  true;

  //bk   chi angles needed for building  HB
  AA[p].chi_required[0] [10] = true;

  //bk   chi angles needed for building 1HG2
  AA[p].chi_required[0] [11] = true;

  //bk   chi angles needed for building 2HG2
  AA[p].chi_required[0] [12] = true;

  //bk   chi angles needed for building 3HG2
  AA[p].chi_required[0] [13] = true;



  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   OG1

  //bk   four atoms that define chi angle  2
  //AA[p].chi_atoms[1] [0] = 1; //   CA
  //AA[p].chi_atoms[1] [1] = 4; //   CB
  //AA[p].chi_atoms[1] [2] = 5; //   OG1
  //AA[p].chi_atoms[1] [3] = 8; //   HG1


  //bk   template coordinates for the amino acid
  //jjh new template from Kevin Karplus

  AA[p].atom[0].icoor[0] = 0.0000; //   N
  AA[p].atom[0].icoor[1] = 16.9170; //   N
  AA[p].atom[0].icoor[2] = 58.2740; //   N
  AA[p].atom[1].icoor[0] = 0.0000; //   CA
  AA[p].atom[1].icoor[1] = 17.7300; //   CA
  AA[p].atom[1].icoor[2] = 59.4840; //   CA
  AA[p].atom[2].icoor[0] = 0.0000; //   C
  AA[p].atom[2].icoor[1] = 16.8580; //   C
  AA[p].atom[2].icoor[2] = 60.7340; //   C
  AA[p].atom[3].icoor[0] = 0.0030; //   O
  AA[p].atom[3].icoor[1] = 15.6300; //   O
  AA[p].atom[3].icoor[2] = 60.6450; //   O
  AA[p].atom[4].icoor[0] = -1.2150; //   CB
  AA[p].atom[4].icoor[1] = 18.6750; //   CB
  AA[p].atom[4].icoor[2] = 59.5300; //   CB
  AA[p].atom[5].icoor[0] = -1.1740; //   OG1
  AA[p].atom[5].icoor[1] = 19.4520; //   OG1
  AA[p].atom[5].icoor[2] = 60.7340; //   OG1
  AA[p].atom[6].icoor[0] = -2.5110; //   CG2
  AA[p].atom[6].icoor[1] = 17.8800; //   CG2
  AA[p].atom[6].icoor[2] = 59.4880; //   CG2
  AA[p].atom[7].icoor[0] = -0.0310; //   H
  AA[p].atom[7].icoor[1] = 15.9120; //   H
  AA[p].atom[7].icoor[2] = 58.3720; //   H
  AA[p].atom[8].icoor[0] = -1.9320; //   HG1
  AA[p].atom[8].icoor[1] = 20.0410; //   HG1
  AA[p].atom[8].icoor[2] = 60.7600; //   HG1
  AA[p].atom[9].icoor[0] = 0.9100; //   HA
  AA[p].atom[9].icoor[1] = 18.3290; //   HA
  AA[p].atom[9].icoor[2] = 59.5260; //   HA
  AA[p].atom[10].icoor[0] = -1.1740; //   HB
  AA[p].atom[10].icoor[1] = 19.3470; //   HB
  AA[p].atom[10].icoor[2] = 58.6730; //   HB
  AA[p].atom[11].icoor[0] = -3.3580; //  1HG2
  AA[p].atom[11].icoor[1] = 18.5650; //  1HG2
  AA[p].atom[11].icoor[2] = 59.5210; //  1HG2
  AA[p].atom[12].icoor[0] = -2.5500; //  2HG2
  AA[p].atom[12].icoor[1] = 17.2970; //  2HG2
  AA[p].atom[12].icoor[2] = 58.5680; //  2HG2
  AA[p].atom[13].icoor[0] = -2.5530; //  3HG2
  AA[p].atom[13].icoor[1] = 17.2090; //  3HG2
  AA[p].atom[13].icoor[2] = 60.3450; //  3HG2

  // atom number for backbone HN
  AA[p].HNpos=7;
  // atom number for backbone HA
  AA[p].HApos=9;

  // number of polar hydrogens
  AA[p].nH_polar_complete=2;
  //atom number for apolar hydrogens
  AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
  AA[p].Hpos_polar_complete[0]=7;
  AA[p].Hpos_polar_complete[1]=8;

  // number of aromatic hydrogens
  AA[p].nH_aromatic_complete=0;

  // number of apolar hydrogens
  AA[p].nH_apolar_complete=5;
  //atom number for apolar hydrogens
  AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
  AA[p].Hpos_apolar_complete[0]=9;
  AA[p].Hpos_apolar_complete[1]=10;
  AA[p].Hpos_apolar_complete[2]=11;
  AA[p].Hpos_apolar_complete[3]=12;
  AA[p].Hpos_apolar_complete[4]=13;

  // number of acceptors
  AA[p].nacceptors=2;
  //acceptor information
  AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
  (AA[p].aBase[1])[0]=5;(AA[p].aBase[1])[1]=8;
  (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;
  (AA[p].accpt_pos[1])[0]=1;(AA[p].accpt_pos[1])[1]=5;

  //Number of Hydrogens connection
  AA[p].nH_hydrogen_connexions=7;
  //atoms hydrogens are connected too
  AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
  AA[p].Hbase[0][0]=7; AA[p].Hbase[0][1]=0;
  AA[p].Hbase[1][0]=8; AA[p].Hbase[1][1]=5;
  AA[p].Hbase[2][0]=9; AA[p].Hbase[2][1]=1;
  AA[p].Hbase[3][0]=10; AA[p].Hbase[3][1]=4;
  AA[p].Hbase[4][0]=11; AA[p].Hbase[4][1]=6;
  AA[p].Hbase[5][0]=12; AA[p].Hbase[5][1]=6;
  AA[p].Hbase[6][0]=13; AA[p].Hbase[6][1]=6;

  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=7;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=9;

  AA[p].atom[4].numHydrogens_atm=1;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=10;

  AA[p].atom[5].numHydrogens_atm=1;
  AA[p].atom[5].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[5].numHydrogens_atm);
  AA[p].atom[5].hydrogens_atm[0]=8;

  AA[p].atom[6].numHydrogens_atm=3;
  AA[p].atom[6].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[6].numHydrogens_atm);
  AA[p].atom[6].hydrogens_atm[0]=11;
  AA[p].atom[6].hydrogens_atm[1]=12;
  AA[p].atom[6].hydrogens_atm[2]=13;

  //FIN DE AMINOACIDO TREONINA
}

/**
* Initializes Valine
*
* @param opt: Rosseta or ICM
*/
void init_VAL(Convention opt)
{
  //AMINOCIDO VALINA
  int k, p =VAL;

  strcpy( AA[p].aa_name3, "VAL" );
  AA[p].aa_name1 = 'V';
  AA[p].mass = 99.134;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = false;
  AA[p].aa_is_nonpolar = true;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = false;

  AA[p].natoms = 16;
  AA[p].nheavyatoms = 7;
  AA[p].nchi = 1;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );
  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;

  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

  AA[p].chi_types = ( int * ) new int[AA[p].nchi];
  AA[p].chi_types[0] = 2;
  AA[p].natoms_EEF1 = 7;

  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " CG1" );
  strcpy( AA[p].atom[6].atom_name, " CG2" );
  strcpy( AA[p].atom[7].atom_name, " H  " );
  strcpy( AA[p].atom[8].atom_name, " HA " );
  strcpy( AA[p].atom[9].atom_name, " HB " );
  strcpy( AA[p].atom[10].atom_name, "1HG1" );
  strcpy( AA[p].atom[11].atom_name, "2HG1" );
  strcpy( AA[p].atom[12].atom_name, "3HG1" );
  strcpy( AA[p].atom[13].atom_name, "1HG2" );
  strcpy( AA[p].atom[14].atom_name, "2HG2" );
  strcpy( AA[p].atom[15].atom_name, "3HG2" );

  strcpy( AA[p].atom[16].atom_name," OXT" );
  strcpy( AA[p].atom[17].atom_name,"2H  " );
  strcpy( AA[p].atom[18].atom_name,"3H  " );

  switch(opt)
  {
    case Rosseta:
    case Sybil:
    AA[p].atom[0].charge = -0.47; //   N
      AA[p].atom[1].charge = 0.07; //   CA
      AA[p].atom[2].charge = 0.51; //   C
      AA[p].atom[3].charge = -0.51; //   O
      AA[p].atom[4].charge = -0.09; //   CB
      AA[p].atom[5].charge = -0.27; //   CG1
      AA[p].atom[6].charge = -0.27; //   CG2
      AA[p].atom[7].charge = 0.31; //   H
      AA[p].atom[8].charge = 0.09; //   HA
      AA[p].atom[9].charge = 0.09; //   HB
      AA[p].atom[10].charge = 0.09; //   1HG1
      AA[p].atom[11].charge = 0.09; //   2HG1
      AA[p].atom[12].charge = 0.09; //   3HG1
      AA[p].atom[13].charge = 0.09; //   1HG2
      AA[p].atom[14].charge = 0.09; //   2HG2
      AA[p].atom[15].charge = 0.09; //   3HG2

      AA[p].atom[16].charge = -0.67; //   0XT
      AA[p].atom[17].charge = 0.33; //  1H
      AA[p].atom[18].charge = 0.33; //  2H
    break;
    case ICM:
      AA[p].atom[0].charge = -0.356; //   N
      AA[p].atom[1].charge = 0.064; //   CA
      AA[p].atom[2].charge = 0.45; //   C
      AA[p].atom[3].charge = -0.384; //   O
      AA[p].atom[4].charge = -0.008; //   CB
      AA[p].atom[5].charge = -0.072; //   CG1
      AA[p].atom[6].charge = -0.072; //   CG2
      AA[p].atom[7].charge = 0.176; //   H
      AA[p].atom[8].charge = 0.02; //   HA
      AA[p].atom[9].charge = 0.016; //   HB
      AA[p].atom[10].charge = 0.025; //   1HG1
      AA[p].atom[11].charge = 0.025; //   2HG1
      AA[p].atom[12].charge = 0.025; //   3HG1
      AA[p].atom[13].charge = 0.025; //   1HG2
      AA[p].atom[14].charge = 0.025; //   2HG2
      AA[p].atom[15].charge = 0.025; //   3HG2

      AA[p].atom[16].charge = -0.67; //   0XT
      AA[p].atom[17].charge = 0.33; //  1H
      AA[p].atom[18].charge = 0.33; //  2H
    break;
  case EEF1:
      AA[p].atom[0].charge = -0.35; //   N
      AA[p].atom[1].charge =  0.10; //   CA
      AA[p].atom[2].charge =  0.55; //   C
      AA[p].atom[3].charge = -0.55; //   O
      AA[p].atom[4].charge =  0.00; //   CB
      AA[p].atom[5].charge =  0.00; //   CG1
      AA[p].atom[6].charge =  0.00; //   CG2
      AA[p].atom[7].charge =  0.25; //   H
      AA[p].atom[8].charge =  0.00; //   HA
      AA[p].atom[9].charge =  0.00; //   HB
      AA[p].atom[10].charge = 0.00; //   1HG1
      AA[p].atom[11].charge = 0.00; //   2HG1
      AA[p].atom[12].charge = 0.00; //   3HG1
      AA[p].atom[13].charge = 0.00; //   1HG2
      AA[p].atom[14].charge = 0.00; //   2HG2
      AA[p].atom[15].charge = 0.00; //   3HG2
    break;
  }
  /*
   1.0 VAL-N:2H    0.33000 HC
   1.0 VAL-N:HA    0.10000 HB
   1.0 VAL-N:3H    0.33000 HC
   1.0 VAL-N:CB   -0.09000 CT1
   1.0 VAL-N:HB    0.09000 HA
   1.0 VAL-N:CG1  -0.27000 CT3
   1.0 VAL-N:1HG1  0.09000 HA
   1.0 VAL-N:2HG1  0.09000 HA
   1.0 VAL-N:3HG1  0.09000 HA
   1.0 VAL-N:CG2  -0.27000 CT3
   1.0 VAL-N:1HG2  0.09000 HA
   1.0 VAL-N:2HG2  0.09000 HA
   1.0 VAL-N:3HG2  0.09000 HA
   1.0 VAL-N:N    -0.30000 NH3
   1.0 VAL-N:C     0.51000 C
   1.0 VAL-N:1H    0.33000 HC
   1.0 VAL-N:O    -0.51000 O
   1.0 VAL-N:CA    0.21000 CT1

   1.0 VAL-C:HA    0.09000 HB
   1.0 VAL-C:CB   -0.09000 CT1
   1.0 VAL-C:HB    0.09000 HA
   1.0 VAL-C:CG1  -0.27000 CT3
   1.0 VAL-C:1HG1  0.09000 HA
   1.0 VAL-C:2HG1  0.09000 HA
   1.0 VAL-C:3HG1  0.09000 HA
   1.0 VAL-C:CG2  -0.27000 CT3
   1.0 VAL-C:O    -0.67000 OC
   1.0 VAL-C:1HG2  0.09000 HA
   1.0 VAL-C:OXT  -0.67000 OC
   1.0 VAL-C:2HG2  0.09000 HA
   1.0 VAL-C:3HG2  0.09000 HA
   1.0 VAL-C:N    -0.47000 NH1
   1.0 VAL-C:C     0.34000 CC
   1.0 VAL-C:HN    0.31000 H
   1.0 VAL-C:CA    0.07000 CT1
*/

  switch (opt) {
  	case EEF1:
	  AA[p].atom[0].fullatom_type =  12; // NH1    N
	  AA[p].atom[1].fullatom_type =   5; // CH1E   CA
	  AA[p].atom[2].fullatom_type =   4; // C      C
	  AA[p].atom[3].fullatom_type =  16; // O      O
	  AA[p].atom[4].fullatom_type =   5; // CH1E   CB
	  AA[p].atom[5].fullatom_type =   7; // CH3E   CG1
	  AA[p].atom[6].fullatom_type =   7; // CH3E   CG2
	  AA[p].atom[7].fullatom_type =   0; // H      H
	  AA[p].atom[8].fullatom_type =  -1; // Hapo   HA
	  AA[p].atom[9].fullatom_type =  -1; // Hapo   HB
	  AA[p].atom[10].fullatom_type = -1; // Hapo   1HG1
	  AA[p].atom[11].fullatom_type = -1; // Hapo   2HG1
	  AA[p].atom[12].fullatom_type = -1; // Hapo   3HG1
	  AA[p].atom[13].fullatom_type = -1; // Hapo   1HG2
	  AA[p].atom[14].fullatom_type = -1; // Hapo   2HG2
	  AA[p].atom[15].fullatom_type = -1; // Hapo   3HG2

	  AA[p].atom[16].fullatom_type = 17; // OC     OXT
	  AA[p].atom[17].fullatom_type =  1; // HC     2H
	  AA[p].atom[18].fullatom_type =  1; // HC     3H
	  break;
  	case Sybil:
  		  AA[p].atom[0].fullatom_type = 40; // Nbb    N
  		  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
  		  AA[p].atom[2].fullatom_type = 10; // CObb   C
  		  AA[p].atom[3].fullatom_type = 56; // OCbb   O
  		  AA[p].atom[4].fullatom_type = 9; // CH1    CB
  		  AA[p].atom[5].fullatom_type = 6; // CH3    CG1
  		  AA[p].atom[6].fullatom_type = 6; // CH3    CG2
  		  AA[p].atom[7].fullatom_type = 66; // HNbb   H
  		  AA[p].atom[8].fullatom_type = 66; // Hapo   HA
  		  AA[p].atom[9].fullatom_type = 66; // Hapo   HB
  		  AA[p].atom[10].fullatom_type = 66; // Hapo  1HG1
  		  AA[p].atom[11].fullatom_type = 66; // Hapo  2HG1
  		  AA[p].atom[12].fullatom_type = 66; // Hapo  3HG1
  		  AA[p].atom[13].fullatom_type = 66; // Hapo  1HG2
  		  AA[p].atom[14].fullatom_type = 66; // Hapo  2HG2
  		  AA[p].atom[15].fullatom_type = 66; // Hapo  3HG2

  		  AA[p].atom[16].fullatom_type = 57; // OXT
  		  AA[p].atom[17].fullatom_type = 66; // Hpol
  		  AA[p].atom[18].fullatom_type = 66; // Hpol
  		  break;
  	default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 3; // CH1    CB
	  AA[p].atom[5].fullatom_type = 5; // CH3    CG1
	  AA[p].atom[6].fullatom_type = 5; // CH3    CG2
	  AA[p].atom[7].fullatom_type = 25; // HNbb   H
	  AA[p].atom[8].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[9].fullatom_type = 23; // Hapo   HB
	  AA[p].atom[10].fullatom_type = 23; // Hapo  1HG1
	  AA[p].atom[11].fullatom_type = 23; // Hapo  2HG1
	  AA[p].atom[12].fullatom_type = 23; // Hapo  3HG1
	  AA[p].atom[13].fullatom_type = 23; // Hapo  1HG2
	  AA[p].atom[14].fullatom_type = 23; // Hapo  2HG2
	  AA[p].atom[15].fullatom_type = 23; // Hapo  3HG2

	  AA[p].atom[16].fullatom_type = 20; // OXT
	  AA[p].atom[17].fullatom_type = 22; // Hpol
	  AA[p].atom[18].fullatom_type = 22; // Hpol
	  break;
  }


  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 7; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 8; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--CG1
  AA[p].atom[4].bonded_neighbor[2] = 6; // CB--CG2
  AA[p].atom[4].bonded_neighbor[3] = 9; // CB--HB
  AA[p].atom[5].nbonded_neighbors = 4; // CG1
  AA[p].atom[5].bonded_neighbor[0] = 4; // CG1--CB
  AA[p].atom[5].bonded_neighbor[1] = 10; // CG1--1HG1
  AA[p].atom[5].bonded_neighbor[2] = 11; // CG1--2HG1
  AA[p].atom[5].bonded_neighbor[3] = 12; // CG1--3HG1
  AA[p].atom[6].nbonded_neighbors = 4; // CG2
  AA[p].atom[6].bonded_neighbor[0] = 4; // CG2--CB
  AA[p].atom[6].bonded_neighbor[1] = 13; // CG2--1HG2
  AA[p].atom[6].bonded_neighbor[2] = 14; // CG2--2HG2
  AA[p].atom[6].bonded_neighbor[3] = 15; // CG2--3HG2
  AA[p].atom[7].nbonded_neighbors = 1; // H
  AA[p].atom[7].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[8].nbonded_neighbors = 1; // HA
  AA[p].atom[8].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[9].nbonded_neighbors = 1; // HB
  AA[p].atom[9].bonded_neighbor[0] = 4; // HB--CB
  AA[p].atom[10].nbonded_neighbors = 1; //1HG1
  AA[p].atom[10].bonded_neighbor[0] = 5; //1HG1--CG1
  AA[p].atom[11].nbonded_neighbors = 1; //2HG1
  AA[p].atom[11].bonded_neighbor[0] = 5; //2HG1--CG1
  AA[p].atom[12].nbonded_neighbors = 1; //3HG1
  AA[p].atom[12].bonded_neighbor[0] = 5; //3HG1--CG1
  AA[p].atom[13].nbonded_neighbors = 1; //1HG2
  AA[p].atom[13].bonded_neighbor[0] = 6; //1HG2--CG2
  AA[p].atom[14].nbonded_neighbors = 1; //2HG2
  AA[p].atom[14].bonded_neighbor[0] = 6; //2HG2--CG2
  AA[p].atom[15].nbonded_neighbors = 1; //3HG2
  AA[p].atom[15].bonded_neighbor[0] = 6; //3HG2--CG2

  AA[p].atom[16].nbonded_neighbors = 1; // OXT
  AA[p].atom[16].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[17].nbonded_neighbors = 1; // 1H
  AA[p].atom[17].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[18].nbonded_neighbors = 1; // 2H
  AA[p].atom[18].bonded_neighbor[0] = 0; // H-N


  for ( k = 0; k < 8; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 6; k < 16; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template atoms used for placing atoms VAL

  //bk   template for building HA
  AA[p].atom[8].ta[0] = 1; //   CA
  AA[p].atom[8].ta[1] = 0; //   N
  AA[p].atom[8].ta[2] = 2; //   C

  //bk   template for building HB
  AA[p].atom[9].ta[0] = 4; //   CB
  AA[p].atom[9].ta[1] = 1; //   CA
  AA[p].atom[9].ta[2] = 5; //   CG1

  //bk   template for building  1HG1
  AA[p].atom[10].ta[0] = 5; //   CG1
  AA[p].atom[10].ta[1] = 4; //   CB
  AA[p].atom[10].ta[2] = 1; //   CA

  //bk   template for building 2HG1
  AA[p].atom[11].ta[0] = 5; //   CG1
  AA[p].atom[11].ta[1] = 4; //   CB
  AA[p].atom[11].ta[2] = 1; //   CA

  //bk   template for building 3HG1
  AA[p].atom[12].ta[0] = 5; //   CG1
  AA[p].atom[12].ta[1] = 4; //   CB
  AA[p].atom[12].ta[2] = 1; //   CA

  //bk   template for building  1HG2
  AA[p].atom[13].ta[0] = 6; //   CG2
  AA[p].atom[13].ta[1] = 4; //   CB
  AA[p].atom[13].ta[2] = 1; //   CA

  //bk   template for building 2HG2
  AA[p].atom[14].ta[0] = 6; //   CG2
  AA[p].atom[14].ta[1] = 4; //   CB
  AA[p].atom[14].ta[2] = 1; //   CA

  //bk   template for building 3HG2
  AA[p].atom[15].ta[0] = 6; //   CG2
  AA[p].atom[15].ta[1] = 4; //   CB
  AA[p].atom[15].ta[2] = 1; //   CA

  // VAL
  //bk   chi angles needed for building  CG1
  AA[p].chi_required[0] [5] = true;

  //bk   chi angles needed for building  CG2
  AA[p].chi_required[0] [6] = true;

  //bk   chi angles needed for building  HB
  AA[p].chi_required[0] [9] = true;

  //bk   chi angles needed for building 1HG1
  AA[p].chi_required[0] [10] = true;

  //bk   chi angles needed for building 2HG1
  AA[p].chi_required[0] [11] = true;

  //bk   chi angles needed for building 3HG1
  AA[p].chi_required[0] [12] = true;

  //bk   chi angles needed for building 1HG2
  AA[p].chi_required[0] [13] = true;

  //bk   chi angles needed for building 2HG2
  AA[p].chi_required[0] [14] = true;

  //bk   chi angles needed for building 3HG2
  AA[p].chi_required[0] [15] = true;


  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   CG1

  //bk   template coordinates for the amino acid
  //jjh new template from Kevin Karplus

  AA[p].atom[0].icoor[0] = 0.0000; //   N
  AA[p].atom[0].icoor[1] = 17.5010; //   N
  AA[p].atom[0].icoor[2] = 61.8970; //   N
  AA[p].atom[1].icoor[0] = 0.0000; //   CA
  AA[p].atom[1].icoor[1] = 16.7860; //   CA
  AA[p].atom[1].icoor[2] = 63.1670; //   CA
  AA[p].atom[2].icoor[0] = 0.0000; //   C
  AA[p].atom[2].icoor[1] = 17.7530; //   C
  AA[p].atom[2].icoor[2] = 64.3440; //   C
  AA[p].atom[3].icoor[0] = -0.0020; //   O
  AA[p].atom[3].icoor[1] = 18.9700; //   O
  AA[p].atom[3].icoor[2] = 64.1600; //   O
  AA[p].atom[4].icoor[0] = 1.2150; //   CB
  AA[p].atom[4].icoor[1] = 15.8470; //   CB
  AA[p].atom[4].icoor[2] = 63.2870; //   CB
  AA[p].atom[5].icoor[0] = 1.1920; //   CG1
  AA[p].atom[5].icoor[1] = 15.1140; //   CG1
  AA[p].atom[5].icoor[2] = 64.6200; //   CG1
  AA[p].atom[6].icoor[0] = 1.2370; //   CG2
  AA[p].atom[6].icoor[1] = 14.8540; //   CG2
  AA[p].atom[6].icoor[2] = 62.1350; //   CG2
  AA[p].atom[7].icoor[0] = 0.0310; //   H
  AA[p].atom[7].icoor[1] = 18.5100; //   H
  AA[p].atom[7].icoor[2] = 61.9150; //   H
  AA[p].atom[8].icoor[0] = -0.9090; //   HA
  AA[p].atom[8].icoor[1] = 16.1960; //   HA
  AA[p].atom[8].icoor[2] = 63.2900; //   HA
  AA[p].atom[9].icoor[0] = 2.1280; //   HB
  AA[p].atom[9].icoor[1] = 16.4380; //   HB
  AA[p].atom[9].icoor[2] = 63.2130; //   HB
  AA[p].atom[10].icoor[0] = 2.0580; //  1HG1
  AA[p].atom[10].icoor[1] = 14.4550; //  1HG1
  AA[p].atom[10].icoor[2] = 64.6880; //  1HG1
  AA[p].atom[11].icoor[0] = 1.2230; //  2HG1
  AA[p].atom[11].icoor[1] = 15.8380; //  2HG1
  AA[p].atom[11].icoor[2] = 65.4340; //  2HG1
  AA[p].atom[12].icoor[0] = 0.2800; //  3HG1
  AA[p].atom[12].icoor[1] = 14.5220; //  3HG1
  AA[p].atom[12].icoor[2] = 64.6940; //  3HG1
  AA[p].atom[13].icoor[0] = 2.1020; //  1HG2
  AA[p].atom[13].icoor[1] = 14.1990; //  1HG2
  AA[p].atom[13].icoor[2] = 62.2350; //  1HG2
  AA[p].atom[14].icoor[0] = 0.3250; //  2HG2
  AA[p].atom[14].icoor[1] = 14.2570; //  2HG2
  AA[p].atom[14].icoor[2] = 62.1530; //  2HG2
  AA[p].atom[15].icoor[0] = 1.2990; //  3HG2
  AA[p].atom[15].icoor[1] = 15.3940; //  3HG2
  AA[p].atom[15].icoor[2] = 61.1900; //  3HG2

  // atom number for backbone HN
  AA[p].HNpos=7;
  // atom number for backbone HA
  AA[p].HApos=8;

  // number of polar hydrogens
  AA[p].nH_polar_complete=1;
  //atom number for apolar hydrogens
  AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
  AA[p].Hpos_polar_complete[0]=7;

  // number of aromatic hydrogens
  AA[p].nH_aromatic_complete=0;

  // number of apolar hydrogens
  AA[p].nH_apolar_complete=8;
  //atom number for apolar hydrogens
  AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
  AA[p].Hpos_apolar_complete[0]=8;
  AA[p].Hpos_apolar_complete[1]=9;
  AA[p].Hpos_apolar_complete[2]=10;
  AA[p].Hpos_apolar_complete[3]=11;
  AA[p].Hpos_apolar_complete[4]=12;
  AA[p].Hpos_apolar_complete[5]=13;
  AA[p].Hpos_apolar_complete[6]=14;
  AA[p].Hpos_apolar_complete[7]=15;

  // number of acceptors
  AA[p].nacceptors=1;
  //acceptor information
  AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
  (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;

  //Number of Hydrogens connection
  AA[p].nH_hydrogen_connexions=9;
  //atoms hydrogens are connected too
  AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
  AA[p].Hbase[0][0]=7; AA[p].Hbase[0][1]=0;
  AA[p].Hbase[1][0]=8; AA[p].Hbase[1][1]=1;
  AA[p].Hbase[2][0]=9; AA[p].Hbase[2][1]=4;
  AA[p].Hbase[3][0]=10; AA[p].Hbase[3][1]=5;
  AA[p].Hbase[4][0]=11; AA[p].Hbase[4][1]=5;
  AA[p].Hbase[5][0]=12; AA[p].Hbase[5][1]=5;
  AA[p].Hbase[6][0]=13; AA[p].Hbase[6][1]=6;
  AA[p].Hbase[7][0]=14; AA[p].Hbase[7][1]=6;
  AA[p].Hbase[8][0]=15; AA[p].Hbase[8][1]=6;

  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=7;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=8;

  AA[p].atom[4].numHydrogens_atm=1;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=9;

  AA[p].atom[5].numHydrogens_atm=3;
  AA[p].atom[5].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[5].numHydrogens_atm);
  AA[p].atom[5].hydrogens_atm[0]=10;
  AA[p].atom[5].hydrogens_atm[1]=11;
  AA[p].atom[5].hydrogens_atm[2]=12;

  AA[p].atom[6].numHydrogens_atm=3;
  AA[p].atom[6].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[6].numHydrogens_atm);
  AA[p].atom[6].hydrogens_atm[0]=13;
  AA[p].atom[6].hydrogens_atm[1]=14;
  AA[p].atom[6].hydrogens_atm[2]=15;


  //FIN DE AMINOACIDO VALINA
}

/**
* Initializes Tritopfane
*
* @param opt: Rosseta or ICM
*/
void init_TRP(Convention opt)
{
  //AMINOACIDO TRIPTÓFANO
  int k, p =TRP;

  strcpy( AA[p].aa_name3, "TRP" );
  AA[p].aa_name1 = 'W';

  AA[p].mass = 186.215;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = false;
  AA[p].aa_is_nonpolar = true;
  AA[p].aa_is_aromatic = true;
  AA[p].aa_is_charged = false;

  AA[p].natoms = 24;
  AA[p].nheavyatoms = 14;
  AA[p].nchi = 2;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );
  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;

  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

  AA[p].chi_types = ( int * ) new int[AA[p].nchi];
  AA[p].chi_types[0] = 2;
  AA[p].chi_types[1] = 4;
  AA[p].natoms_EEF1 = 15;

  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " CG " );
  strcpy( AA[p].atom[6].atom_name, " CD1" );
  strcpy( AA[p].atom[7].atom_name, " CD2" );
  strcpy( AA[p].atom[8].atom_name, " NE1" );
  strcpy( AA[p].atom[9].atom_name, " CE2" );
  strcpy( AA[p].atom[10].atom_name, " CE3" );
  strcpy( AA[p].atom[11].atom_name, " CZ2" );
  strcpy( AA[p].atom[12].atom_name, " CZ3" );
  strcpy( AA[p].atom[13].atom_name, " CH2" );
  strcpy( AA[p].atom[14].atom_name, " H  " );
  strcpy( AA[p].atom[15].atom_name, " HE1" );
  strcpy( AA[p].atom[16].atom_name, " HD1" );
  strcpy( AA[p].atom[17].atom_name, " HZ2" );
  strcpy( AA[p].atom[18].atom_name, " HH2" );
  strcpy( AA[p].atom[19].atom_name, " HZ3" );
  strcpy( AA[p].atom[20].atom_name, " HE3" );
  strcpy( AA[p].atom[21].atom_name, " HA " );
  strcpy( AA[p].atom[22].atom_name, "2HB " );
  strcpy( AA[p].atom[23].atom_name, "3HB " );

  strcpy( AA[p].atom[24].atom_name," OXT" );
  strcpy( AA[p].atom[25].atom_name,"2H  " );
  strcpy( AA[p].atom[26].atom_name,"3H  " );

  switch(opt)
  {
    case Rosseta:
    case Sybil:
      AA[p].atom[0].charge = -0.47; //    N
      AA[p].atom[1].charge = 0.07; //    CA
      AA[p].atom[2].charge = 0.51; //    C
      AA[p].atom[3].charge = -0.51; //    O
      AA[p].atom[4].charge = -0.18; //    CB
      AA[p].atom[5].charge = -0.03; //    CG
      AA[p].atom[6].charge = 0.035; //    CD1
      AA[p].atom[7].charge = -0.02; //    CD2
      AA[p].atom[8].charge = -0.61; //    NE1
      AA[p].atom[9].charge = 0.13; //    CE2
      AA[p].atom[10].charge = -0.115; //    CE3
      AA[p].atom[11].charge = -0.115; //    CZ2
      AA[p].atom[12].charge = -0.115; //    CZ3
      AA[p].atom[13].charge = -0.115; //    CH2
      AA[p].atom[14].charge = 0.31; //    H
      AA[p].atom[15].charge = 0.38; //    HE1
      AA[p].atom[16].charge = 0.115; //    HD1
      AA[p].atom[17].charge = 0.115; //    HZ2
      AA[p].atom[18].charge = 0.115; //    HH2
      AA[p].atom[19].charge = 0.115; //    HZ3
      AA[p].atom[20].charge = 0.115; //    HE3
      AA[p].atom[21].charge = 0.09; //    HA
      AA[p].atom[22].charge = 0.09; //    2HB
      AA[p].atom[23].charge = 0.09; //    3HB

      AA[p].atom[24].charge = -0.67; //   0XT
      AA[p].atom[25].charge = 0.33; //  1H
      AA[p].atom[26].charge = 0.33; //  2H
    break;
    case ICM:
      AA[p].atom[0].charge = -0.356; //    N
      AA[p].atom[1].charge = 0.064; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = -0.03; //    CB
      AA[p].atom[5].charge = -0.025; //    CG
      AA[p].atom[6].charge = 0.085; //    CD1
      AA[p].atom[7].charge = 0.005; //    CD2
      AA[p].atom[8].charge = -0.28; //    NE1
      AA[p].atom[9].charge = 0.15; //    CE2
      AA[p].atom[10].charge = -0.08; //    CE3
      AA[p].atom[11].charge = -0.06; //    CZ2
      AA[p].atom[12].charge = -0.02; //    CZ3
      AA[p].atom[13].charge = -0.005; //    CH2
      AA[p].atom[14].charge = 0.176; //    H
      AA[p].atom[15].charge = 0.13; //    HE1
      AA[p].atom[16].charge = 0.000; //    HD1
      AA[p].atom[17].charge = 0.03; //    HZ2
      AA[p].atom[18].charge = 0.01; //    HH2
      AA[p].atom[19].charge = 0.01; //    HZ3
      AA[p].atom[20].charge = 0.08; //    HE3
      AA[p].atom[21].charge = 0.02; //    HA
      AA[p].atom[22].charge = 0.015; //    2HB
      AA[p].atom[23].charge = 0.015; //    3HB

      AA[p].atom[24].charge = -0.67; //   0XT
      AA[p].atom[25].charge = 0.33; //  1H
      AA[p].atom[26].charge = 0.33; //  2H
    break;
   case EEF1:
      AA[p].atom[0].charge = -0.35; //    N
      AA[p].atom[1].charge = 0.10; //    CA
      AA[p].atom[2].charge = 0.55; //    C
      AA[p].atom[3].charge = -0.55; //    O
      AA[p].atom[4].charge = 0.00; //    CB
      AA[p].atom[5].charge = -0.03; //    CG
      AA[p].atom[6].charge = 0.06; //    CD1
      AA[p].atom[7].charge = 0.10; //    CD2
      AA[p].atom[8].charge = -0.36; //    NE1
      AA[p].atom[9].charge = -0.04; //    CE2
      AA[p].atom[10].charge = -0.03; //    CE3
      AA[p].atom[11].charge = 0.00; //    CZ2
      AA[p].atom[12].charge = 0.00; //    CZ3
      AA[p].atom[13].charge = 0.00; //    CH2
      AA[p].atom[14].charge = 0.25; //    H
      AA[p].atom[15].charge = 0.30; //    HE1
      AA[p].atom[16].charge = 0.00; //    HD1
      AA[p].atom[17].charge = 0.00; //    HZ2
      AA[p].atom[18].charge = 0.00; //    HH2
      AA[p].atom[19].charge = 0.00; //    HZ3
      AA[p].atom[20].charge = 0.00; //    HE3
      AA[p].atom[21].charge = 0.00; //    HA
      AA[p].atom[22].charge = 0.00; //    2HB
      AA[p].atom[23].charge = 0.00; //    3HB
    break;
  }

  /*
   1.0 TRP-N:2H    0.33000 HC
   1.0 TRP-N:CZ3  -0.11500 CA
   1.0 TRP-N:HA    0.10000 HB
   1.0 TRP-N:3H    0.33000 HC
   1.0 TRP-N:HZ3   0.11500 HP
   1.0 TRP-N:CB   -0.18000 CT2
   1.0 TRP-N:CZ2  -0.11500 CA
   1.0 TRP-N:HB1   0.09000 HA
   1.0 TRP-N:1HB   0.09000 HA
   1.0 TRP-N:CG   -0.03000 CY
   1.0 TRP-N:CD1   0.03500 CA
   1.0 TRP-N:HD1   0.11500 HP
   1.0 TRP-N:HZ2   0.11500 HP
   1.0 TRP-N:NE1  -0.61000 NY
   1.0 TRP-N:CH2  -0.11500 CA
   1.0 TRP-N:HE1   0.38000 H
   1.0 TRP-N:HH2   0.11500 HP
   1.0 TRP-N:CE2   0.13000 CPT
   1.0 TRP-N:C     0.51000 C
   1.0 TRP-N:CD2  -0.02000 CPT
   1.0 TRP-N:N    -0.30000 NH3
   1.0 TRP-N:O    -0.51000 O
   1.0 TRP-N:CE3  -0.11500 CA
   1.0 TRP-N:1H    0.33000 HC
   1.0 TRP-N:HE3   0.11500 HP
   1.0 TRP-N:CA    0.21000 CT1

   1.0 TRP-C:CZ3  -0.11500 CA
   1.0 TRP-C:HA    0.09000 HB
   1.0 TRP-C:HZ3   0.11500 HP
   1.0 TRP-C:CB   -0.18000 CT2
   1.0 TRP-C:CZ2  -0.11500 CA
   1.0 TRP-C:HB1   0.09000 HA
   1.0 TRP-C:1HB   0.09000 HA
   1.0 TRP-C:CG   -0.03000 CY
   1.0 TRP-C:CD1   0.03500 CA
   1.0 TRP-C:HD1   0.11500 HP
   1.0 TRP-C:HZ2   0.11500 HP
   1.0 TRP-C:NE1  -0.61000 NY
   1.0 TRP-C:O    -0.67000 OC
   1.0 TRP-C:CH2  -0.11500 CA
   1.0 TRP-C:HE1   0.38000 H
   1.0 TRP-C:OXT  -0.67000 OC
   1.0 TRP-C:HH2   0.11500 HP
   1.0 TRP-C:CE2   0.13000 CPT
   1.0 TRP-C:C     0.34000 CC
   1.0 TRP-C:CD2  -0.02000 CPT
   1.0 TRP-C:N    -0.47000 NH1
   1.0 TRP-C:CE3  -0.11500 CA
   1.0 TRP-C:HN    0.31000 H
   1.0 TRP-C:HE3   0.11500 HP
   1.0 TRP-C:CA    0.07000 CT1
*/

  switch (opt) {
  	case EEF1:
	  AA[p].atom[0].fullatom_type = 12; // NH1    N
	  AA[p].atom[1].fullatom_type =  5; // CH1E   CA
	  AA[p].atom[2].fullatom_type =  4; // C      C
	  AA[p].atom[3].fullatom_type = 16; // O      O
	  AA[p].atom[4].fullatom_type = 6; // CH2E    CB
	  AA[p].atom[5].fullatom_type = 24; // CR     CG
	  AA[p].atom[6].fullatom_type = 8; // CR1E    CD1
	  AA[p].atom[7].fullatom_type = 24; // CR     CD2
	  AA[p].atom[8].fullatom_type = 12; // NH1    NE1
	  AA[p].atom[9].fullatom_type = 24; // CR     CE2
	  AA[p].atom[10].fullatom_type = 8; // CR1E   CE3
	  AA[p].atom[11].fullatom_type = 8; // CR1E   CZ2
	  AA[p].atom[12].fullatom_type = 8; // CR1E   CZ3
	  AA[p].atom[13].fullatom_type = 8; // CR1E   CH2
	  AA[p].atom[14].fullatom_type = 0; // H      H
	  AA[p].atom[15].fullatom_type = 0; // H      HE1
	  AA[p].atom[16].fullatom_type = -1; // Haro  HD1
	  AA[p].atom[17].fullatom_type = -1; // Haro  HZ2
	  AA[p].atom[18].fullatom_type = -1; // Haro  HH2
	  AA[p].atom[19].fullatom_type = -1; // Haro  HZ3
	  AA[p].atom[20].fullatom_type = -1; // Haro  HE3
	  AA[p].atom[21].fullatom_type = -1; // Hapo  HA
	  AA[p].atom[22].fullatom_type = -1; // Hapo  2HB
	  AA[p].atom[23].fullatom_type = -1; // Hapo  3HB

	  AA[p].atom[24].fullatom_type = 17; // OC    OXT
	  AA[p].atom[25].fullatom_type =  1; // HC    2H
	  AA[p].atom[26].fullatom_type =  1; // HC    3H
	  break;
  	case Sybil:
	  AA[p].atom[0].fullatom_type = 40; // Nbb    N
	  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 10; // CObb   C
	  AA[p].atom[3].fullatom_type = 56; // OCbb   O
	  AA[p].atom[4].fullatom_type = 8; // CH2    CB
	  AA[p].atom[5].fullatom_type = 11; // aroC   CG
	  AA[p].atom[6].fullatom_type = 11; // aroC   CD1
	  AA[p].atom[7].fullatom_type = 12; // aroC   CD2
	  AA[p].atom[8].fullatom_type = 44; // Ntrp   NE1
	  AA[p].atom[9].fullatom_type = 12; // aroC   CE2
	  AA[p].atom[10].fullatom_type = 12; // aroC   CE3
	  AA[p].atom[11].fullatom_type = 12; // aroC   CZ2
	  AA[p].atom[12].fullatom_type = 12; // aroC   CZ3
	  AA[p].atom[13].fullatom_type = 12; // aroC   CH2
	  AA[p].atom[14].fullatom_type = 66; // HNbb   H
	  AA[p].atom[15].fullatom_type = 66; // Hpol   HE1
	  AA[p].atom[16].fullatom_type = 66; // Haro   HD1
	  AA[p].atom[17].fullatom_type = 66; // Haro   HZ2
	  AA[p].atom[18].fullatom_type = 66; // Haro   HH2
	  AA[p].atom[19].fullatom_type = 66; // Haro   HZ3
	  AA[p].atom[20].fullatom_type = 66; // Haro   HE3
	  AA[p].atom[21].fullatom_type = 66; // Hapo   HA
	  AA[p].atom[22].fullatom_type = 66; // Hapo  2HB
	  AA[p].atom[23].fullatom_type = 66; // Hapo  3HB

	  AA[p].atom[24].fullatom_type = 57; // OXT
	  AA[p].atom[25].fullatom_type = 66; // Hpol
	  AA[p].atom[26].fullatom_type = 66; // Hpol
	  break;
	default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 6; // aroC   CG
	  AA[p].atom[6].fullatom_type = 6; // aroC   CD1
	  AA[p].atom[7].fullatom_type = 6; // aroC   CD2
	  AA[p].atom[8].fullatom_type = 7; // Ntrp   NE1
	  AA[p].atom[9].fullatom_type = 6; // aroC   CE2
	  AA[p].atom[10].fullatom_type = 6; // aroC   CE3
	  AA[p].atom[11].fullatom_type = 6; // aroC   CZ2
	  AA[p].atom[12].fullatom_type = 6; // aroC   CZ3
	  AA[p].atom[13].fullatom_type = 6; // aroC   CH2
	  AA[p].atom[14].fullatom_type = 25; // HNbb   H
	  AA[p].atom[15].fullatom_type = 22; // Hpol   HE1
	  AA[p].atom[16].fullatom_type = 24; // Haro   HD1
	  AA[p].atom[17].fullatom_type = 24; // Haro   HZ2
	  AA[p].atom[18].fullatom_type = 24; // Haro   HH2
	  AA[p].atom[19].fullatom_type = 24; // Haro   HZ3
	  AA[p].atom[20].fullatom_type = 24; // Haro   HE3
	  AA[p].atom[21].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[22].fullatom_type = 23; // Hapo  2HB
	  AA[p].atom[23].fullatom_type = 23; // Hapo  3HB

	  AA[p].atom[24].fullatom_type = 20; // OXT
	  AA[p].atom[25].fullatom_type = 22; // Hpol
	  AA[p].atom[26].fullatom_type = 22; // Hpol
	  break;
  }

  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 14; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 21; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--CG
  AA[p].atom[4].bonded_neighbor[2] = 22; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 23; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 3; // CG
  AA[p].atom[5].bonded_neighbor[0] = 4; // CG--CB
  AA[p].atom[5].bonded_neighbor[1] = 6; // CG--CD1
  AA[p].atom[5].bonded_neighbor[2] = 7; // CG--CD2
  AA[p].atom[6].nbonded_neighbors = 3; // CD1
  AA[p].atom[6].bonded_neighbor[0] = 5; // CD1--CG
  AA[p].atom[6].bonded_neighbor[1] = 8; // CD1--NE1
  AA[p].atom[6].bonded_neighbor[2] = 16; // CD1--HD1
  AA[p].atom[7].nbonded_neighbors = 3; // CD2
  AA[p].atom[7].bonded_neighbor[0] = 5; // CD2--CG
  AA[p].atom[7].bonded_neighbor[1] = 9; // CD2--CE2
  AA[p].atom[7].bonded_neighbor[2] = 10; // CD2--CE3
  AA[p].atom[8].nbonded_neighbors = 3; // NE1
  AA[p].atom[8].bonded_neighbor[0] = 6; // NE1--CD1
  AA[p].atom[8].bonded_neighbor[1] = 9; // NE1--CE2
  AA[p].atom[8].bonded_neighbor[2] = 15; // NE1--HE1
  AA[p].atom[9].nbonded_neighbors = 3; // CE2
  AA[p].atom[9].bonded_neighbor[0] = 7; // CE2--CD2
  AA[p].atom[9].bonded_neighbor[1] = 8; // CE2--NE1
  AA[p].atom[9].bonded_neighbor[2] = 11; // CE2--CZ2
  AA[p].atom[10].nbonded_neighbors = 3; // CE3
  AA[p].atom[10].bonded_neighbor[0] = 7; // CE3--CD2
  AA[p].atom[10].bonded_neighbor[1] = 12; // CE3--CZ3
  AA[p].atom[10].bonded_neighbor[2] = 20; // CE3--HE3
  AA[p].atom[11].nbonded_neighbors = 3; // CZ2
  AA[p].atom[11].bonded_neighbor[0] = 9; // CZ2--CE2
  AA[p].atom[11].bonded_neighbor[1] = 13; // CZ2--CH2
  AA[p].atom[11].bonded_neighbor[2] = 17; // CZ2--HZ2
  AA[p].atom[12].nbonded_neighbors = 3; // CZ3
  AA[p].atom[12].bonded_neighbor[0] = 10; // CZ3--CE3
  AA[p].atom[12].bonded_neighbor[1] = 13; // CZ3--CH2
  AA[p].atom[12].bonded_neighbor[2] = 19; // CZ3--HZ3
  AA[p].atom[13].nbonded_neighbors = 3; // CH2
  AA[p].atom[13].bonded_neighbor[0] = 11; // CH2--CZ2
  AA[p].atom[13].bonded_neighbor[1] = 12; // CH2--CZ3
  AA[p].atom[13].bonded_neighbor[2] = 18; // CH2--HH2
  AA[p].atom[14].nbonded_neighbors = 1; // H
  AA[p].atom[14].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[15].nbonded_neighbors = 1; // HE1
  AA[p].atom[15].bonded_neighbor[0] = 8; // HE1--NE1
  AA[p].atom[16].nbonded_neighbors = 1; // HD1
  AA[p].atom[16].bonded_neighbor[0] = 6; // HD1--CD1
  AA[p].atom[17].nbonded_neighbors = 1; // HZ2
  AA[p].atom[17].bonded_neighbor[0] = 11; // HZ2--CZ2
  AA[p].atom[18].nbonded_neighbors = 1; // HH2
  AA[p].atom[18].bonded_neighbor[0] = 13; // HH2--CH2
  AA[p].atom[19].nbonded_neighbors = 1; // HZ3
  AA[p].atom[19].bonded_neighbor[0] = 12; // HZ3--CZ3
  AA[p].atom[20].nbonded_neighbors = 1; // HE3
  AA[p].atom[20].bonded_neighbor[0] = 10; // HE3--CE3
  AA[p].atom[21].nbonded_neighbors = 1; // HA
  AA[p].atom[21].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[22].nbonded_neighbors = 1; //2HB
  AA[p].atom[22].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[23].nbonded_neighbors = 1; //3HB
  AA[p].atom[23].bonded_neighbor[0] = 4; //3HB--CB

  AA[p].atom[24].nbonded_neighbors = 1; // OXT
  AA[p].atom[24].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[25].nbonded_neighbors = 1; // 1H
  AA[p].atom[25].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[26].nbonded_neighbors = 1; // 2H
  AA[p].atom[26].bonded_neighbor[0] = 0; // H-N


  for ( k = 0; k < 15; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 15; k < 24; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template atoms used for placing atoms TRP

  //bk   template for building HE1
  AA[p].atom[15].ta[0] = 8; //   NE1
  AA[p].atom[15].ta[1] = 6; //   CD1
  AA[p].atom[15].ta[2] = 9; //   CE2

  //bk   template for building HD1
  AA[p].atom[16].ta[0] = 6; //   CD1
  AA[p].atom[16].ta[1] = 5; //   CG
  AA[p].atom[16].ta[2] = 8; //   NE1

  //bk   template for building  HZ2
  AA[p].atom[17].ta[0] = 11; //   CZ2
  AA[p].atom[17].ta[1] = 9; //   CE2
  AA[p].atom[17].ta[2] = 13; //   CH2

  //bk   template for building HH2
  AA[p].atom[18].ta[0] = 13; //   CH2
  AA[p].atom[18].ta[1] = 11; //   CZ2
  AA[p].atom[18].ta[2] = 12; //   CZ3

  //bk   template for building HZ3
  AA[p].atom[19].ta[0] = 12; //   CZ3
  AA[p].atom[19].ta[1] = 10; //   CE3
  AA[p].atom[19].ta[2] = 13; //   CH2

  //bk   template for building HE3
  AA[p].atom[20].ta[0] = 10; //   CE3
  AA[p].atom[20].ta[1] = 7; //   CD2
  AA[p].atom[20].ta[2] = 12; //   CZ3

  //bk   template for building HA
  AA[p].atom[21].ta[0] = 1; //   CA
  AA[p].atom[21].ta[1] = 0; //   N
  AA[p].atom[21].ta[2] = 2; //   C

  //bk   template for building 2HB
  AA[p].atom[22].ta[0] = 4; //   CB
  AA[p].atom[22].ta[1] = 1; //   CA
  AA[p].atom[22].ta[2] = 5; //   CG

  //bk   template for building 3HB
  AA[p].atom[23].ta[0] = 4; //   CB
  AA[p].atom[23].ta[1] = 1; //   CA
  AA[p].atom[23].ta[2] = 5; //   CG

  // TRP
  //bk   chi angles needed for building  CG
  AA[p].chi_required[0] [5] = true;

  //bk   chi angles needed for building  CD1
  AA[p].chi_required[0] [6] = true;
  AA[p].chi_required[1] [6] = true;

  //bk   chi angles needed for building  CD2
  AA[p].chi_required[0] [7] = true;
  AA[p].chi_required[1] [7] = true;

  //bk   chi angles needed for building  NE1
  AA[p].chi_required[0] [8] = true;
  AA[p].chi_required[1] [8] = true;

  //bk   chi angles needed for building  CE2
  AA[p].chi_required[0] [9] = true;
  AA[p].chi_required[1] [9] = true;

  //bk   chi angles needed for building  CE3
  AA[p].chi_required[0] [10] = true;
  AA[p].chi_required[1] [10] = true;

  //bk   chi angles needed for building  CZ2
  AA[p].chi_required[0] [11] = true;
  AA[p].chi_required[1] [11] = true;

  //bk   chi angles needed for building  CZ3
  AA[p].chi_required[0] [12] = true;
  AA[p].chi_required[1] [12] = true;

  //bk   chi angles needed for building  CH2
  AA[p].chi_required[0] [13] = true;
  AA[p].chi_required[1] [13] = true;

  //bk   chi angles needed for building  HE1
  AA[p].chi_required[0] [15] = true;
  AA[p].chi_required[1] [15] = true;

  //bk   chi angles needed for building  HD1
  AA[p].chi_required[0] [16] = true;
  AA[p].chi_required[1] [16] = true;


  //bk   chi angles needed for building  HZ2
  AA[p].chi_required[0] [17] = true;
  AA[p].chi_required[1] [17] = true;


  //bk   chi angles needed for building  HH2
  AA[p].chi_required[0] [18] = true;
  AA[p].chi_required[1] [18] = true;


  //bk   chi angles needed for building  HZ3
  AA[p].chi_required[0] [19] = true;
  AA[p].chi_required[1] [19] = true;


  //bk   chi angles needed for building  HE3
  AA[p].chi_required[0] [20] = true;
  AA[p].chi_required[1] [20] = true;

  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [22] = true;

  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [23] = true;


  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   CG

  //bk   four atoms that define chi angle  2
  AA[p].chi_atoms[1] [0] = 1; //   CA
  AA[p].chi_atoms[1] [1] = 4; //   CB
  AA[p].chi_atoms[1] [2] = 5; //   CG
  AA[p].chi_atoms[1] [3] = 6; //   CD1


  //bk   template coordinates for the amino acid
  //jjh new template from Kevin Karplus

  AA[p].atom[0].icoor[0] = 0.0000; //   N
  AA[p].atom[0].icoor[1] = 17.2040; //   N
  AA[p].atom[0].icoor[2] = 65.5540; //   N
  AA[p].atom[1].icoor[0] = 0.0000; //   CA
  AA[p].atom[1].icoor[1] = 18.0180; //   CA
  AA[p].atom[1].icoor[2] = 66.7640; //   CA
  AA[p].atom[2].icoor[0] = 0.0000; //   C
  AA[p].atom[2].icoor[1] = 17.1480; //   C
  AA[p].atom[2].icoor[2] = 68.0140; //   C
  AA[p].atom[3].icoor[0] = 0.0020; //   O
  AA[p].atom[3].icoor[1] = 15.9210; //   O
  AA[p].atom[3].icoor[2] = 67.9270; //   O
  AA[p].atom[4].icoor[0] = -1.2100; //   CB
  AA[p].atom[4].icoor[1] = 18.9540; //   CB
  AA[p].atom[4].icoor[2] = 66.7770; //   CB
  AA[p].atom[5].icoor[0] = -1.2790; //   CG
  AA[p].atom[5].icoor[1] = 19.8260; //   CG
  AA[p].atom[5].icoor[2] = 67.9940; //   CG
  AA[p].atom[6].icoor[0] = -2.2370; //   CD1
  AA[p].atom[6].icoor[1] = 20.7520; //   CD1
  AA[p].atom[6].icoor[2] = 68.2800; //   CD1
  AA[p].atom[7].icoor[0] = -0.3540; //   CD2
  AA[p].atom[7].icoor[1] = 19.8530; //   CD2
  AA[p].atom[7].icoor[2] = 69.0880; //   CD2
  AA[p].atom[8].icoor[0] = -1.9670; //   NE1
  AA[p].atom[8].icoor[1] = 21.3560; //   NE1
  AA[p].atom[8].icoor[2] = 69.4830; //   NE1
  AA[p].atom[9].icoor[0] = -0.8150; //   CE2
  AA[p].atom[9].icoor[1] = 20.8190; //   CE2
  AA[p].atom[9].icoor[2] = 70.0000; //   CE2
  AA[p].atom[10].icoor[0] = 0.8210; //   CE3
  AA[p].atom[10].icoor[1] = 19.1510; //   CE3
  AA[p].atom[10].icoor[2] = 69.3840; //   CE3
  AA[p].atom[11].icoor[0] = -0.1480; //   CZ2
  AA[p].atom[11].icoor[1] = 21.1040; //   CZ2
  AA[p].atom[11].icoor[2] = 71.1810; //   CZ2
  AA[p].atom[12].icoor[0] = 1.4890; //   CZ3
  AA[p].atom[12].icoor[1] = 19.4360; //   CZ3
  AA[p].atom[12].icoor[2] = 70.5690; //   CZ3
  AA[p].atom[13].icoor[0] = 1.0180; //   CH2
  AA[p].atom[13].icoor[1] = 20.3840; //   CH2
  AA[p].atom[13].icoor[2] = 71.4420; //   CH2
  AA[p].atom[14].icoor[0] = -0.0380; //   H
  AA[p].atom[14].icoor[1] = 16.1990; //   H
  AA[p].atom[14].icoor[2] = 65.6530; //   H
  AA[p].atom[15].icoor[0] = -2.5240; //   HE1
  AA[p].atom[15].icoor[1] = 22.0780; //   HE1
  AA[p].atom[15].icoor[2] = 69.9170; //   HE1
  AA[p].atom[16].icoor[0] = -3.0370; //   HD1
  AA[p].atom[16].icoor[1] = 20.8740; //   HD1
  AA[p].atom[16].icoor[2] = 67.5520; //   HD1
  AA[p].atom[17].icoor[0] = -0.5570; //   HZ2
  AA[p].atom[17].icoor[1] = 21.8650; //   HZ2
  AA[p].atom[17].icoor[2] = 71.8460; //   HZ2
  AA[p].atom[18].icoor[0] = 1.5720; //   HH2
  AA[p].atom[18].icoor[1] = 20.5770; //   HH2
  AA[p].atom[18].icoor[2] = 72.3610; //   HH2
  AA[p].atom[19].icoor[0] = 2.4030; //   HZ3
  AA[p].atom[19].icoor[1] = 18.8830; //   HZ3
  AA[p].atom[19].icoor[2] = 70.7870; //   HZ3
  AA[p].atom[20].icoor[0] = 1.2370; //   HE3
  AA[p].atom[20].icoor[1] = 18.3870; //   HE3
  AA[p].atom[20].icoor[2] = 68.7280; //   HE3
  AA[p].atom[21].icoor[0] = 0.9070; //   HA
  AA[p].atom[21].icoor[1] = 18.6210; //   HA
  AA[p].atom[21].icoor[2] = 66.8040; //   HA
  AA[p].atom[22].icoor[0] = -1.1800; //  2HB
  AA[p].atom[22].icoor[1] = 19.6220; //  2HB
  AA[p].atom[22].icoor[2] = 65.9160; //  2HB
  AA[p].atom[23].icoor[0] = -2.1330; //  3HB
  AA[p].atom[23].icoor[1] = 18.3750; //  3HB
  AA[p].atom[23].icoor[2] = 66.7550; //  3HB

  // atom number for backbone HN
  AA[p].HNpos=14;
  // atom number for backbone HA
  AA[p].HApos=21;

  // number of polar hydrogens
  AA[p].nH_polar_complete=2;
  //atom number for apolar hydrogens
  AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
  AA[p].Hpos_polar_complete[0]=14;
  AA[p].Hpos_polar_complete[1]=15;

  // number of aromatic hydrogens
  AA[p].nH_aromatic_complete=5;
  AA[p].Hpos_aromatic_complete=(int*)malloc(sizeof(int)*AA[p].nH_aromatic_complete);
  AA[p].Hpos_aromatic_complete[0]=16;
  AA[p].Hpos_aromatic_complete[1]=17;
  AA[p].Hpos_aromatic_complete[2]=18;
  AA[p].Hpos_aromatic_complete[3]=19;
  AA[p].Hpos_aromatic_complete[4]=20;

  // number of apolar hydrogens
  AA[p].nH_apolar_complete=3;
  //atom number for apolar hydrogens
  AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
  AA[p].Hpos_apolar_complete[0]=21;
  AA[p].Hpos_apolar_complete[1]=22;
  AA[p].Hpos_apolar_complete[2]=23;

  // number of acceptors
  AA[p].nacceptors=1;
  //acceptor information
  AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
  (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;

  //Number of Hydrogens connection
  AA[p].nH_hydrogen_connexions=10;
  //atoms hydrogens are connected too
  AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
  AA[p].Hbase[0][0]=14; AA[p].Hbase[0][1]=0;
  AA[p].Hbase[1][0]=15; AA[p].Hbase[1][1]=8;
  AA[p].Hbase[2][0]=16; AA[p].Hbase[2][1]=6;
  AA[p].Hbase[3][0]=17; AA[p].Hbase[3][1]=11;
  AA[p].Hbase[4][0]=18; AA[p].Hbase[4][1]=13;
  AA[p].Hbase[5][0]=19; AA[p].Hbase[5][1]=12;
  AA[p].Hbase[6][0]=20; AA[p].Hbase[6][1]=10;
  AA[p].Hbase[7][0]=21; AA[p].Hbase[7][1]=1;
  AA[p].Hbase[8][0]=22; AA[p].Hbase[8][1]=4;
  AA[p].Hbase[9][0]=23; AA[p].Hbase[9][1]=4;


  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=16;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=21;

  AA[p].atom[4].numHydrogens_atm=2;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=22;
  AA[p].atom[4].hydrogens_atm[1]=23;

  AA[p].atom[6].numHydrogens_atm=1;
  AA[p].atom[6].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[6].numHydrogens_atm);
  AA[p].atom[6].hydrogens_atm[0]=16;

  AA[p].atom[8].numHydrogens_atm=1;
  AA[p].atom[8].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[8].numHydrogens_atm);
  AA[p].atom[8].hydrogens_atm[0]=15;

  AA[p].atom[10].numHydrogens_atm=1;
  AA[p].atom[10].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[10].numHydrogens_atm);
  AA[p].atom[10].hydrogens_atm[0]=20;

  AA[p].atom[11].numHydrogens_atm=1;
  AA[p].atom[11].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[11].numHydrogens_atm);
  AA[p].atom[11].hydrogens_atm[0]=17;

  AA[p].atom[12].numHydrogens_atm=1;
  AA[p].atom[12].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[12].numHydrogens_atm);
  AA[p].atom[12].hydrogens_atm[0]=19;

  AA[p].atom[13].numHydrogens_atm=1;
  AA[p].atom[13].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[13].numHydrogens_atm);
  AA[p].atom[13].hydrogens_atm[0]=18;

  //FIN DE AMINOACIDO TRIPTÓFANO

}
/**
* Initializes Tyrosine
*
* @param opt: Rosseta or ICM
*/
void init_TYR(Convention opt)
{
  //AMINOACIDO TIROSINA
  int k, p =TYR;

  strcpy( AA[p].aa_name3, "TYR" );
  AA[p].aa_name1 = 'Y';

  AA[p].mass = 163.178;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = false;
  AA[p].aa_is_nonpolar = true;
  AA[p].aa_is_aromatic = true;
  AA[p].aa_is_charged = false;

  AA[p].natoms = 21;
  AA[p].nheavyatoms = 12;
  AA[p].nchi = 2;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );
  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;

  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

  AA[p].chi_types = ( int * ) new int[AA[p].nchi];
  AA[p].chi_types[0] = 2;
  AA[p].chi_types[1] = 4;
  AA[p].natoms_EEF1 = 13;

  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " CG " );
  strcpy( AA[p].atom[6].atom_name, " CD1" );
  strcpy( AA[p].atom[7].atom_name, " CD2" );
  strcpy( AA[p].atom[8].atom_name, " CE1" );
  strcpy( AA[p].atom[9].atom_name, " CE2" );
  strcpy( AA[p].atom[10].atom_name, " CZ " );
  strcpy( AA[p].atom[11].atom_name, " OH " );
  strcpy( AA[p].atom[12].atom_name, " H  " );
  strcpy( AA[p].atom[13].atom_name, " HH " );
  strcpy( AA[p].atom[14].atom_name, " HD1" );
  strcpy( AA[p].atom[15].atom_name, " HE1" );
  strcpy( AA[p].atom[16].atom_name, " HE2" );
  strcpy( AA[p].atom[17].atom_name, " HD2" );
  strcpy( AA[p].atom[18].atom_name, " HA " );
  strcpy( AA[p].atom[19].atom_name, "2HB " );
  strcpy( AA[p].atom[20].atom_name, "3HB " );

  strcpy( AA[p].atom[21].atom_name," OXT" );
  strcpy( AA[p].atom[22].atom_name,"2H  " );
  strcpy( AA[p].atom[23].atom_name,"3H  " );

  switch(opt)
  {
    case Rosseta:
    case Sybil:
      AA[p].atom[0].charge = -0.47; //    N
      AA[p].atom[1].charge = 0.07; //    CA
      AA[p].atom[2].charge = 0.51; //    C
      AA[p].atom[3].charge = -0.51; //    O
      AA[p].atom[4].charge = -0.18; //    CB
      AA[p].atom[5].charge = 0.00; //    CG
      AA[p].atom[6].charge = -0.115; //    CD1
      AA[p].atom[7].charge = -0.115; //    CD2
      AA[p].atom[8].charge = -0.115; //    CE1
      AA[p].atom[9].charge = -0.115; //    CE2
      AA[p].atom[10].charge = 0.11; //    CZ
      AA[p].atom[11].charge = -0.54; //    OH
      AA[p].atom[12].charge = 0.31; //    H
      AA[p].atom[13].charge = 0.43; //    HH
      AA[p].atom[14].charge = 0.115; //    HD1
      AA[p].atom[15].charge = 0.115; //    HE1
      AA[p].atom[16].charge = 0.115; //    HE2
      AA[p].atom[17].charge = 0.115; //    HD2
      AA[p].atom[18].charge = 0.09; //    HA
      AA[p].atom[19].charge = 0.09; //   2HB
      AA[p].atom[20].charge = 0.09; //   3HB

      AA[p].atom[21].charge = -0.67; //   0XT
      AA[p].atom[22].charge = 0.33; //  1H
      AA[p].atom[23].charge = 0.33; //  2H
    break;
    case ICM:
      AA[p].atom[0].charge = -0.356; //    N
      AA[p].atom[1].charge = 0.064; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = -0.04; //    CB
      AA[p].atom[5].charge = 0.02; //    CG
      AA[p].atom[6].charge = -0.01; //    CD1
      AA[p].atom[7].charge = -0.01; //    CD2
      AA[p].atom[8].charge = -0.06; //    CE1
      AA[p].atom[9].charge = -0.06; //    CE2
      AA[p].atom[10].charge = 0.225; //    CZ
      AA[p].atom[11].charge = -0.33; //    OH
      AA[p].atom[12].charge = 0.176; //    H
      AA[p].atom[13].charge = 0.165; //    HH
      AA[p].atom[14].charge = 0.01; //    HD1
      AA[p].atom[15].charge = 0.03; //    HE1
      AA[p].atom[16].charge = 0.03; //    HE2
      AA[p].atom[17].charge = 0.01; //    HD2
      AA[p].atom[18].charge = 0.02; //    HA
      AA[p].atom[19].charge = 0.025; //   2HB
      AA[p].atom[20].charge = 0.025; //   3HB

      AA[p].atom[21].charge = -0.67; //   0XT
      AA[p].atom[22].charge = 0.33; //  1H
      AA[p].atom[23].charge = 0.33; //  2H

     break;
     case EEF1:
      AA[p].atom[0].charge = -0.35; //    N
      AA[p].atom[1].charge =  0.10; //    CA
      AA[p].atom[2].charge =  0.55; //    C
      AA[p].atom[3].charge = -0.55; //    O
      AA[p].atom[4].charge =  0.00; //    CB
      AA[p].atom[5].charge =  0.00; //    CG
      AA[p].atom[6].charge =  0.00; //    CD1
      AA[p].atom[7].charge =  0.00; //    CD2
      AA[p].atom[8].charge =  0.00; //    CE1
      AA[p].atom[9].charge =  0.00; //    CE2
      AA[p].atom[10].charge = 0.25; //    CZ
      AA[p].atom[11].charge =-0.65; //    OH
      AA[p].atom[12].charge = 0.25; //    H
      AA[p].atom[13].charge = 0.40; //    HH
      AA[p].atom[14].charge = 0.00; //    HD1
      AA[p].atom[15].charge = 0.00; //    HE1
      AA[p].atom[16].charge = 0.00; //    HE2
      AA[p].atom[17].charge = 0.00; //    HD2
      AA[p].atom[18].charge = 0.00; //    HA
      AA[p].atom[19].charge = 0.00; //   2HB
      AA[p].atom[20].charge = 0.00; //   3HB
      break;
  }
  /*
   1.0 TYR-N:2H    0.33000 HC
   1.0 TYR-N:HD2   0.11500 HP
   1.0 TYR-N:HA    0.10000 HB
   1.0 TYR-N:3H    0.33000 HC
   1.0 TYR-N:CE2  -0.11500 CA
   1.0 TYR-N:CB   -0.18000 CT2
   1.0 TYR-N:HE2   0.11500 HP
   1.0 TYR-N:HB1   0.09000 HA
   1.0 TYR-N:1HB   0.09000 HA
   1.0 TYR-N:CG    0.00000 CA
   1.0 TYR-N:CD1  -0.11500 CA
   1.0 TYR-N:HD1   0.11500 HP
   1.0 TYR-N:C     0.51000 C
   1.0 TYR-N:CE1  -0.11500 CA
   1.0 TYR-N:O    -0.51000 O
   1.0 TYR-N:HE1   0.11500 HP
   1.0 TYR-N:CZ    0.11000 CA
   1.0 TYR-N:OH   -0.54000 OH1
   1.0 TYR-N:N    -0.30000 NH3
   1.0 TYR-N:HH    0.43000 H
   1.0 TYR-N:1H    0.33000 HC
   1.0 TYR-N:CD2  -0.11500 CA
   1.0 TYR-N:CA    0.21000 CT1

   1.0 TYR-C:HD2   0.11500 HP
   1.0 TYR-C:HA    0.09000 HB
   1.0 TYR-C:CE2  -0.11500 CA
   1.0 TYR-C:CB   -0.18000 CT2
   1.0 TYR-C:HE2   0.11500 HP
   1.0 TYR-C:HB1   0.09000 HA
   1.0 TYR-C:1HB   0.09000 HA
   1.0 TYR-C:CG    0.00000 CA
   1.0 TYR-C:CD1  -0.11500 CA
   1.0 TYR-C:HD1   0.11500 HP
   1.0 TYR-C:C     0.34000 CC
   1.0 TYR-C:CE1  -0.11500 CA
   1.0 TYR-C:O    -0.67000 OC
   1.0 TYR-C:HE1   0.11500 HP
   1.0 TYR-C:OXT  -0.67000 OC
   1.0 TYR-C:CZ    0.11000 CA
   1.0 TYR-C:OH   -0.54000 OH1
   1.0 TYR-C:N    -0.47000 NH1
   1.0 TYR-C:HH    0.43000 H
   1.0 TYR-C:HN    0.31000 H
   1.0 TYR-C:CD2  -0.11500 CA
   1.0 TYR-C:CA    0.07000 CT1
*/

  switch (opt) {
  	case EEF1:
	  AA[p].atom[0].fullatom_type =  12; // NH1    N
	  AA[p].atom[1].fullatom_type =   5; // CH1E   CA
	  AA[p].atom[2].fullatom_type =   4; // C   C
	  AA[p].atom[3].fullatom_type =  16; // O   O
	  AA[p].atom[4].fullatom_type =   6; // CH2E    CB
	  AA[p].atom[5].fullatom_type =  24; // CR   CG
	  AA[p].atom[6].fullatom_type =   8; // CR1E   CD1
	  AA[p].atom[7].fullatom_type =   8; // CR1E   CD2
	  AA[p].atom[8].fullatom_type =   8; // CR1E   CE1
	  AA[p].atom[9].fullatom_type =   8; // CR1E   CE2
	  AA[p].atom[10].fullatom_type =  4; // C   CZ
	  AA[p].atom[11].fullatom_type = 18; // OH1     OH
	  AA[p].atom[12].fullatom_type =  0; // H   H
	  AA[p].atom[13].fullatom_type =  0; // H   HH
	  AA[p].atom[14].fullatom_type = -1; // Haro   HD1
	  AA[p].atom[15].fullatom_type = -1; // Haro   HE1
	  AA[p].atom[16].fullatom_type = -1; // Haro   HE2
	  AA[p].atom[17].fullatom_type = -1; // Haro   HD2
	  AA[p].atom[18].fullatom_type = -1; // Hapo   HA
	  AA[p].atom[19].fullatom_type = -1; // Hapo  2HB
	  AA[p].atom[20].fullatom_type = -1; // Hapo  3HB

	  AA[p].atom[21].fullatom_type = 17; //OC OXT
	  AA[p].atom[22].fullatom_type =  1; //HC Hpol
	  AA[p].atom[23].fullatom_type =  1; //HC Hpol
	  break;
  	case Sybil:
	  AA[p].atom[0].fullatom_type = 40; // Nbb    N
	  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 10; // CObb   C
	  AA[p].atom[3].fullatom_type = 56; // OCbb   O
	  AA[p].atom[4].fullatom_type = 8; // CH2    CB
	  AA[p].atom[5].fullatom_type = 12; // aroC   CG
	  AA[p].atom[6].fullatom_type = 12; // aroC   CD1
	  AA[p].atom[7].fullatom_type = 12; // aroC   CD2
	  AA[p].atom[8].fullatom_type = 12; // aroC   CE1
	  AA[p].atom[9].fullatom_type = 12; // aroC   CE2
	  AA[p].atom[10].fullatom_type = 12; // aroC   CZ
	  AA[p].atom[11].fullatom_type = 55; // OH     OH
	  AA[p].atom[12].fullatom_type = 66; // HNbb   H
	  AA[p].atom[13].fullatom_type = 66; // Hpol   HH
	  AA[p].atom[14].fullatom_type = 66; // Haro   HD1
	  AA[p].atom[15].fullatom_type = 66; // Haro   HE1
	  AA[p].atom[16].fullatom_type = 66; // Haro   HE2
	  AA[p].atom[17].fullatom_type = 66; // Haro   HD2
	  AA[p].atom[18].fullatom_type = 66; // Hapo   HA
	  AA[p].atom[19].fullatom_type = 66; // Hapo  2HB
	  AA[p].atom[20].fullatom_type = 66; // Hapo  3HB

	  AA[p].atom[21].fullatom_type = 57; // OXT
	  AA[p].atom[22].fullatom_type = 66; // Hpol
	  AA[p].atom[23].fullatom_type = 66; // Hpol
	  break;
  	default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 6; // aroC   CG
	  AA[p].atom[6].fullatom_type = 6; // aroC   CD1
	  AA[p].atom[7].fullatom_type = 6; // aroC   CD2
	  AA[p].atom[8].fullatom_type = 6; // aroC   CE1
	  AA[p].atom[9].fullatom_type = 6; // aroC   CE2
	  AA[p].atom[10].fullatom_type = 6; // aroC   CZ
	  AA[p].atom[11].fullatom_type = 13; // OH     OH
	  AA[p].atom[12].fullatom_type = 25; // HNbb   H
	  AA[p].atom[13].fullatom_type = 22; // Hpol   HH
	  AA[p].atom[14].fullatom_type = 24; // Haro   HD1
	  AA[p].atom[15].fullatom_type = 24; // Haro   HE1
	  AA[p].atom[16].fullatom_type = 24; // Haro   HE2
	  AA[p].atom[17].fullatom_type = 24; // Haro   HD2
	  AA[p].atom[18].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[19].fullatom_type = 23; // Hapo  2HB
	  AA[p].atom[20].fullatom_type = 23; // Hapo  3HB

	  AA[p].atom[21].fullatom_type = 20; // OXT
	  AA[p].atom[22].fullatom_type = 22; // Hpol
	  AA[p].atom[23].fullatom_type = 22; // Hpol
	  break;
  }

  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 12; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 18; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--CG
  AA[p].atom[4].bonded_neighbor[2] = 19; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 20; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 3; // CG
  AA[p].atom[5].bonded_neighbor[0] = 4; // CG--CB
  AA[p].atom[5].bonded_neighbor[1] = 6; // CG--CD1
  AA[p].atom[5].bonded_neighbor[2] = 7; // CG--CD2
  AA[p].atom[6].nbonded_neighbors = 3; // CD1
  AA[p].atom[6].bonded_neighbor[0] = 5; // CD1--CG
  AA[p].atom[6].bonded_neighbor[1] = 8; // CD1--CE1
  AA[p].atom[6].bonded_neighbor[2] = 14; // CD1--HD1
  AA[p].atom[7].nbonded_neighbors = 3; // CD2
  AA[p].atom[7].bonded_neighbor[0] = 5; // CD2--CG
  AA[p].atom[7].bonded_neighbor[1] = 9; // CD2--CE2
  AA[p].atom[7].bonded_neighbor[2] = 17; // CD2--HD2
  AA[p].atom[8].nbonded_neighbors = 3; // CE1
  AA[p].atom[8].bonded_neighbor[0] = 6; // CE1--CD1
  AA[p].atom[8].bonded_neighbor[1] = 10; // CE1--CZ
  AA[p].atom[8].bonded_neighbor[2] = 15; // CE1--HE1
  AA[p].atom[9].nbonded_neighbors = 3; // CE2
  AA[p].atom[9].bonded_neighbor[0] = 7; // CE2--CD2
  AA[p].atom[9].bonded_neighbor[1] = 10; // CE2--CZ
  AA[p].atom[9].bonded_neighbor[2] = 16; // CE2--HE2
  AA[p].atom[10].nbonded_neighbors = 3; // CZ
  AA[p].atom[10].bonded_neighbor[0] = 8; // CZ--CE1
  AA[p].atom[10].bonded_neighbor[1] = 9; // CZ--CE2
  AA[p].atom[10].bonded_neighbor[2] = 11; // CZ--OH
  AA[p].atom[11].nbonded_neighbors = 2; // OH
  AA[p].atom[11].bonded_neighbor[0] = 10; // OH--CZ
  AA[p].atom[11].bonded_neighbor[1] = 13; // OH--HH
  AA[p].atom[12].nbonded_neighbors = 1; // H
  AA[p].atom[12].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[13].nbonded_neighbors = 1; // HH
  AA[p].atom[13].bonded_neighbor[0] = 11; // HH--OH
  AA[p].atom[14].nbonded_neighbors = 1; // HD1
  AA[p].atom[14].bonded_neighbor[0] = 6; // HD1--CD1
  AA[p].atom[15].nbonded_neighbors = 1; // HE1
  AA[p].atom[15].bonded_neighbor[0] = 8; // HE1--CE1
  AA[p].atom[16].nbonded_neighbors = 1; // HE2
  AA[p].atom[16].bonded_neighbor[0] = 9; // HE2--CE2
  AA[p].atom[17].nbonded_neighbors = 1; // HD2
  AA[p].atom[17].bonded_neighbor[0] = 7; // HD2--CD2
  AA[p].atom[18].nbonded_neighbors = 1; // HA
  AA[p].atom[18].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[19].nbonded_neighbors = 1; //2HB
  AA[p].atom[19].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[20].nbonded_neighbors = 1; //3HB
  AA[p].atom[20].bonded_neighbor[0] = 4; //3HB--CB

  AA[p].atom[21].nbonded_neighbors = 1; // OXT
  AA[p].atom[21].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[22].nbonded_neighbors = 1; // 1H
  AA[p].atom[22].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[23].nbonded_neighbors = 1; // 2H
  AA[p].atom[23].bonded_neighbor[0] = 0; // H-N


  for ( k = 0; k < 13; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 6; k < 21; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template atoms used for placing atoms TYR

  //bk   template for building HH
  AA[p].atom[13].ta[0] = 11; //   OH
  AA[p].atom[13].ta[1] = 10; //   CZ
  AA[p].atom[13].ta[2] = 8; //   CE1

  //bk   template for building HD1
  AA[p].atom[14].ta[0] = 6; //   CD1
  AA[p].atom[14].ta[1] = 5; //   CG
  AA[p].atom[14].ta[2] = 8; //   CE1

  //bk   template for building  HE1
  AA[p].atom[15].ta[0] = 8; //   CE1
  AA[p].atom[15].ta[1] = 6; //   CD1
  AA[p].atom[15].ta[2] = 10; //   CZ

  //bk   template for building HE2
  AA[p].atom[16].ta[0] = 9; //   CE2
  AA[p].atom[16].ta[1] = 7; //   CD2
  AA[p].atom[16].ta[2] = 10; //   CZ

  //bk   template for building HD2
  AA[p].atom[17].ta[0] = 7; //   CD2
  AA[p].atom[17].ta[1] = 5; //   CG
  AA[p].atom[17].ta[2] = 9; //   CE2

  //bk   template for building HA
  AA[p].atom[18].ta[0] = 1; //   CA
  AA[p].atom[18].ta[1] = 0; //   N
  AA[p].atom[18].ta[2] = 2; //   C

  //bk   template for building 2HB
  AA[p].atom[19].ta[0] = 4; //   CB
  AA[p].atom[19].ta[1] = 1; //   CA
  AA[p].atom[19].ta[2] = 5; //   CG

  //bk   template for building 3HB
  AA[p].atom[20].ta[0] = 4; //   CB
  AA[p].atom[20].ta[1] = 1; //   CA
  AA[p].atom[20].ta[2] = 5; //   CG


  //bk   chi angles needed for building  CG
  AA[p].chi_required[0] [5] = true;

  //bk   chi angles needed for building  CD1
  AA[p].chi_required[0] [6] = true;
  AA[p].chi_required[1] [6] = true;

  //bk   chi angles needed for building  CD2
  AA[p].chi_required[0] [7] = true;
  AA[p].chi_required[1] [7] = true;

  //bk   chi angles needed for building  CE1
  AA[p].chi_required[0] [8] = true;
  AA[p].chi_required[1] [8] = true;

  //bk   chi angles needed for building  CE2
  AA[p].chi_required[0] [9] = true;
  AA[p].chi_required[1] [9] = true;

  //bk   chi angles needed for building  CZ
  AA[p].chi_required[0] [10] = true;
  AA[p].chi_required[1] [10] = true;

  //bk   chi angles needed for building  OH
  AA[p].chi_required[0] [11] = true;
  AA[p].chi_required[1] [11] = true;


  //bk   chi angles needed for building  HH
  AA[p].chi_required[0] [13] = true;
  AA[p].chi_required[1] [13] = true;
  //        AA[p].chi_required[2][13] =  true;


  //bk   chi angles needed for building  HD1
  AA[p].chi_required[0] [14] = true;
  AA[p].chi_required[1] [14] = true;

  //bk   chi angles needed for building  HE1
  AA[p].chi_required[0] [15] = true;
  AA[p].chi_required[1] [15] = true;

  //bk   chi angles needed for building  HE2
  AA[p].chi_required[0] [16] = true;
  AA[p].chi_required[1] [16] = true;

  //bk   chi angles needed for building  HD2
  AA[p].chi_required[0] [17] = true;
  AA[p].chi_required[1] [17] = true;

  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [19] = true;

  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [20] = true;


  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   CG

  //bk   four atoms that define chi angle  2
  AA[p].chi_atoms[1] [0] = 1; //   CA
  AA[p].chi_atoms[1] [1] = 4; //   CB
  AA[p].chi_atoms[1] [2] = 5; //   CG
  AA[p].chi_atoms[1] [3] = 6; //   CD1

  //bk   four atoms that define chi angle  3
  //AA[p].chi_atoms[2] [0] = 9; //   CE2
  //AA[p].chi_atoms[2] [1] = 10; //   CZ
  //AA[p].chi_atoms[2] [2] = 11; //   OH
  //AA[p].chi_atoms[2] [3] = 13; //   HH


  //bk   template coordinates for the amino acid
  //jjh new template from Kevin Karplus

  AA[p].atom[0].icoor[0] = 0.0000; //   N
  AA[p].atom[0].icoor[1] = 17.7910; //   N
  AA[p].atom[0].icoor[2] = 69.1770; //   N
  AA[p].atom[1].icoor[0] = 0.0000; //   CA
  AA[p].atom[1].icoor[1] = 17.0770; //   CA
  AA[p].atom[1].icoor[2] = 70.4470; //   CA
  AA[p].atom[2].icoor[0] = 0.0000; //   C
  AA[p].atom[2].icoor[1] = 18.0440; //   C
  AA[p].atom[2].icoor[2] = 71.6230; //   C
  AA[p].atom[3].icoor[0] = -0.0010; //   O
  AA[p].atom[3].icoor[1] = 19.2610; //   O
  AA[p].atom[3].icoor[2] = 71.4380; //   O
  AA[p].atom[4].icoor[0] = 1.2090; //   CB
  AA[p].atom[4].icoor[1] = 16.1430; //   CB
  AA[p].atom[4].icoor[2] = 70.5360; //   CB
  AA[p].atom[5].icoor[0] = 1.2850; //   CG
  AA[p].atom[5].icoor[1] = 15.3600; //   CG
  AA[p].atom[5].icoor[2] = 71.8280; //   CG
  AA[p].atom[6].icoor[0] = 2.3230; //   CD1
  AA[p].atom[6].icoor[1] = 14.4700; //   CD1
  AA[p].atom[6].icoor[2] = 72.0620; //   CD1
  AA[p].atom[7].icoor[0] = 0.3190; //   CD2
  AA[p].atom[7].icoor[1] = 15.5150; //   CD2
  AA[p].atom[7].icoor[2] = 72.8110; //   CD2
  AA[p].atom[8].icoor[0] = 2.3970; //   CE1
  AA[p].atom[8].icoor[1] = 13.7520; //   CE1
  AA[p].atom[8].icoor[2] = 73.2400; //   CE1
  AA[p].atom[9].icoor[0] = 0.3830; //   CE2
  AA[p].atom[9].icoor[1] = 14.8030; //   CE2
  AA[p].atom[9].icoor[2] = 73.9930; //   CE2
  AA[p].atom[10].icoor[0] = 1.4240; //   CZ
  AA[p].atom[10].icoor[1] = 13.9220; //   CZ
  AA[p].atom[10].icoor[2] = 74.2040; //   CZ
  AA[p].atom[11].icoor[0] = 1.4930; //   OH
  AA[p].atom[11].icoor[1] = 13.2110; //   OH
  AA[p].atom[11].icoor[2] = 75.3800; //   OH
  AA[p].atom[12].icoor[0] = 0.0390; //   H
  AA[p].atom[12].icoor[1] = 18.8010; //   H
  AA[p].atom[12].icoor[2] = 69.1950; //   H
  AA[p].atom[13].icoor[0] = 0.7680; //   HH
  AA[p].atom[13].icoor[1] = 13.4050; //   HH
  AA[p].atom[13].icoor[2] = 75.9790; //   HH
  AA[p].atom[14].icoor[0] = 3.0880; //   HD1
  AA[p].atom[14].icoor[1] = 14.3410; //   HD1
  AA[p].atom[14].icoor[2] = 71.2960; //   HD1
  AA[p].atom[15].icoor[0] = 3.2180; //   HE1
  AA[p].atom[15].icoor[1] = 13.0560; //   HE1
  AA[p].atom[15].icoor[2] = 73.4100; //   HE1
  AA[p].atom[16].icoor[0] = -0.3870; //   HE2
  AA[p].atom[16].icoor[1] = 14.9380; //   HE2
  AA[p].atom[16].icoor[2] = 74.7530; //   HE2
  AA[p].atom[17].icoor[0] = -0.5010; //   HD2
  AA[p].atom[17].icoor[1] = 16.2120; //   HD2
  AA[p].atom[17].icoor[2] = 72.6380; //   HD2
  AA[p].atom[18].icoor[0] = -0.9080; //   HA
  AA[p].atom[18].icoor[1] = 16.4800; //   HA
  AA[p].atom[18].icoor[2] = 70.5360; //   HA
  AA[p].atom[19].icoor[0] = 1.1450; //  2HB
  AA[p].atom[19].icoor[1] = 15.4510; //  2HB
  AA[p].atom[19].icoor[2] = 69.6950; //  2HB
  AA[p].atom[20].icoor[0] = 2.1010; //  3HB
  AA[p].atom[20].icoor[1] = 16.7600; //  3HB
  AA[p].atom[20].icoor[2] = 70.4320; //  3HB

  // atom number for backbone HN
  AA[p].HNpos=12;
  // atom number for backbone HA
  AA[p].HApos=18;

  // number of polar hydrogens
  AA[p].nH_polar_complete=2;
  //atom number for apolar hydrogens
  AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
  AA[p].Hpos_polar_complete[0]=12;
  AA[p].Hpos_polar_complete[1]=13;

  // number of aromatic hydrogens
  AA[p].nH_aromatic_complete=4;
  AA[p].Hpos_aromatic_complete=(int*)malloc(sizeof(int)*AA[p].nH_aromatic_complete);
  AA[p].Hpos_aromatic_complete[0]=14;
  AA[p].Hpos_aromatic_complete[1]=15;
  AA[p].Hpos_aromatic_complete[2]=16;
  AA[p].Hpos_aromatic_complete[3]=17;

  // number of apolar hydrogens
  AA[p].nH_apolar_complete=3;
  //atom number for apolar hydrogens
  AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
  AA[p].Hpos_apolar_complete[0]=18;
  AA[p].Hpos_apolar_complete[1]=19;
  AA[p].Hpos_apolar_complete[2]=20;

  // number of acceptors
  AA[p].nacceptors=2;
  //acceptor information
  AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
  (AA[p].aBase[1])[0]=11;(AA[p].aBase[1])[1]=13;
  (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;
  (AA[p].accpt_pos[1])[0]=1;(AA[p].accpt_pos[1])[1]=11;

  //Number of Hydrogens connection
  AA[p].nH_hydrogen_connexions=9;
  //atoms hydrogens are connected too
  AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
  AA[p].Hbase[0][0]=12; AA[p].Hbase[0][1]=0;
  AA[p].Hbase[1][0]=13; AA[p].Hbase[1][1]=11;
  AA[p].Hbase[2][0]=14; AA[p].Hbase[2][1]=6;
  AA[p].Hbase[3][0]=15; AA[p].Hbase[3][1]=8;
  AA[p].Hbase[4][0]=16; AA[p].Hbase[4][1]=9;
  AA[p].Hbase[5][0]=17; AA[p].Hbase[5][1]=7;
  AA[p].Hbase[6][0]=18; AA[p].Hbase[6][1]=1;
  AA[p].Hbase[7][0]=19; AA[p].Hbase[7][1]=4;
  AA[p].Hbase[8][0]=20; AA[p].Hbase[8][1]=4;

  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=12;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=18;

  AA[p].atom[4].numHydrogens_atm=2;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=19;
  AA[p].atom[4].hydrogens_atm[1]=20;

  AA[p].atom[6].numHydrogens_atm=1;
  AA[p].atom[6].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[6].numHydrogens_atm);
  AA[p].atom[6].hydrogens_atm[0]=14;

  AA[p].atom[7].numHydrogens_atm=1;
  AA[p].atom[7].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[7].numHydrogens_atm);
  AA[p].atom[7].hydrogens_atm[0]=17;

  AA[p].atom[8].numHydrogens_atm=1;
  AA[p].atom[8].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[8].numHydrogens_atm);
  AA[p].atom[8].hydrogens_atm[0]=15;

  AA[p].atom[9].numHydrogens_atm=1;
  AA[p].atom[9].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[9].numHydrogens_atm);
  AA[p].atom[9].hydrogens_atm[0]=16;

  AA[p].atom[11].numHydrogens_atm=1;
  AA[p].atom[11].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[11].numHydrogens_atm);
  AA[p].atom[11].hydrogens_atm[0]=13;

  //FIN DE AMINOACIDO TIROSINA
}

//FIN DE AMINOACIDO CYSTEINA
/**
* Initializes Aspartate
*
* @param opt: Rosseta or ICM
*/
void init_ASH(Convention opt)
{
  //AMINOACIDO NEUTRAL ASPARTATO
  int k, p = ASH;

  strcpy( AA[p].aa_name3, "ASH" );
  AA[p].aa_name1 = 'D';
  AA[p].mass = 114.083;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = true;
  AA[p].aa_is_nonpolar = false;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = true;

  AA[p].natoms = 13; //NEW ADDED
  AA[p].nheavyatoms = 8;
  AA[p].nchi = 2;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );

  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;
  /* for ( int chino = 0; chino < AA[p].nchi; chino++ ) for ( int i = 0 ; i < AA[p].natoms; i++)
  if (AA[p].chi_required[chino][i]) printf( "--->%d %d true\n", chino, i); else  printf("--->%d %d false\n", chino, i); */

  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms +3) * sizeof( t_aa_atom_p ) );



  strcpy( AA[p].atom[0].atom_name, " N  " );AA[p].chi_types = ( int * ) new int[AA[p].nchi];
  AA[p].chi_types[0] = 2;
  AA[p].chi_types[1] = 4;

  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " CG " );
  strcpy( AA[p].atom[6].atom_name, " OD1" );
  strcpy( AA[p].atom[7].atom_name, " OD2" );
  strcpy( AA[p].atom[8].atom_name, " H  " );
  strcpy( AA[p].atom[9].atom_name, " HA " );
  strcpy( AA[p].atom[10].atom_name, "2HB " );
  strcpy( AA[p].atom[11].atom_name, "3HB " );
  strcpy( AA[p].atom[12].atom_name, "2HD " ); //NEW ATOM

  strcpy( AA[p].atom[13].atom_name," OXT" );
  strcpy( AA[p].atom[14].atom_name,"2H  " );
  strcpy( AA[p].atom[15].atom_name,"3H  " );


  switch(opt)
  {
    //    like CHARMM ASPP
    case Rosseta:
    case Sybil:
      AA[p].atom[0].charge = -0.47; //    N  <--- Nt(-0.30)
      AA[p].atom[1].charge = 0.07; //    CA  <--- Nt(0.21)
      AA[p].atom[2].charge = 0.51; //    C  <--- Ct(0.34)
      AA[p].atom[3].charge = -0.51; //    O  <--- Ct(-0.67)
      AA[p].atom[4].charge = -0.21; //    CB
      AA[p].atom[5].charge = 0.75; //    CG
      AA[p].atom[6].charge = -0.55; //    OD1
      AA[p].atom[7].charge = -0.61; //    OD2
      AA[p].atom[8].charge = 0.31; //    H  <--- Nt(0.33)
      AA[p].atom[9].charge = 0.09; //    HA  <--- Nt(0.10)
      AA[p].atom[10].charge = 0.09; //    2HB
      AA[p].atom[11].charge = 0.09; //    3HB
      AA[p].atom[12].charge = 0.44; //    2HD NEW ADDED - like CHARMM ASPP

      AA[p].atom[13].charge = -0.67; //   0XT
      AA[p].atom[14].charge = 0.33; //  1H
      AA[p].atom[15].charge = 0.33; //  2H
      break;
		//A REVISAR
    case ICM:
      AA[p].atom[0].charge = -0.356; //    N
      AA[p].atom[1].charge = 0.06; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = -0.17; //    CB
      AA[p].atom[5].charge = 0.5; //    CG
      AA[p].atom[6].charge = -0.57; //    OD1
      AA[p].atom[7].charge = -0.57; //    OD2
      AA[p].atom[8].charge = 0.176; //    H
      AA[p].atom[9].charge = 0.024; //    HA
      AA[p].atom[10].charge = 0.02; //    2HB
      AA[p].atom[11].charge = 0.02; //    3HB
      AA[p].atom[12].charge = 0.02; //    2HD NEW ADDED BUT UNKNOWN

      AA[p].atom[13].charge = -0.67; //   0XT
      AA[p].atom[14].charge = 0.33; //  1H
      AA[p].atom[15].charge = 0.33; //  2H

    break;
    default:
      break;

  }

  /*   Nt/Ct - Charge Table

   1.0 ASP-N:2H    0.33000 HC
   1.0 ASP-N:HA    0.10000 HB
   1.0 ASP-N:3H    0.33000 HC
   1.0 ASP-N:CB   -0.28000 CT2
   1.0 ASP-N:HB1   0.09000 HA
   1.0 ASP-N:1HB   0.09000 HA
   1.0 ASP-N:CG    0.62000 CC
   1.0 ASP-N:OD1  -0.76000 OC
   1.0 ASP-N:OD2  -0.76000 OC
   1.0 ASP-N:C     0.51000 C
   1.0 ASP-N:O    -0.51000 O
   1.0 ASP-N:N    -0.30000 NH3
   1.0 ASP-N:1H    0.33000 HC
   1.0 ASP-N:CA    0.21000 CT1

   1.0 ASP-C:HA    0.09000 HB
   1.0 ASP-C:CB   -0.28000 CT2
   1.0 ASP-C:HB1   0.09000 HA
   1.0 ASP-C:1HB   0.09000 HA
   1.0 ASP-C:CG    0.62000 CC
   1.0 ASP-C:OD1  -0.76000 OC
   1.0 ASP-C:OD2  -0.76000 OC
   1.0 ASP-C:C     0.34000 CC
   1.0 ASP-C:O    -0.67000 OC
   1.0 ASP-C:OXT  -0.67000 OC
   1.0 ASP-C:N    -0.47000 NH1
   1.0 ASP-C:HN    0.31000 H
   1.0 ASP-C:CA    0.07000 CT1
  */

  switch (opt) {
  case Sybil:
      AA[p].atom[0].fullatom_type = 40; // Nbb    N
   	  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
   	  AA[p].atom[2].fullatom_type = 10; // CObb   C
   	  AA[p].atom[3].fullatom_type = 56; // OCbb   O
   	  AA[p].atom[4].fullatom_type = 8; // CH2    CB
   	  AA[p].atom[5].fullatom_type = 15; // COO    CG
   	  AA[p].atom[6].fullatom_type = 57; // OOC    OD1
   	  AA[p].atom[7].fullatom_type = 57; // OOC    OD2
  	  AA[p].atom[8].fullatom_type = 66; // HNbb   H
  	  AA[p].atom[9].fullatom_type = 66; // Hapo   HA
  	  AA[p].atom[10].fullatom_type = 66; // Hapo  2HB
  	  AA[p].atom[11].fullatom_type = 66; // Hapo  3HB
  	  AA[p].atom[12].fullatom_type = 66; // Hapo  2HD NEW ADDED BUT UNKNOWN

  	  AA[p].atom[13].fullatom_type = 57; // OXT
  	  AA[p].atom[14].fullatom_type = 66; // Hpol
  	  AA[p].atom[15].fullatom_type = 66; // Hpol

  	  break;

  default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 2; // COO    CG
	  AA[p].atom[6].fullatom_type = 15; // OOC    OD1
	  AA[p].atom[7].fullatom_type = 15; // OOC    OD2
	  AA[p].atom[8].fullatom_type = 25; // HNbb   H
	  AA[p].atom[9].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[10].fullatom_type = 23; // Hapo  2HB
	  AA[p].atom[11].fullatom_type = 23; // Hapo  3HB
	  AA[p].atom[12].fullatom_type = 23; // Hapo  2HD NEW ADDED BUT UNKNOWN

	  AA[p].atom[13].fullatom_type = 20; // OXT
	  AA[p].atom[14].fullatom_type = 22; // Hpol
	  AA[p].atom[15].fullatom_type = 22; // Hpol

	  break;
  }

  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 8; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 9; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--CG
  AA[p].atom[4].bonded_neighbor[2] = 10; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 11; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 3; // CG
  AA[p].atom[5].bonded_neighbor[0] = 4; // CG--CB
  AA[p].atom[5].bonded_neighbor[1] = 6; // CG--OD1
  AA[p].atom[5].bonded_neighbor[2] = 7; // CG--OD2
  AA[p].atom[6].nbonded_neighbors = 1; // OD1
  AA[p].atom[6].bonded_neighbor[0] = 5; // OD1--CG
  AA[p].atom[7].nbonded_neighbors = 2; // OD2
  AA[p].atom[7].bonded_neighbor[0] = 5; // OD2--CG
  AA[p].atom[7].bonded_neighbor[0] = 12; // OD2--2HD NEW ADDED
  AA[p].atom[8].nbonded_neighbors = 1; // H
  AA[p].atom[8].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[9].nbonded_neighbors = 1; // HA
  AA[p].atom[9].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[10].nbonded_neighbors = 1; //2HB
  AA[p].atom[10].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[11].nbonded_neighbors = 1; //3HB
  AA[p].atom[11].bonded_neighbor[0] = 4; //3HB--CB
	//NEW ADDED
  AA[p].atom[12].nbonded_neighbors = 1; //3HB
  AA[p].atom[12].bonded_neighbor[0] = 7; //2HD--OD2


  AA[p].atom[13].nbonded_neighbors = 1; // OXT
  AA[p].atom[13].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[14].nbonded_neighbors = 1; // 1H
  AA[p].atom[14].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[15].nbonded_neighbors = 1; // 2H
  AA[p].atom[15].bonded_neighbor[0] = 0; // H-N



  // A REVISAR
  for ( k = 0; k < 9; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 9; k < 12; k++ )
    AA[p].atom[k].hastemplate = true;

  //bk   template for building  HA
  AA[p].atom[9].ta[0] = 1; //  CA
  AA[p].atom[9].ta[1] = 0; //  N
  AA[p].atom[9].ta[2] = 2; //  C

  //bk   template for building 2HB
  AA[p].atom[10].ta[0] = 4; //  CB
  AA[p].atom[10].ta[1] = 1; //  CA
  AA[p].atom[10].ta[2] = 5; //  CG

  //bk   template for building 3HB
  AA[p].atom[11].ta[0] = 4; //  CB
  AA[p].atom[11].ta[1] = 1; //  CA
  AA[p].atom[11].ta[2] = 5; //  CG



  //bk   chi angles required to build atoms ASP
  //bk   chi angles needed for building  CG
  AA[p].chi_required[0] [5] = true;

  //bk   chi angles needed for building  OD1
  AA[p].chi_required[0] [6] = true;
  AA[p].chi_required[1] [6] = true;

  //bk   chi angles needed for building  OD2
  AA[p].chi_required[0] [7] = true;
  AA[p].chi_required[1] [7] = true;

  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [10] = true;

  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [11] = true;



  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   CG

  //bk   four atoms that define chi angle  2
  AA[p].chi_atoms[1] [0] = 1; //   CA
  AA[p].chi_atoms[1] [1] = 4; //   CB
  AA[p].chi_atoms[1] [2] = 5; //   CG
  AA[p].chi_atoms[1] [3] = 6; //   OD1

//FIN A REVISAR

  AA[p].atom[0].icoor[0] = 0.0000; //   N
  AA[p].atom[0].icoor[1] = 5.2000; //   N
  AA[p].atom[0].icoor[2] = 9.6580; //   N
  AA[p].atom[1].icoor[0] = 0.0000; //   CA
  AA[p].atom[1].icoor[1] = 6.4570; //   CA
  AA[p].atom[1].icoor[2] = 10.3940; //   CA
  AA[p].atom[2].icoor[0] = 0.0000; //   C
  AA[p].atom[2].icoor[1] = 6.2160; //   C
  AA[p].atom[2].icoor[2] = 11.8980; //   C
  AA[p].atom[3].icoor[0] = 0.0010; //   O
  AA[p].atom[3].icoor[1] = 5.0720; //   O
  AA[p].atom[3].icoor[2] = 12.3530; //   O
  AA[p].atom[4].icoor[0] = -1.2080; //   CB
  AA[p].atom[4].icoor[1] = 7.3100; //   CB
  AA[p].atom[4].icoor[2] = 9.9990; //   CB
  AA[p].atom[5].icoor[0] = -1.2950; //   CG
  AA[p].atom[5].icoor[1] = 8.6490; //   CG
  AA[p].atom[5].icoor[2] = 10.7190; //   CG
  AA[p].atom[6].icoor[0] = -2.2210; //   OD1
  AA[p].atom[6].icoor[1] = 9.3810; //   OD1
  AA[p].atom[6].icoor[2] = 10.4610; //   OD1
  AA[p].atom[7].icoor[0] = -0.5300; //   OD2
  AA[p].atom[7].icoor[1] = 8.8620; //   OD2
  AA[p].atom[7].icoor[2] = 11.6290; //   OD2
  AA[p].atom[8].icoor[0] = -0.0400; //   H
  AA[p].atom[8].icoor[1] = 4.3370; //   H
  AA[p].atom[8].icoor[2] = 10.1830; //   H
  AA[p].atom[9].icoor[0] = 0.9090; //   HA
  AA[p].atom[9].icoor[1] = 7.0160; //   HA
  AA[p].atom[9].icoor[2] = 10.1700; //   HA
  AA[p].atom[10].icoor[0] = -1.2940; //   2HB
  AA[p].atom[10].icoor[1] = 7.4700; //   2HB
  AA[p].atom[10].icoor[2] = 8.9240; //   2HB
  AA[p].atom[11].icoor[0] = -2.0180; //   3HB
  AA[p].atom[11].icoor[1] = 6.6630; //   3HB
  AA[p].atom[11].icoor[2] = 10.3370; //   3HB
// NEW ADDED BUT UNKNOWN
  AA[p].atom[12].icoor[0] = 0; //   2HD
  AA[p].atom[12].icoor[1] = 0; //   2HD
  AA[p].atom[12].icoor[2] = 0; //   2HD

  // atom number for backbone HN
  AA[p].HNpos=8;
  // atom number for backbone HA
  AA[p].HApos=9;

  // number of polar hydrogens
  AA[p].nH_polar_complete=2; //NEW ADDED
  // atom numbers for polar H
  AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
  AA[p].Hpos_polar_complete[0]=8;
  AA[p].Hpos_polar_complete[1]=12; //NEW ADDED

  // number of aromatic hydrogens
  AA[p].nH_aromatic_complete=0;

  // number of apolar hydrogens
  AA[p].nH_apolar_complete=3;
  //atom number for apolar hydrogens
  AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
  AA[p].Hpos_apolar_complete[0]=9;
  AA[p].Hpos_apolar_complete[1]=10;
  AA[p].Hpos_apolar_complete[2]=11;


  //A REVISAR
  // number of acceptors
  AA[p].nacceptors=3;
  //acceptor information
  AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
  (AA[p].aBase[1])[0]=6;(AA[p].aBase[1])[1]=5;
  (AA[p].aBase[2])[0]=7;(AA[p].aBase[2])[1]=5;
  (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;
  (AA[p].accpt_pos[1])[0]=1;(AA[p].accpt_pos[1])[1]=6;
  (AA[p].accpt_pos[2])[0]=2;(AA[p].accpt_pos[2])[1]=7;

  //Number of Hydrogens connection
  AA[p].nH_hydrogen_connexions=5; //NEW ADDED
  //atoms hydrogens are connected too
  AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
  AA[p].Hbase[0][0]=8; AA[p].Hbase[0][1]=0;
  AA[p].Hbase[1][0]=9; AA[p].Hbase[1][1]=1;
  AA[p].Hbase[2][0]=10; AA[p].Hbase[2][1]=4;
  AA[p].Hbase[3][0]=11; AA[p].Hbase[3][1]=4;
  AA[p].Hbase[4][0]=12; AA[p].Hbase[4][1]=7; // NEW ADDED OD2-2HD


  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=8;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=9;

  AA[p].atom[4].numHydrogens_atm=2;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=10;
  AA[p].atom[4].hydrogens_atm[0]=11;

  //NEW ADDED
  AA[p].atom[7].numHydrogens_atm=1;
  AA[p].atom[7].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[7].hydrogens_atm[0]=12;


}



//FIN DE AMINOACIDO ALANINA
/**
* Initializes SS-bounded Cysteine
*
* @param opt: Rosseta or ICM
*/
void init_CYX(Convention opt)
{
  //AMINOACIDO SS-bounded CYSTEINA
  int k, p =CYX;

  strcpy( AA[p].aa_name3, "CYX" );
  AA[p].aa_name1 = 'C';
  AA[p].mass = 103.146;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = false;
  AA[p].aa_is_nonpolar = true;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = false;

  AA[p].natoms = 10; //NEW ADDED
  AA[p].nheavyatoms = 6;
  AA[p].nchi = 1;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );


  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;

  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

  AA[p].chi_types = ( int * ) new int[AA[p].nchi];
  AA[p].chi_types[0] = 2;

  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " SG " );
  strcpy( AA[p].atom[6].atom_name, " H  " );
  strcpy( AA[p].atom[7].atom_name, " HA " );
  strcpy( AA[p].atom[8].atom_name, "2HB " );
  strcpy( AA[p].atom[9].atom_name, "3HB " );

  strcpy( AA[p].atom[10].atom_name," OXT" );
  strcpy( AA[p].atom[11].atom_name,"2H  " );
  strcpy( AA[p].atom[12].atom_name,"3H  " );

  switch(opt)
  {
  //MODIFIED BY DISU
  case Rosseta:
  case Sybil:
	AA[p].atom[0].charge = -0.47; //   N  <--- Nt(-0.30)
    AA[p].atom[1].charge = 0.07; //   CA  <--- Nt(0.21)
    AA[p].atom[2].charge = 0.51; //   C  <--- Ct(0.34)
    AA[p].atom[3].charge = -0.51; //   O  <--- Ct(-0.67)
    AA[p].atom[4].charge = -0.10; //   CB
    AA[p].atom[5].charge = -0.08; //   SG
    AA[p].atom[6].charge = 0.31; //   H  <--- Nt(0.33)
    AA[p].atom[7].charge = 0.09; //   HA  <--- Nt(0.10)
    AA[p].atom[8].charge = 0.09; //   2HB
    AA[p].atom[9].charge = 0.09; //   3HB

    AA[p].atom[10].charge = -0.67; //   0XT
    AA[p].atom[11].charge = 0.33; //  1H
    AA[p].atom[12].charge = 0.33; //  2H
    break;
  //A REVISAR
  case ICM:
    AA[p].atom[0].charge = -0.356; //   N
    AA[p].atom[1].charge = 0.064; //   CA
    AA[p].atom[2].charge = 0.45; //   C
    AA[p].atom[3].charge = -0.384; //   O
    AA[p].atom[4].charge = -0.105; //   CB
    AA[p].atom[5].charge = 0.015; //   SG
    AA[p].atom[6].charge = 0.176; //   H
    AA[p].atom[7].charge = 0.02; //   HA
    AA[p].atom[8].charge = 0.055; //   2HB
    AA[p].atom[9].charge = 0.055; //   3HB

    AA[p].atom[10].charge = -0.67; //   0XT
    AA[p].atom[11].charge = 0.33; //  1H
    AA[p].atom[12].charge = 0.33; //  2H
    break;
	default:
	  break;
  }

//CARGAS INCORRECTAS
/*AA[p].atom[0].charge = 0.5; //   N
AA[p].atom[1].charge = 0.5; //   CA
AA[p].atom[2].charge = 0.5; //   C
AA[p].atom[3].charge = 0.5; //   O
AA[p].atom[4].charge = 0.5; //   CB
AA[p].atom[5].charge = 0.5; //   SG
AA[p].atom[6].charge = 0.5; //   H
AA[p].atom[7].charge = 0.5; //   HA
AA[p].atom[8].charge = 0.5; //   2HB
AA[p].atom[9].charge = 0.5; //   3HB
AA[p].atom[10].charge = 0.5; //   HG

AA[p].atom[11].charge = 0.5; //   0XT
AA[p].atom[12].charge = 0.5; //  1H
AA[p].atom[13].charge = 0.5; //  2H
*/

  /*   Nt/Ct - Charge Table

   1.0 CYS-N:2H    0.33000 HC
   1.0 CYS-N:HA    0.10000 HB
   1.0 CYS-N:3H    0.33000 HC
   1.0 CYS-N:CB   -0.11000 CT2
   1.0 CYS-N:HB1   0.09000 HA
   1.0 CYS-N:1HB   0.09000 HA
   1.0 CYS-N:SG   -0.23000 S
   1.0 CYS-N:HG1   0.16000 HS
   1.0 CYS-N:C     0.51000 C
   1.0 CYS-N:O    -0.51000 O
   1.0 CYS-N:N    -0.30000 NH3
   1.0 CYS-N:1H    0.33000 HC
   1.0 CYS-N:CA    0.21000 CT1

   1.0 CYS-C:HA    0.09000 HB
   1.0 CYS-C:CB   -0.11000 CT2
   1.0 CYS-C:HB1   0.09000 HA
   1.0 CYS-C:1HB   0.09000 HA
   1.0 CYS-C:SG   -0.23000 S
   1.0 CYS-C:HG1   0.16000 HS
   1.0 CYS-C:C     0.34000 CC
   1.0 CYS-C:O    -0.67000 OC
   1.0 CYS-C:OXT  -0.67000 OC
   1.0 CYS-C:N    -0.47000 NH1
   1.0 CYS-C:HN    0.31000 H
   1.0 CYS-C:CA    0.07000 CT1
*/
  switch (opt) {
  case Sybil:

	  AA[p].atom[0].fullatom_type = 40; // Nbb    N
	  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 10; // CObb   C
      AA[p].atom[3].fullatom_type = 56; // OCbb   O
	  AA[p].atom[4].fullatom_type = 8; // CH2    CB
	  AA[p].atom[5].fullatom_type = 63; // S      SG
	  AA[p].atom[6].fullatom_type = 66; // HNbb   H
	  AA[p].atom[7].fullatom_type = 66; // Hapo   HA
	  AA[p].atom[8].fullatom_type = 66; // Hapo   2HB
	  AA[p].atom[9].fullatom_type = 66; // Hapo   3HB

	  AA[p].atom[10].fullatom_type = 57; // OXT
	  AA[p].atom[11].fullatom_type = 66; // Hpol
	  AA[p].atom[12].fullatom_type = 66; // Hpol
	  break;
  default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 16; // S      SG
	  AA[p].atom[6].fullatom_type = 25; // HNbb   H
	  AA[p].atom[7].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[8].fullatom_type = 23; // Hapo   2HB
	  AA[p].atom[9].fullatom_type = 23; // Hapo   3HB

	  AA[p].atom[10].fullatom_type = 20; // OXT
	  AA[p].atom[11].fullatom_type = 22; // Hpol
	  AA[p].atom[12].fullatom_type = 22; // Hpol
	  break;
  }

  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 6; // N--H (3H) HNbb
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 7; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--SG
  AA[p].atom[4].bonded_neighbor[2] = 8; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 9; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 1; // SG
  AA[p].atom[5].bonded_neighbor[0] = 4; // SG--CB
  AA[p].atom[6].nbonded_neighbors = 1; // H
  AA[p].atom[6].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[7].nbonded_neighbors = 1; // HA
  AA[p].atom[7].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[8].nbonded_neighbors = 1; //2HB
  AA[p].atom[8].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[9].nbonded_neighbors = 1; //3HB
  AA[p].atom[9].bonded_neighbor[0] = 4; //3HB--CB

  AA[p].atom[10].nbonded_neighbors = 1; // OXT
  AA[p].atom[10].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[11].nbonded_neighbors = 1; // 1H
  AA[p].atom[11].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[12].nbonded_neighbors = 1; // 2H
  AA[p].atom[12].bonded_neighbor[0] = 0; // H-N



  for ( k = 0; k < 7; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 7; k < 11; k++ )
    AA[p].atom[k].hastemplate = true;

  //bk   template for building  HA
  AA[p].atom[7].ta[0] = 1; //  CA
  AA[p].atom[7].ta[1] = 0; //  N
  AA[p].atom[7].ta[2] = 2; //  C

  //bk   template for building 2HB
  AA[p].atom[8].ta[0] = 4; //  CB
  AA[p].atom[8].ta[1] = 1; //  CA
  AA[p].atom[8].ta[2] = 5; //  SG

  //bk   template for building 3HB
  AA[p].atom[9].ta[0] = 4; //  CB
  AA[p].atom[9].ta[1] = 1; //  CA
  AA[p].atom[9].ta[2] = 5; //  SG


  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   SG

  //bk   chi angles required to build atoms CYS
  //bk   chi angles needed for building  SG
  AA[p].chi_required[0] [5] = true;

  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [8] = true;

  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [9] = true;


  //bk   template coordinates for the amino acid
  AA[p].atom[0].icoor[0] = 0.0000; //   N
  AA[p].atom[0].icoor[1] = 3.8070; //   N
  AA[p].atom[0].icoor[2] = 6.2240; //   N
  AA[p].atom[1].icoor[0] = 0.0000; //   CA
  AA[p].atom[1].icoor[1] = 3.7830; //   CA
  AA[p].atom[1].icoor[2] = 7.6880; //   CA
  AA[p].atom[2].icoor[0] = 0.0000; //   C
  AA[p].atom[2].icoor[1] = 5.1680; //   C
  AA[p].atom[2].icoor[2] = 8.3290; //   C
  AA[p].atom[3].icoor[0] = 0.0050; //   O
  AA[p].atom[3].icoor[1] = 6.1860; //   O
  AA[p].atom[3].icoor[2] = 7.6360; //   O
  AA[p].atom[4].icoor[0] = 1.3120; //   CB
  AA[p].atom[4].icoor[1] = 3.0690; //   CB
  AA[p].atom[4].icoor[2] = 8.0140; //   CB
  AA[p].atom[5].icoor[0] = 1.6350; //   SG
  AA[p].atom[5].icoor[1] = 2.8650; //   SG
  AA[p].atom[5].icoor[2] = 9.7820; //   SG
  AA[p].atom[6].icoor[0] = -0.0250; //   H
  AA[p].atom[6].icoor[1] = 4.7040; //   H
  AA[p].atom[6].icoor[2] = 5.7600; //   H
  AA[p].atom[7].icoor[0] = -0.8240; //   HA
  AA[p].atom[7].icoor[1] = 3.2090; //   HA
  AA[p].atom[7].icoor[2] = 8.1120; //   HA
  AA[p].atom[8].icoor[0] = 1.3110; //   2HB
  AA[p].atom[8].icoor[1] = 2.0650; //   2HB
  AA[p].atom[8].icoor[2] = 7.5890; //   2HB
  AA[p].atom[9].icoor[0] = 2.1560; //   3HB
  AA[p].atom[9].icoor[1] = 3.6320; //   3HB
  AA[p].atom[9].icoor[2] = 7.6160; //   3HB

  // atom number for backbone HN
  AA[p].HNpos=6;
  // atom number for backbone HA
  AA[p].HApos=7;

  // number of polar hydrogens
  AA[p].nH_polar_complete=1;
  // atom numbers for polar H
  AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
  AA[p].Hpos_polar_complete[0]=6;

  // number of aromatic hydrogens
  AA[p].nH_aromatic_complete=0;

  // number of apolar hydrogens
  AA[p].nH_apolar_complete=3;
  //atom number for apolar hydrogens
  AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
  AA[p].Hpos_apolar_complete[0]=7;
  AA[p].Hpos_apolar_complete[1]=8;
  AA[p].Hpos_apolar_complete[2]=9;

  //A REVISAR
  // number of acceptors
  AA[p].nacceptors=1;
  //acceptor information
  AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
  (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;

  //Number of Hydrogens connection
  AA[p].nH_hydrogen_connexions=4;
  //atoms hydrogens are connected too
  AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
  AA[p].Hbase[0][0]=6; AA[p].Hbase[0][1]=0;
  AA[p].Hbase[1][0]=7; AA[p].Hbase[1][1]=1;
  AA[p].Hbase[2][0]=8; AA[p].Hbase[2][1]=4;
  AA[p].Hbase[3][0]=9; AA[p].Hbase[3][1]=4;

  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=6;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=7;

  AA[p].atom[4].numHydrogens_atm=2;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=8;
  AA[p].atom[4].hydrogens_atm[1]=9;


}

//FIN DE AMINOACIDO SS-bonded CYSTEINA

/**
* Initializes Cysteine N
*
* @param opt: Rosseta or ICM
*/
//NO SE USA
void init_CYM(Convention opt)
{
  //AMINOACIDO CYSTEINA NEGATIVA
  int k, p =CYM;

  strcpy( AA[p].aa_name3, "CYM" );
  AA[p].aa_name1 = 'C';
  AA[p].mass = 103.146;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = false;
  AA[p].aa_is_nonpolar = true;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = false;

  AA[p].natoms = 10; //NEW ADDED
  AA[p].nheavyatoms = 6;
  AA[p].nchi = 1;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );


  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;

  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

  AA[p].chi_types = ( int * ) new int[AA[p].nchi];
  AA[p].chi_types[0] = 2;

  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " SG " );
  strcpy( AA[p].atom[6].atom_name, " H  " );
  strcpy( AA[p].atom[7].atom_name, " HA " );
  strcpy( AA[p].atom[8].atom_name, "2HB " );
  strcpy( AA[p].atom[9].atom_name, "3HB " );

  strcpy( AA[p].atom[10].atom_name," OXT" );
  strcpy( AA[p].atom[11].atom_name,"2H  " );
  strcpy( AA[p].atom[12].atom_name,"3H  " );
// OJO COMENTADO MOMENTANEAMENTE
  switch(opt)
  {
  // A REVISAR
  case Rosseta:
  case Sybil:
	AA[p].atom[0].charge = -0.47; //   N  <--- Nt(-0.30)
    AA[p].atom[1].charge = 0.07; //   CA  <--- Nt(0.21)
    AA[p].atom[2].charge = 0.51; //   C  <--- Ct(0.34)
    AA[p].atom[3].charge = -0.51; //   O  <--- Ct(-0.67)
    AA[p].atom[4].charge = -0.11; //   CB
    AA[p].atom[5].charge = -0.23; //   SG
    AA[p].atom[6].charge = 0.31; //   H  <--- Nt(0.33)
    AA[p].atom[7].charge = 0.09; //   HA  <--- Nt(0.10)
    AA[p].atom[8].charge = 0.09; //   2HB
    AA[p].atom[9].charge = 0.09; //   3HB

    AA[p].atom[10].charge = -0.67; //   0XT
    AA[p].atom[11].charge = 0.33; //  1H
    AA[p].atom[12].charge = 0.33; //  2H
    break;
	// A REVISAR
  case ICM:
    AA[p].atom[0].charge = -0.356; //   N
    AA[p].atom[1].charge = 0.064; //   CA
    AA[p].atom[2].charge = 0.45; //   C
    AA[p].atom[3].charge = -0.384; //   O
    AA[p].atom[4].charge = -0.105; //   CB
    AA[p].atom[5].charge = 0.015; //   SG
    AA[p].atom[6].charge = 0.176; //   H
    AA[p].atom[7].charge = 0.02; //   HA
    AA[p].atom[8].charge = 0.055; //   2HB
    AA[p].atom[9].charge = 0.055; //   3HB

    AA[p].atom[10].charge = -0.67; //   0XT
    AA[p].atom[11].charge = 0.33; //  1H
    AA[p].atom[12].charge = 0.33; //  2H
    break;
  default:
    AA[p].atom[0].charge = -0.356; //   N
    AA[p].atom[1].charge = 0.064; //   CA
    AA[p].atom[2].charge = 0.45; //   C
    AA[p].atom[3].charge = -0.384; //   O
    AA[p].atom[4].charge = -0.105; //   CB
    AA[p].atom[5].charge = 0.015; //   SG
    AA[p].atom[6].charge = 0.176; //   H
    AA[p].atom[7].charge = 0.02; //   HA
    AA[p].atom[8].charge = 0.055; //   2HB
    AA[p].atom[9].charge = 0.055; //   3HB

    AA[p].atom[10].charge = -0.67; //   0XT
    AA[p].atom[11].charge = 0.33; //  1H
    AA[p].atom[12].charge = 0.33; //  2H
    break;
  }

//CARGAS INCORRECTAS
/*AA[p].atom[0].charge = 0.5; //   N
AA[p].atom[1].charge = 0.5; //   CA
AA[p].atom[2].charge = 0.5; //   C
AA[p].atom[3].charge = 0.5; //   O
AA[p].atom[4].charge = 0.5; //   CB
AA[p].atom[5].charge = 0.5; //   SG
AA[p].atom[6].charge = 0.5; //   H
AA[p].atom[7].charge = 0.5; //   HA
AA[p].atom[8].charge = 0.5; //   2HB
AA[p].atom[9].charge = 0.5; //   3HB
AA[p].atom[10].charge = 0.5; //   HG

AA[p].atom[11].charge = 0.5; //   0XT
AA[p].atom[12].charge = 0.5; //  1H
AA[p].atom[13].charge = 0.5; //  2H
*/

  /*   Nt/Ct - Charge Table

   1.0 CYS-N:2H    0.33000 HC
   1.0 CYS-N:HA    0.10000 HB
   1.0 CYS-N:3H    0.33000 HC
   1.0 CYS-N:CB   -0.11000 CT2
   1.0 CYS-N:HB1   0.09000 HA
   1.0 CYS-N:1HB   0.09000 HA
   1.0 CYS-N:SG   -0.23000 S
   1.0 CYS-N:HG1   0.16000 HS
   1.0 CYS-N:C     0.51000 C
   1.0 CYS-N:O    -0.51000 O
   1.0 CYS-N:N    -0.30000 NH3
   1.0 CYS-N:1H    0.33000 HC
   1.0 CYS-N:CA    0.21000 CT1

   1.0 CYS-C:HA    0.09000 HB
   1.0 CYS-C:CB   -0.11000 CT2
   1.0 CYS-C:HB1   0.09000 HA
   1.0 CYS-C:1HB   0.09000 HA
   1.0 CYS-C:SG   -0.23000 S
   1.0 CYS-C:HG1   0.16000 HS
   1.0 CYS-C:C     0.34000 CC
   1.0 CYS-C:O    -0.67000 OC
   1.0 CYS-C:OXT  -0.67000 OC
   1.0 CYS-C:N    -0.47000 NH1
   1.0 CYS-C:HN    0.31000 H
   1.0 CYS-C:CA    0.07000 CT1
*/
  switch (opt) {
  case Sybil:
	  AA[p].atom[0].fullatom_type = 40; // Nbb    N
	  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 10; // CObb   C
      AA[p].atom[3].fullatom_type = 56; // OCbb   O
	  AA[p].atom[4].fullatom_type = 8; // CH2    CB
	  AA[p].atom[5].fullatom_type = 63; // S      SG
 	  AA[p].atom[6].fullatom_type = 66; // HNbb   H
 	  AA[p].atom[7].fullatom_type = 66; // Hapo   HA
 	  AA[p].atom[8].fullatom_type = 66; // Hapo   2HB
 	  AA[p].atom[9].fullatom_type = 66; // Hapo   3HB
 	  AA[p].atom[10].fullatom_type = 57; // OXT
 	  AA[p].atom[11].fullatom_type = 66; // Hpol
 	  AA[p].atom[12].fullatom_type = 66; // Hpol
 	  break;
  default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 16; // S      SG
	  AA[p].atom[6].fullatom_type = 25; // HNbb   H
	  AA[p].atom[7].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[8].fullatom_type = 23; // Hapo   2HB
	  AA[p].atom[9].fullatom_type = 23; // Hapo   3HB

	  AA[p].atom[10].fullatom_type = 20; // OXT
	  AA[p].atom[11].fullatom_type = 22; // Hpol
	  AA[p].atom[12].fullatom_type = 22; // Hpol
	  break;
  }

  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 6; // N--H (3H) HNbb
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 7; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--SG
  AA[p].atom[4].bonded_neighbor[2] = 8; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 9; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 1; // SG
  AA[p].atom[5].bonded_neighbor[0] = 4; // SG--CB
  AA[p].atom[6].nbonded_neighbors = 1; // H
  AA[p].atom[6].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[7].nbonded_neighbors = 1; // HA
  AA[p].atom[7].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[8].nbonded_neighbors = 1; //2HB
  AA[p].atom[8].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[9].nbonded_neighbors = 1; //3HB
  AA[p].atom[9].bonded_neighbor[0] = 4; //3HB--CB

  AA[p].atom[10].nbonded_neighbors = 1; // OXT
  AA[p].atom[10].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[11].nbonded_neighbors = 1; // 1H
  AA[p].atom[11].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[12].nbonded_neighbors = 1; // 2H
  AA[p].atom[12].bonded_neighbor[0] = 0; // H-N



  for ( k = 0; k < 7; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 7; k < 11; k++ )
    AA[p].atom[k].hastemplate = true;

  //bk   template for building  HA
  AA[p].atom[7].ta[0] = 1; //  CA
  AA[p].atom[7].ta[1] = 0; //  N
  AA[p].atom[7].ta[2] = 2; //  C

  //bk   template for building 2HB
  AA[p].atom[8].ta[0] = 4; //  CB
  AA[p].atom[8].ta[1] = 1; //  CA
  AA[p].atom[8].ta[2] = 5; //  SG

  //bk   template for building 3HB
  AA[p].atom[9].ta[0] = 4; //  CB
  AA[p].atom[9].ta[1] = 1; //  CA
  AA[p].atom[9].ta[2] = 5; //  SG


  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   SG

  //bk   chi angles required to build atoms CYS
  //bk   chi angles needed for building  SG
  AA[p].chi_required[0] [5] = true;

  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [8] = true;

  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [9] = true;



  //bk   template coordinates for the amino acid
  AA[p].atom[0].icoor[0] = 0.0000; //   N
  AA[p].atom[0].icoor[1] = 3.8070; //   N
  AA[p].atom[0].icoor[2] = 6.2240; //   N
  AA[p].atom[1].icoor[0] = 0.0000; //   CA
  AA[p].atom[1].icoor[1] = 3.7830; //   CA
  AA[p].atom[1].icoor[2] = 7.6880; //   CA
  AA[p].atom[2].icoor[0] = 0.0000; //   C
  AA[p].atom[2].icoor[1] = 5.1680; //   C
  AA[p].atom[2].icoor[2] = 8.3290; //   C
  AA[p].atom[3].icoor[0] = 0.0050; //   O
  AA[p].atom[3].icoor[1] = 6.1860; //   O
  AA[p].atom[3].icoor[2] = 7.6360; //   O
  AA[p].atom[4].icoor[0] = 1.3120; //   CB
  AA[p].atom[4].icoor[1] = 3.0690; //   CB
  AA[p].atom[4].icoor[2] = 8.0140; //   CB
  AA[p].atom[5].icoor[0] = 1.6350; //   SG
  AA[p].atom[5].icoor[1] = 2.8650; //   SG
  AA[p].atom[5].icoor[2] = 9.7820; //   SG
  AA[p].atom[6].icoor[0] = -0.0250; //   H
  AA[p].atom[6].icoor[1] = 4.7040; //   H
  AA[p].atom[6].icoor[2] = 5.7600; //   H
  AA[p].atom[7].icoor[0] = -0.8240; //   HA
  AA[p].atom[7].icoor[1] = 3.2090; //   HA
  AA[p].atom[7].icoor[2] = 8.1120; //   HA
  AA[p].atom[8].icoor[0] = 1.3110; //   2HB
  AA[p].atom[8].icoor[1] = 2.0650; //   2HB
  AA[p].atom[8].icoor[2] = 7.5890; //   2HB
  AA[p].atom[9].icoor[0] = 2.1560; //   3HB
  AA[p].atom[9].icoor[1] = 3.6320; //   3HB
  AA[p].atom[9].icoor[2] = 7.6160; //   3HB

  // atom number for backbone HN
  AA[p].HNpos=6;
  // atom number for backbone HA
  AA[p].HApos=7;

  // number of polar hydrogens
  AA[p].nH_polar_complete=1;
  // atom numbers for polar H
  AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
  AA[p].Hpos_polar_complete[0]=6;

  // number of aromatic hydrogens
  AA[p].nH_aromatic_complete=0;

  // number of apolar hydrogens
  AA[p].nH_apolar_complete=3;
  //atom number for apolar hydrogens
  AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
  AA[p].Hpos_apolar_complete[0]=7;
  AA[p].Hpos_apolar_complete[1]=8;
  AA[p].Hpos_apolar_complete[2]=9;

  //A REVISAR
  // number of acceptors
  AA[p].nacceptors=1;
  //acceptor information
  AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
  (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;

  //Number of Hydrogens connection
  AA[p].nH_hydrogen_connexions=4;
  //atoms hydrogens are connected too
  AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
  AA[p].Hbase[0][0]=6; AA[p].Hbase[0][1]=0;
  AA[p].Hbase[1][0]=7; AA[p].Hbase[1][1]=1;
  AA[p].Hbase[2][0]=8; AA[p].Hbase[2][1]=4;
  AA[p].Hbase[3][0]=9; AA[p].Hbase[3][1]=4;

  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=6;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=7;

  AA[p].atom[4].numHydrogens_atm=2;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=8;
  AA[p].atom[4].hydrogens_atm[1]=9;

}
//FIN DE AMINOACIDO CYSTEINA NEGATIVA

/**
* Initializes Glutamate NEUTRAL
*
* @param opt: Rosseta or ICM
*/
void init_GLH(Convention opt)
{

  //AMINOACIDO GLUTAMATO NEUTRAL
  int k, p = GLH;


  strcpy( AA[p].aa_name3, "GLH" );
  AA[p].aa_name1 = 'E';
  AA[p].mass = 128.11;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = true;
  AA[p].aa_is_nonpolar = false;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = true;

  AA[p].natoms = 16; // NEW ADDED
  AA[p].nheavyatoms = 9;
  AA[p].nchi = 3;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );

  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;

  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

// Bug 23/11/2009
//AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );
AA[p].chi_types = ( int * ) new int[AA[p].nchi];
AA[p].chi_types[0] = 2;
AA[p].chi_types[1] = 2;
AA[p].chi_types[2] = 4;


  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " CG " );
  strcpy( AA[p].atom[6].atom_name, " CD " );
  strcpy( AA[p].atom[7].atom_name, " OE1" );
  strcpy( AA[p].atom[8].atom_name, " OE2" );
  strcpy( AA[p].atom[9].atom_name, " H  " );
  strcpy( AA[p].atom[10].atom_name, " HA " );
  strcpy( AA[p].atom[11].atom_name, "2HB " );
  strcpy( AA[p].atom[12].atom_name, "3HB " );
  strcpy( AA[p].atom[13].atom_name, "2HG " );
  strcpy( AA[p].atom[14].atom_name, "3HG " );
  strcpy( AA[p].atom[15].atom_name, "2HE " );//NEW ADDED


  strcpy( AA[p].atom[16].atom_name," OXT" );
  strcpy( AA[p].atom[17].atom_name,"2H  " );
  strcpy( AA[p].atom[18].atom_name,"3H  " );


  switch (opt)
  {
    //NEW ADDED. MODIFIED BY GLUP
    case Rosseta:
    case Sybil:
      AA[p].atom[0].charge = -0.47; //    N
      AA[p].atom[1].charge = 0.07; //    CA
      AA[p].atom[2].charge = 0.51; //    C
      AA[p].atom[3].charge = -0.51; //    O
      AA[p].atom[4].charge = -0.18; //    CB
      AA[p].atom[5].charge = -0.21; //    CG
      AA[p].atom[6].charge = 0.75; //    CD
      AA[p].atom[7].charge = -0.55; //    OE1
      AA[p].atom[8].charge = -0.61; //    OE2
      AA[p].atom[9].charge = 0.31; //    H
      AA[p].atom[10].charge = 0.09; //    HA
      AA[p].atom[11].charge = 0.09; //    2HB
      AA[p].atom[12].charge = 0.09; //    3HB
      AA[p].atom[13].charge = 0.09; //    2HG
      AA[p].atom[14].charge = 0.09; //    3HG
      AA[p].atom[15].charge = 0.09; //    2HE NEW ADDED

      AA[p].atom[16].charge = -0.67; //   0XT
      AA[p].atom[17].charge = 0.33; //  1H
      AA[p].atom[18].charge = 0.33; //  2H
      break;
    //A REVISAR
    case ICM:
      AA[p].atom[0].charge = -0.356; //    N
      AA[p].atom[1].charge = 0.064; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = -0.12; //    CB
      AA[p].atom[5].charge = -0.17; //    CG
      AA[p].atom[6].charge = 0.5; //    CD
      AA[p].atom[7].charge = -0.57; //    OE1
      AA[p].atom[8].charge = -0.57; //    OE2
      AA[p].atom[9].charge = 0.176; //    H
      AA[p].atom[10].charge = 0.02; //    HA
      AA[p].atom[11].charge = 0.02; //    2HB
      AA[p].atom[12].charge = 0.02; //    3HB
      AA[p].atom[13].charge = -0.04; //    2HG
      AA[p].atom[14].charge = -0.04; //    3HG
      AA[p].atom[15].charge = -0.04; //    2HE NEW ADDED

      AA[p].atom[16].charge = -0.67; //   0XT
      AA[p].atom[17].charge = 0.33; //  1H
      AA[p].atom[18].charge = 0.33; //  2H
      break;
    default:
      AA[p].atom[0].charge = -0.356; //    N
      AA[p].atom[1].charge = 0.064; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = -0.12; //    CB
      AA[p].atom[5].charge = -0.17; //    CG
      AA[p].atom[6].charge = 0.5; //    CD
      AA[p].atom[7].charge = -0.57; //    OE1
      AA[p].atom[8].charge = -0.57; //    OE2
      AA[p].atom[9].charge = 0.176; //    H
      AA[p].atom[10].charge = 0.02; //    HA
      AA[p].atom[11].charge = 0.02; //    2HB
      AA[p].atom[12].charge = 0.02; //    3HB
      AA[p].atom[13].charge = -0.04; //    2HG
      AA[p].atom[14].charge = -0.04; //    3HG
      AA[p].atom[15].charge = -0.04; //    2HE NEW ADDED

      AA[p].atom[16].charge = -0.67; //   0XT
      AA[p].atom[17].charge = 0.33; //  1H
      AA[p].atom[18].charge = 0.33; //  2H
      break;
  }
/*
   1.0 GLU-N:2H    0.33000 HC
   1.0 GLU-N:HA    0.10000 HB
   1.0 GLU-N:3H    0.33000 HC
   1.0 GLU-N:CB   -0.18000 CT2
   1.0 GLU-N:HB1   0.09000 HA
   1.0 GLU-N:1HB   0.09000 HA
   1.0 GLU-N:CG   -0.28000 CT2
   1.0 GLU-N:HG1   0.09000 HA
   1.0 GLU-N:1HG   0.09000 HA
   1.0 GLU-N:CD    0.62000 CC
   1.0 GLU-N:OE1  -0.76000 OC
   1.0 GLU-N:OE2  -0.76000 OC
   1.0 GLU-N:C     0.51000 C
   1.0 GLU-N:N    -0.30000 NH3
   1.0 GLU-N:O    -0.51000 O
   1.0 GLU-N:1H    0.33000 HC
   1.0 GLU-N:CA    0.21000 CT1

   1.0 GLU-C:HA    0.09000 HB
   1.0 GLU-C:CB   -0.18000 CT2
   1.0 GLU-C:HB1   0.09000 HA
   1.0 GLU-C:1HB   0.09000 HA
   1.0 GLU-C:CG   -0.28000 CT2
   1.0 GLU-C:HG1   0.09000 HA
   1.0 GLU-C:1HG   0.09000 HA
   1.0 GLU-C:CD    0.62000 CC
   1.0 GLU-C:O    -0.67000 OC
   1.0 GLU-C:OE1  -0.76000 OC
   1.0 GLU-C:OXT  -0.67000 OC
   1.0 GLU-C:OE2  -0.76000 OC
   1.0 GLU-C:C     0.34000 CC
   1.0 GLU-C:N    -0.47000 NH1
   1.0 GLU-C:HN    0.31000 H
   1.0 GLU-C:CA    0.07000 CT1
*/
  switch (opt) {
  case Sybil:
	  AA[p].atom[0].fullatom_type = 40; // Nbb    N
  	  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
  	  AA[p].atom[2].fullatom_type = 10; // CObb   C
  	  AA[p].atom[3].fullatom_type = 56; // OCbb   O
  	  AA[p].atom[4].fullatom_type = 8; // CH2    CB
  	  AA[p].atom[5].fullatom_type = 8; // CH2    CG
  	  AA[p].atom[6].fullatom_type = 11; // COO    CD
  	  AA[p].atom[7].fullatom_type = 57; // OOC    OE1
  	  AA[p].atom[8].fullatom_type = 57; // OOC    OE2
  	  AA[p].atom[9].fullatom_type = 66; // HNbb   H
  	  AA[p].atom[10].fullatom_type = 66; // Hapo   HA
  	  AA[p].atom[11].fullatom_type = 66; // Hapo  2HB
  	  AA[p].atom[12].fullatom_type = 66; // Hapo  3HB
  	  AA[p].atom[13].fullatom_type = 66; // Hapo  2HG
  	  AA[p].atom[14].fullatom_type = 66; // Hapo  3HG
  	  AA[p].atom[15].fullatom_type = 66; // Hapo  2HE NEW ADDED BUT UNKNOWN

  	  AA[p].atom[16].fullatom_type = 57; // OXT
  	  AA[p].atom[17].fullatom_type = 66; // Hpol
  	  AA[p].atom[18].fullatom_type = 66; // Hpol
  	  break;
  default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 4; // CH2    CG
	  AA[p].atom[6].fullatom_type = 2; // COO    CD
	  AA[p].atom[7].fullatom_type = 15; // OOC    OE1
	  AA[p].atom[8].fullatom_type = 15; // OOC    OE2
	  AA[p].atom[9].fullatom_type = 25; // HNbb   H
	  AA[p].atom[10].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[11].fullatom_type = 23; // Hapo  2HB
	  AA[p].atom[12].fullatom_type = 23; // Hapo  3HB
	  AA[p].atom[13].fullatom_type = 23; // Hapo  2HG
	  AA[p].atom[14].fullatom_type = 23; // Hapo  3HG
	  AA[p].atom[15].fullatom_type = 23; // Hapo  2HE NEW ADDED BUT UNKNOWN

	  AA[p].atom[16].fullatom_type = 20; // OXT
	  AA[p].atom[17].fullatom_type = 22; // Hpol
	  AA[p].atom[18].fullatom_type = 22; // Hpol
	  break;
  }


  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 9; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 10; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--CG
  AA[p].atom[4].bonded_neighbor[2] = 11; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 12; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 4; // CG
  AA[p].atom[5].bonded_neighbor[0] = 4; // CG--CB
  AA[p].atom[5].bonded_neighbor[1] = 6; // CG--CD
  AA[p].atom[5].bonded_neighbor[2] = 13; // CG--2HG
  AA[p].atom[5].bonded_neighbor[3] = 14; // CG--3HG
  AA[p].atom[6].nbonded_neighbors = 3; // CD
  AA[p].atom[6].bonded_neighbor[0] = 5; // CD--CG
  AA[p].atom[6].bonded_neighbor[1] = 7; // CD--OE1
  AA[p].atom[6].bonded_neighbor[2] = 8; // CD--OE2
  AA[p].atom[7].nbonded_neighbors = 1; // OE1
  AA[p].atom[7].bonded_neighbor[0] = 6; // OE1--CD
  AA[p].atom[8].nbonded_neighbors = 2; // OE2
  AA[p].atom[8].bonded_neighbor[0] = 6; // OE2--CD
  AA[p].atom[8].bonded_neighbor[0] = 15; // OE2--HE2 NEW ADDED
  AA[p].atom[9].nbonded_neighbors = 1; // H
  AA[p].atom[9].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[10].nbonded_neighbors = 1; // HA
  AA[p].atom[10].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[11].nbonded_neighbors = 1; //2HB
  AA[p].atom[11].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[12].nbonded_neighbors = 1; //3HB
  AA[p].atom[12].bonded_neighbor[0] = 4; //3HB--CB
  AA[p].atom[13].nbonded_neighbors = 1; //2HG
  AA[p].atom[13].bonded_neighbor[0] = 5; //2HG--CG
  AA[p].atom[14].nbonded_neighbors = 1; //3HG
  AA[p].atom[14].bonded_neighbor[0] = 5; //3HG--CG
  //NEW ADDED
  AA[p].atom[15].nbonded_neighbors = 1; //2HE
  AA[p].atom[15].bonded_neighbor[0] = 8; //2HE--OE2


  AA[p].atom[16].nbonded_neighbors = 1; // OXT
  AA[p].atom[16].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[17].nbonded_neighbors = 1; // 1H
  AA[p].atom[17].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[18].nbonded_neighbors = 1; // 2H
  AA[p].atom[18].bonded_neighbor[0] = 0; // H-N


  //A REVISAR CHIS
  for ( k = 0; k < 10; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 10; k < 15; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template for building  HA
  AA[p].atom[10].ta[0] = 1; //  CA
  AA[p].atom[10].ta[1] = 0; //  N
  AA[p].atom[10].ta[2] = 2; //  C

  //bk   template for building 2HB
  AA[p].atom[11].ta[0] = 4; //  CB
  AA[p].atom[11].ta[1] = 1; //  CA
  AA[p].atom[11].ta[2] = 5; //  CG

  //bk   template for building 3HB
  AA[p].atom[12].ta[0] = 4; //  CB
  AA[p].atom[12].ta[1] = 1; //  CA
  AA[p].atom[12].ta[2] = 5; //  CG

  //bk   template for building 2HG
  AA[p].atom[13].ta[0] = 5; //  CG
  AA[p].atom[13].ta[1] = 4; //  CB
  AA[p].atom[13].ta[2] = 6; //  CD

  //bk   template for building 3HG
  AA[p].atom[14].ta[0] = 5; //  CG
  AA[p].atom[14].ta[1] = 4; //  CB
  AA[p].atom[14].ta[2] = 6; //  CD


  //bk   chi angles required to build atoms GLU
  //bk   chi angles needed for building  CG
  AA[p].chi_required[0] [5] = true;

  //bk   chi angles needed for building  CD
  AA[p].chi_required[0] [6] = true;
  AA[p].chi_required[1] [6] = true;

  //bk   chi angles needed for building  OE1
  AA[p].chi_required[0] [7] = true;
  AA[p].chi_required[1] [7] = true;
  AA[p].chi_required[2] [7] = true;

  //bk   chi angles needed for building  OE2
  AA[p].chi_required[0] [8] = true;
  AA[p].chi_required[1] [8] = true;
  AA[p].chi_required[2] [8] = true;

  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [11] = true;

  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [12] = true;

  //bk   chi angles needed for building 2HG
  AA[p].chi_required[0] [13] = true;
  AA[p].chi_required[1] [13] = true;

  //bk   chi angles needed for building 3HG
  AA[p].chi_required[0] [14] = true;
  AA[p].chi_required[1] [14] = true;




  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   CG
  //bk   four atoms that define chi angle  2
  AA[p].chi_atoms[1] [0] = 1; //   CA
  AA[p].chi_atoms[1] [1] = 4; //   CB
  AA[p].chi_atoms[1] [2] = 5; //   CG
  AA[p].chi_atoms[1] [3] = 6; //   CD
  //bk   four atoms that define chi angle  3
  AA[p].chi_atoms[2] [0] = 4; //   CB
  AA[p].chi_atoms[2] [1] = 5; //   CG
  AA[p].chi_atoms[2] [2] = 6; //   CD
  AA[p].chi_atoms[2] [3] = 7; //   OE1


  AA[p].atom[0].icoor[0] = 0.0000; //   N
  AA[p].atom[0].icoor[1] = 7.3000; //   N
  AA[p].atom[0].icoor[2] = 12.6660; //   N
  AA[p].atom[1].icoor[0] = 0.0000; //   CA
  AA[p].atom[1].icoor[1] = 7.2080; //   CA
  AA[p].atom[1].icoor[2] = 14.1210; //   CA
  AA[p].atom[2].icoor[0] = 0.0000; //   C
  AA[p].atom[2].icoor[1] = 8.5910; //   C
  AA[p].atom[2].icoor[2] = 14.7610; //   C
  AA[p].atom[3].icoor[0] = -0.0030; //   O
  AA[p].atom[3].icoor[1] = 9.6070; //   O
  AA[p].atom[3].icoor[2] = 14.0660; //   O
  AA[p].atom[4].icoor[0] = 1.2080; //   CB
  AA[p].atom[4].icoor[1] = 6.4040; //   CB
  AA[p].atom[4].icoor[2] = 14.6070; //   CB
  AA[p].atom[5].icoor[0] = 1.2950; //   CG
  AA[p].atom[5].icoor[1] = 6.2520; //   CG
  AA[p].atom[5].icoor[2] = 16.1190; //   CG
  AA[p].atom[6].icoor[0] = 2.5040; //   CD
  AA[p].atom[6].icoor[1] = 5.4530; //   CD
  AA[p].atom[6].icoor[2] = 16.5190; //   CD
  AA[p].atom[7].icoor[0] = 2.6990; //   OE1
  AA[p].atom[7].icoor[1] = 5.2540; //   OE1
  AA[p].atom[7].icoor[2] = 17.6940; //   OE1
  AA[p].atom[8].icoor[0] = 3.1710; //   OE2
  AA[p].atom[8].icoor[1] = 4.9460; //   OE2
  AA[p].atom[8].icoor[2] = 15.6480; //   OE2
  AA[p].atom[9].icoor[0] = 0.0400; //   H
  AA[p].atom[9].icoor[1] = 8.2170; //   H
  AA[p].atom[9].icoor[2] = 12.2440; //   H
  AA[p].atom[10].icoor[0] = -0.9090; //   HA
  AA[p].atom[10].icoor[1] = 6.7100; //   HA
  AA[p].atom[10].icoor[2] = 14.4600; //   HA
  AA[p].atom[11].icoor[0] = 1.1400; //  2HB
  AA[p].atom[11].icoor[1] = 5.4180; //  2HB
  AA[p].atom[11].icoor[2] = 14.1470; //  2HB
  AA[p].atom[12].icoor[0] = 2.0990; //  3HB
  AA[p].atom[12].icoor[1] = 6.9150; //  3HB
  AA[p].atom[12].icoor[2] = 14.2420; //  3HB
  AA[p].atom[13].icoor[0] = 1.2970; //  2HG
  AA[p].atom[13].icoor[1] = 7.2030; //  2HG
  AA[p].atom[13].icoor[2] = 16.6510; //  2HG
  AA[p].atom[14].icoor[0] = 0.3930; //  3HG
  AA[p].atom[14].icoor[1] = 5.6970; //  3HG
  AA[p].atom[14].icoor[2] = 16.3730; //  3HG
	//NEW ADDED BUT UNKNOWN
  AA[p].atom[15].icoor[0] = 0; //  2HE
  AA[p].atom[15].icoor[1] = 0; //  2HE
  AA[p].atom[15].icoor[2] = 0; //  2HE


// atom number for backbone HN
 AA[p].HNpos=9;
 // atom number for backbone HA
 AA[p].HApos=10;

 // number of polar hydrogens
 AA[p].nH_polar_complete=2;
 // atom numbers for polar H
 AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
 AA[p].Hpos_polar_complete[0]=9;
 AA[p].Hpos_polar_complete[1]=15; //NEW ADDED

 // number of aromatic hydrogens
 AA[p].nH_aromatic_complete=0;

 // number of apolar hydrogens
 AA[p].nH_apolar_complete=5;
 //atom number for apolar hydrogens
 AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
 AA[p].Hpos_apolar_complete[0]=10;
 AA[p].Hpos_apolar_complete[1]=11;
 AA[p].Hpos_apolar_complete[2]=12;
 AA[p].Hpos_apolar_complete[3]=13;
 AA[p].Hpos_apolar_complete[4]=14;

	//A REVISAR
 // number of acceptors
 AA[p].nacceptors=3;
 //acceptor information
 AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
 AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
 (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
 (AA[p].aBase[1])[0]=7;(AA[p].aBase[1])[1]=6;
 (AA[p].aBase[2])[0]=8;(AA[p].aBase[2])[1]=6;
 (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;
 (AA[p].accpt_pos[1])[0]=1;(AA[p].accpt_pos[1])[1]=7;
 (AA[p].accpt_pos[2])[0]=2;(AA[p].accpt_pos[2])[1]=8;

 //Number of Hydrogens connection
 AA[p].nH_hydrogen_connexions=7;
 //atoms hydrogens are connected too
 AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
 AA[p].Hbase[0][0]=9; AA[p].Hbase[0][1]=0;
 AA[p].Hbase[1][0]=10; AA[p].Hbase[1][1]=1;
 AA[p].Hbase[2][0]=11; AA[p].Hbase[2][1]=4;
 AA[p].Hbase[3][0]=12; AA[p].Hbase[3][1]=4;
 AA[p].Hbase[4][0]=13; AA[p].Hbase[4][1]=5;
 AA[p].Hbase[5][0]=14; AA[p].Hbase[5][1]=5;
 AA[p].Hbase[6][0]=15; AA[p].Hbase[6][1]=8; //NEW ADDED

  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=9;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=10;

  AA[p].atom[4].numHydrogens_atm=2;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=11;
  AA[p].atom[4].hydrogens_atm[1]=12;

  AA[p].atom[5].numHydrogens_atm=2;
  AA[p].atom[5].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[5].numHydrogens_atm);
  AA[p].atom[5].hydrogens_atm[0]=13;
  AA[p].atom[5].hydrogens_atm[1]=14;

  AA[p].atom[8].numHydrogens_atm=1;
  AA[p].atom[8].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[5].numHydrogens_atm);
  AA[p].atom[8].hydrogens_atm[0]=15;

}

//FIN DE AMINOACIDO GLUTAMATO NEUTRAL

/**
* Initializes Histidine POSITIVE
*
* @param opt: Rosseta or ICM
*/
void init_HIP(Convention opt)
{
  //AMINOACIDO HISTIDINE POSITIVE
  int k, p =HIP;


  strcpy( AA[p].aa_name3, "HIP" ); // HIE (Protonated at NE2, neutral)
  AA[p].aa_name1 = 'H';            // Note: HID (Protonated at ND1, neutral)
  AA[p].mass = 137.143;            // <---- Note: HIP = HSC (Di-Protonated at ND1 & NE2, +)

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = true;
  AA[p].aa_is_nonpolar = false;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = false;

  AA[p].natoms = 18; //NEW ADDED
  AA[p].nheavyatoms = 10;
  AA[p].nchi = 2;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );
  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;

  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

  AA[p].chi_types = ( int * ) new int[AA[p].nchi];
  AA[p].chi_types[0] = 2;
  AA[p].chi_types[1] = 4;
  AA[p].natoms_EEF1 = 12;


  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " CG " );
  strcpy( AA[p].atom[6].atom_name, " ND1" ); // <--- Protonated here!
  strcpy( AA[p].atom[7].atom_name, " CD2" );
  strcpy( AA[p].atom[8].atom_name, " CE1" );
  strcpy( AA[p].atom[9].atom_name, " NE2" ); // <--- Protonated here!
  strcpy( AA[p].atom[10].atom_name," H  " );
  strcpy( AA[p].atom[11].atom_name," HE2" ); // <--- The proton.
  strcpy( AA[p].atom[12].atom_name," HD1" ); // <--- The proton.
  strcpy( AA[p].atom[13].atom_name," HA " );
  strcpy( AA[p].atom[14].atom_name,"2HB " );
  strcpy( AA[p].atom[15].atom_name,"3HB " );
  strcpy( AA[p].atom[16].atom_name," HE1" );
  strcpy( AA[p].atom[17].atom_name," HD2" );

  strcpy( AA[p].atom[18].atom_name," OXT" );
  strcpy( AA[p].atom[19].atom_name,"2H  " );
  strcpy( AA[p].atom[20].atom_name,"3H  " );


  switch(opt)
  {
  	// NEW ADDED. MODYFIED BY HSP
    case Rosseta:
    case Sybil:
      AA[p].atom[0].charge = -0.47; //    N
      AA[p].atom[1].charge = 0.07; //    CA
      AA[p].atom[2].charge = 0.51; //    C
      AA[p].atom[3].charge = -0.51; //    O
      AA[p].atom[4].charge = -0.05; //    CB
      AA[p].atom[5].charge = 0.19; //    CG
      AA[p].atom[6].charge = -0.51; //    ND1
      AA[p].atom[7].charge =  0.19; //    CD2
      AA[p].atom[8].charge = 0.32; //    CE1  (¿HE1?)
      AA[p].atom[9].charge = -0.51; //    NE2
      AA[p].atom[10].charge = 0.31; //    H
      AA[p].atom[11].charge = 0.44; //    HE2
      AA[p].atom[12].charge = 0.44; //    HD1 //NEW ADDED
      AA[p].atom[13].charge = 0.09; //    HA
      AA[p].atom[14].charge = 0.09; //    2HB
      AA[p].atom[15].charge = 0.09; //    3HB
      AA[p].atom[16].charge = 0.18; //    HE1
      AA[p].atom[17].charge = 0.13; //    HD2

      AA[p].atom[18].charge = -0.67; //   0XT
      AA[p].atom[19].charge = 0.33; //  1H
      AA[p].atom[20].charge = 0.33; //  2H
      break;

    //A REVISAR
    case ICM:
      AA[p].atom[0].charge = -0.356; //    N
      AA[p].atom[1].charge = 0.064; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = -0.05; //    CB
      AA[p].atom[5].charge = 0.15; //    CG
      AA[p].atom[6].charge = -0.2; //    ND1
      AA[p].atom[7].charge = 0.1; //    CD2
      AA[p].atom[8].charge = 0.275; //    CE1
      AA[p].atom[9].charge = -0.2; //    NE2
      AA[p].atom[10].charge = 0.176; //    H
      AA[p].atom[11].charge = 0.265; //    HE2
      AA[p].atom[12].charge = 0.135; //    HD1 // NEW ADDED
      AA[p].atom[13].charge = 0.02; //    HA
      AA[p].atom[14].charge = 0.065; //    2HB
      AA[p].atom[15].charge = 0.065; //    3HB
      AA[p].atom[16].charge = 0.15; //    HE1
      AA[p].atom[17].charge = 0.135; //    HD2

      AA[p].atom[18].charge = -0.67; //   0XT
      AA[p].atom[18].charge = 0.33; //  1H
      AA[p].atom[20].charge = 0.33; //  2H
      break;
   case EEF1:
      AA[p].atom[0].charge = -0.35; //    N
      AA[p].atom[1].charge = 0.10; //    CA
      AA[p].atom[2].charge = 0.55; //    C
      AA[p].atom[3].charge = -0.55; //    O
      AA[p].atom[4].charge = 0.10; //    CB
      AA[p].atom[5].charge = 0.15; //    CG
      AA[p].atom[6].charge = -0.30; //    ND1
      AA[p].atom[7].charge = 0.20; //    CD2
      AA[p].atom[8].charge = 0.45; //    CE1
      AA[p].atom[9].charge = -0.30; //    NE2
      AA[p].atom[10].charge = 0.25; //    H
      AA[p].atom[11].charge = 0.35; //    HE2
      AA[p].atom[12].charge = 0.35; //    HD1
      AA[p].atom[13].charge = 0.00; //    HA
      AA[p].atom[14].charge = 0.00; //    2HB
      AA[p].atom[15].charge = 0.00; //    3HB
      AA[p].atom[16].charge = 0.00; //    HE1
      AA[p].atom[17].charge = 0.00; //    HD2
      break;
  }
  /*
   1.0 HSE-N:2H    0.33000 HC
   1.0 HSE-N:O    -0.51000 O
   1.0 HSE-N:HA    0.10000 HB
   1.0 HSE-N:3H    0.33000 HC
   1.0 HSE-N:CB   -0.08000 CT2
   1.0 HSE-N:HB1   0.09000 HA
   1.0 HSE-N:HB2   0.09000 HA
   1.0 HSE-N:ND1  -0.70000 NR2
   1.0 HSE-N:CG    0.22000 CPH1
   1.0 HSE-N:CE1   0.25000 CPH2
   1.0 HSE-N:HE1   0.13000 HR1
   1.0 HSE-N:NE2  -0.36000 NR1
   1.0 HSE-N:HE2   0.32000 H
   1.0 HSE-N:CD2  -0.05000 CPH1
   1.0 HSE-N:N    -0.30000 NH3
   1.0 HSE-N:HD2   0.09000 HR3
   1.0 HSE-N:1H    0.33000 HC
   1.0 HSE-N:C     0.51000 C
   1.0 HSE-N:CA    0.21000 CT1

   1.0 HSE-C:HA    0.09000 HB
   1.0 HSE-C:CB   -0.08000 CT2
   1.0 HSE-C:HB1   0.09000 HA
   1.0 HSE-C:HB2   0.09000 HA
   1.0 HSE-C:ND1  -0.70000 NR2
   1.0 HSE-C:CG    0.22000 CPH1
   1.0 HSE-C:CE1   0.25000 CPH2
   1.0 HSE-C:HE1   0.13000 HR1
   1.0 HSE-C:O    -0.67000 OC
   1.0 HSE-C:NE2  -0.36000 NR1
   1.0 HSE-C:OXT  -0.67000 OC
   1.0 HSE-C:HE2   0.32000 H
   1.0 HSE-C:CD2  -0.05000 CPH1
   1.0 HSE-C:N    -0.47000 NH1
   1.0 HSE-C:HD2   0.09000 HR3
   1.0 HSE-C:HN    0.31000 H
   1.0 HSE-C:C     0.34000 CC
   1.0 HSE-C:CA    0.07000 CT1
*/

  switch (opt) {
  	case EEF1:
	  AA[p].atom[0].fullatom_type =  12; // NH1    N
	  AA[p].atom[1].fullatom_type =   5; // CH1E   CA
	  AA[p].atom[2].fullatom_type =   4; // C      C
	  AA[p].atom[3].fullatom_type =  16; // O      O
	  AA[p].atom[4].fullatom_type =   6; // CH2E   CB
	  AA[p].atom[5].fullatom_type =  24; // CR     CG
	  AA[p].atom[6].fullatom_type =  12; // NH1    ND1
	  AA[p].atom[7].fullatom_type =   8; // CR1E   CD2
	  AA[p].atom[8].fullatom_type =   8; // CR1E   CE1
	  AA[p].atom[9].fullatom_type =  12; // NH1    NE2
	  AA[p].atom[10].fullatom_type =  0; // H      H
	  AA[p].atom[11].fullatom_type =  0; // H      HE2
	  AA[p].atom[12].fullatom_type =  0; // H      HD1 NEW ADDED BUT UNKNOWN
	  AA[p].atom[13].fullatom_type = -1; // Hapo   HA
	  AA[p].atom[14].fullatom_type = -1; // Hapo   2HB
	  AA[p].atom[15].fullatom_type = -1; // Hapo   3HB
	  AA[p].atom[16].fullatom_type = -1; // Hapo   HE1
	  AA[p].atom[17].fullatom_type = -1; // Hapo   HD2

	  AA[p].atom[18].fullatom_type = 17; //OC OXT
	  AA[p].atom[19].fullatom_type =  1; //HC Hpol
	  AA[p].atom[20].fullatom_type =  1; //HC Hpol
	  break;
  	case Sybil:


	  AA[p].atom[0].fullatom_type = 40; // Nbb    N
	  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 10; // CObb   C
	  AA[p].atom[3].fullatom_type = 56; // OCbb   O
	  AA[p].atom[4].fullatom_type = 8; // CH2    CB
	  AA[p].atom[5].fullatom_type = 11; // aroC   CG
	  AA[p].atom[6].fullatom_type = 44; // Nhis   ND1
	  AA[p].atom[7].fullatom_type = 11; // aroC   CD2
	  AA[p].atom[8].fullatom_type = 11; // aroC   CE1
	  AA[p].atom[9].fullatom_type = 42; // Ntrp   NE2
	  AA[p].atom[10].fullatom_type = 66; // HNbb   H
	  AA[p].atom[11].fullatom_type = 66; // Hpol   HE2
	  AA[p].atom[12].fullatom_type = 66; // Hapo   HD1 NEW ADDED BUT UNKNOWN
	  AA[p].atom[13].fullatom_type = 66; // Hapo   HA
	  AA[p].atom[14].fullatom_type = 66; // Hapo  2HB
	  AA[p].atom[15].fullatom_type = 66; // Hapo  3HB
	  AA[p].atom[16].fullatom_type = 66; // Hapo   HE1
	  AA[p].atom[17].fullatom_type = 66; // Hapo   HD2

	  AA[p].atom[18].fullatom_type = 57; // OXT
	  AA[p].atom[19].fullatom_type = 66; // Hpol
	  AA[p].atom[20].fullatom_type = 66; // Hpol
      break;
      default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 6; // aroC   CG
	  AA[p].atom[6].fullatom_type = 8; // Nhis   ND1
	  AA[p].atom[7].fullatom_type = 6; // aroC   CD2
	  AA[p].atom[8].fullatom_type = 6; // aroC   CE1
	  AA[p].atom[9].fullatom_type = 7; // Ntrp   NE2
	  AA[p].atom[10].fullatom_type = 25; // HNbb   H
	  AA[p].atom[11].fullatom_type = 22; // Hpol   HE2
	  AA[p].atom[12].fullatom_type = 23; // Hapo   HD1 NEW ADDED BUT UNKNOWN
	  AA[p].atom[13].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[14].fullatom_type = 23; // Hapo  2HB
	  AA[p].atom[15].fullatom_type = 23; // Hapo  3HB
	  AA[p].atom[16].fullatom_type = 23; // Hapo   HE1
	  AA[p].atom[17].fullatom_type = 23; // Hapo   HD2

	  AA[p].atom[18].fullatom_type = 20; // OXT
	  AA[p].atom[19].fullatom_type = 22; // Hpol
	  AA[p].atom[20].fullatom_type = 22; // Hpol
      break;
  }

  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 10; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 13; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--CG
  AA[p].atom[4].bonded_neighbor[2] = 13; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 14; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 3; // CG
  AA[p].atom[5].bonded_neighbor[0] = 4; // CG--CB
  AA[p].atom[5].bonded_neighbor[1] = 6; // CG--ND1
  AA[p].atom[5].bonded_neighbor[2] = 7; // CG--CD2
  AA[p].atom[6].nbonded_neighbors = 3; // ND1
  AA[p].atom[6].bonded_neighbor[0] = 5; // ND1--CG
  AA[p].atom[6].bonded_neighbor[1] = 8; // ND1--CE1
  AA[p].atom[6].bonded_neighbor[1] = 12; // ND1--HD1 NEW ADDED
  AA[p].atom[7].nbonded_neighbors = 3; // CD2
  AA[p].atom[7].bonded_neighbor[0] = 5; // CD2--CG
  AA[p].atom[7].bonded_neighbor[1] = 9; // CD2--NE2
  AA[p].atom[7].bonded_neighbor[2] = 17; // CD2--HD2
  AA[p].atom[8].nbonded_neighbors = 3; // CE1
  AA[p].atom[8].bonded_neighbor[0] = 6; // CE1--ND1
  AA[p].atom[8].bonded_neighbor[1] = 9; // CE1--NE2
  AA[p].atom[8].bonded_neighbor[2] = 16; // CE1--HE1
  AA[p].atom[9].nbonded_neighbors = 3; // NE2
  AA[p].atom[9].bonded_neighbor[0] = 7; // NE2--CD2
  AA[p].atom[9].bonded_neighbor[1] = 8; // NE2--CE1
  AA[p].atom[9].bonded_neighbor[2] = 11; // NE2--HE2
  AA[p].atom[10].nbonded_neighbors = 1; // H
  AA[p].atom[10].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[11].nbonded_neighbors = 1; // HE2
  AA[p].atom[11].bonded_neighbor[0] = 9; // HE2--NE2
  AA[p].atom[12].nbonded_neighbors = 1; // HD1
  AA[p].atom[12].bonded_neighbor[0] = 6; // HD1--ND1 NEW ADDED
  AA[p].atom[13].nbonded_neighbors = 1; // HA
  AA[p].atom[13].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[14].nbonded_neighbors = 1; //2HB
  AA[p].atom[14].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[15].nbonded_neighbors = 1; //3HB
  AA[p].atom[15].bonded_neighbor[0] = 4; //3HB--CB
  AA[p].atom[16].nbonded_neighbors = 1; // HE1
  AA[p].atom[16].bonded_neighbor[0] = 8; // HE1--CE1
  AA[p].atom[17].nbonded_neighbors = 1; // HD2
  AA[p].atom[17].bonded_neighbor[0] = 7; // HD2--CD2


  AA[p].atom[18].nbonded_neighbors = 1; // OXT
  AA[p].atom[18].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[19].nbonded_neighbors = 1; // 1H
  AA[p].atom[19].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[20].nbonded_neighbors = 1; // 2H
  AA[p].atom[20].bonded_neighbor[0] = 0; // H-N


  for ( k = 0; k < 11; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 11; k < 17; k++ )
    AA[p].atom[k].hastemplate = true;

	//A REVISAR TEMPLATES Y CHIS
  //bk   template for building  HE2
  AA[p].atom[11].ta[0] = 9; //   NE2
  AA[p].atom[11].ta[1] = 7; //   CD2
  AA[p].atom[11].ta[2] = 8; //   CE1

  //bk   template for building  HA
  AA[p].atom[13].ta[0] = 1; //   CA
  AA[p].atom[13].ta[1] = 0; //   N
  AA[p].atom[13].ta[2] = 2; //   C

  //bk   template for building  2HB
  AA[p].atom[14].ta[0] = 4; //   CB
  AA[p].atom[14].ta[1] = 1; //   CA
  AA[p].atom[14].ta[2] = 5; //   CG

  //bk   template for building  3HB
  AA[p].atom[15].ta[0] = 4; //   CB
  AA[p].atom[15].ta[1] = 1; //   CA
  AA[p].atom[15].ta[2] = 5; //   CG

  //bk   template for building  HE1
  AA[p].atom[16].ta[0] = 8; //   CE1
  AA[p].atom[16].ta[1] = 6; //   ND1
  AA[p].atom[16].ta[2] = 9; //  NE2

  //bk   template for building  HD2
  AA[p].atom[17].ta[0] = 7; //   CD2
  AA[p].atom[17].ta[1] = 5; //   CG
  AA[p].atom[17].ta[2] = 9; //   NE2




  //bk   chi angles required to build atoms HIS
  //bk   chi angles needed for building  CG:HE1
  AA[p].chi_required[0] [5] = true;

  //bk   chi angles needed for building  ND1
  AA[p].chi_required[0] [6] = true;
  AA[p].chi_required[1] [6] = true;

  //bk   chi angles needed for building  CD2
  AA[p].chi_required[0] [7] = true;
  AA[p].chi_required[1] [7] = true;

  //bk   chi angles needed for building  CE1
  AA[p].chi_required[0] [8] = true;
  AA[p].chi_required[1] [8] = true;

  //bk   chi angles needed for building  NE2
  AA[p].chi_required[0] [9] = true;
  AA[p].chi_required[1] [9] = true;

  //bk   chi angles needed for building  HE2
  AA[p].chi_required[0] [11] = true;
  AA[p].chi_required[1] [11] = true;

  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [14] = true;

  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [15] = true;

  //bk   chi angles needed for building  HE1
  AA[p].chi_required[0] [16] = true;
  AA[p].chi_required[1] [16] = true;

  //bk   chi angles needed for building  HD2
  AA[p].chi_required[0] [17] = true;
  AA[p].chi_required[1] [17] = true;

  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   CG
  //bk   four atoms that define chi angle  2
  AA[p].chi_atoms[1] [0] = 1; //   CA
  AA[p].chi_atoms[1] [1] = 4; //   CB
  AA[p].chi_atoms[1] [2] = 5; //   CG:HE1
  AA[p].chi_atoms[1] [3] = 6; //   ND1



  AA[p].atom[0].icoor[0] = 0.0000; //   N
  AA[p].atom[0].icoor[1] = 12.0140; //   N
  AA[p].atom[0].icoor[2] = 22.5220; //   N
  AA[p].atom[1].icoor[0] = 0.0000; //   CA
  AA[p].atom[1].icoor[1] = 13.2300; //   CA
  AA[p].atom[1].icoor[2] = 23.3160; //   CA
  AA[p].atom[2].icoor[0] = 0.0000; //   C
  AA[p].atom[2].icoor[1] = 12.6360; //   C
  AA[p].atom[2].icoor[2] = 24.7310; //   C
  AA[p].atom[3].icoor[0] = 0.1540; //   O
  AA[p].atom[3].icoor[1] = 11.4260; //   O:HE1
  AA[p].atom[3].icoor[2] = 24.9130; //   O
  AA[p].atom[4].icoor[0] = -1.2080; //   CB
  AA[p].atom[4].icoor[1] = 14.1050; //   CB
  AA[p].atom[4].icoor[2] = 22.9660; //   CB
  AA[p].atom[5].icoor[0] = -1.2740; //   CG
  AA[p].atom[5].icoor[1] = 15.3800; //   CG
  AA[p].atom[5].icoor[2] = 23.7480; //   CG
  AA[p].atom[6].icoor[0] = -2.2840; //   ND1
  AA[p].atom[6].icoor[1] = 16.3040; //   ND1
  AA[p].atom[6].icoor[2] = 23.5800; //   ND1
  AA[p].atom[7].icoor[0] = -0.4560; //   CD2
  AA[p].atom[7].icoor[1] = 15.8850; //   CD2
  AA[p].atom[7].icoor[2] = 24.7010; //   CD2
  AA[p].atom[8].icoor[0] = -2.0840; //   CE1
  AA[p].atom[8].icoor[1] = 17.3230; //   CE1:HE1
  AA[p].atom[8].icoor[2] = 24.3980; //   CE1
  AA[p].atom[9].icoor[0] = -0.9820; //   NE2
  AA[p].atom[9].icoor[1] = 17.0930; //   NE2
  AA[p].atom[9].icoor[2] = 25.0880; //   NE2
  AA[p].atom[10].icoor[0] = -0.0410; //   H
  AA[p].atom[10].icoor[1] = 11.1270; //   H
  AA[p].atom[10].icoor[2] = 23.0030; //   H
  AA[p].atom[11].icoor[0] = -0.5830; //   HE2
  AA[p].atom[11].icoor[1] = 17.7010; //   HE2
  AA[p].atom[11].icoor[2] = 25.7890; //   HE2
  AA[p].atom[12].icoor[0] = 0.0; //   HD1 NEW ADDED BUT UNKNOWN
  AA[p].atom[12].icoor[1] = 0.0; //   HD1
  AA[p].atom[12].icoor[2] = 0.0; //   HD1
  AA[p].atom[13].icoor[0] = 0.9100; //   HA
  AA[p].atom[13].icoor[1] = 13.7960; //   HA
  AA[p].atom[13].icoor[2] = 23.1190; //   HA
  AA[p].atom[14].icoor[0] = -1.1760; //  2HB
  AA[p].atom[14].icoor[1] = 14.3870; //  2HB:HE1
  AA[p].atom[14].icoor[2] = 21.9130; //  2HB
  AA[p].atom[15].icoor[0] = -2.1330; //  3HB
  AA[p].atom[15].icoor[1] = 13.5660; //  3HB
  AA[p].atom[15].icoor[2] = 23.1690; //  3HB
  AA[p].atom[16].icoor[0] = -2.7780; //   HE1
  AA[p].atom[16].icoor[1] = 18.1630; //   HE1
  AA[p].atom[16].icoor[2] = 24.4180; //   HE1
  AA[p].atom[17].icoor[0] = 0.4620; //   HD2
  AA[p].atom[17].icoor[1] = 15.5170; //   HD2
  AA[p].atom[17].icoor[2] = 25.1600; //   HD2

  // atom number for backbone HN
  AA[p].HNpos=10;
  // atom number for backbone HA
  AA[p].HApos=13;

  // number of polar hydrogens
  AA[p].nH_polar_complete=3;
  // atom numbers for polar H
  AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
  AA[p].Hpos_polar_complete[0]=10;
  AA[p].Hpos_polar_complete[1]=11;
  AA[p].Hpos_polar_complete[2]=12; //NEW ADDED

  // number of aromatic hydrogens
  AA[p].nH_aromatic_complete=0;

  // number of apolar hydrogens
  AA[p].nH_apolar_complete=5;
  //atom number for apolar hydrogens
  AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
  AA[p].Hpos_apolar_complete[0]=14;
  AA[p].Hpos_apolar_complete[1]=15;
  AA[p].Hpos_apolar_complete[2]=16;
  AA[p].Hpos_apolar_complete[3]=17;
  AA[p].Hpos_apolar_complete[4]=12;

	//A REVISAR
  // number of acceptors
  AA[p].nacceptors=2;
  //acceptor information
  AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
  (AA[p].aBase[1])[0]=6;(AA[p].aBase[1])[1]=5;
  (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;
  (AA[p].accpt_pos[1])[0]=1;(AA[p].accpt_pos[1])[1]=6;
  //OJO. HAY UNA RELACION ABASE2(6,8)

  //Number of Hydrogens connection
  AA[p].nH_hydrogen_connexions=8;
  //atoms hydrogens are connected too
  AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
  AA[p].Hbase[0][0]=10; AA[p].Hbase[0][1]=0;
  AA[p].Hbase[1][0]=11; AA[p].Hbase[1][1]=9;
  AA[p].Hbase[2][0]=14; AA[p].Hbase[2][1]=1;
  AA[p].Hbase[3][0]=15; AA[p].Hbase[3][1]=4;
  AA[p].Hbase[4][0]=16; AA[p].Hbase[4][1]=4;
  AA[p].Hbase[5][0]=17; AA[p].Hbase[5][1]=8;
  AA[p].Hbase[6][0]=12; AA[p].Hbase[6][1]=7;
  AA[p].Hbase[7][0]=13; AA[p].Hbase[7][1]=6; //NEW ADDED

  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=10;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=14;

  AA[p].atom[4].numHydrogens_atm=2;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=15;
  AA[p].atom[4].hydrogens_atm[1]=16;

  AA[p].atom[7].numHydrogens_atm=1;
  AA[p].atom[7].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[7].numHydrogens_atm);
  AA[p].atom[7].hydrogens_atm[0]=12;

  AA[p].atom[8].numHydrogens_atm=1;
  AA[p].atom[8].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[8].numHydrogens_atm);
  AA[p].atom[8].hydrogens_atm[0]=17;

  AA[p].atom[9].numHydrogens_atm=1;
  AA[p].atom[9].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[9].numHydrogens_atm);
  AA[p].atom[9].hydrogens_atm[0]=11;

	//NEW ADDED
  AA[p].atom[6].numHydrogens_atm=1;
  AA[p].atom[6].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[9].numHydrogens_atm);
  AA[p].atom[6].hydrogens_atm[0]=13;


  //FIN DE AMINOACIDO HISTIDINE
}

/**
* Initializes Histidine NEUTRAL PROTON HD1
*
* @param opt: Rosseta or ICM
*/
void init_HID(Convention opt)
{
  //AMINOACIDO HISTIDINE
  int k, p =HID;


  strcpy( AA[p].aa_name3, "HID" ); // HIE (Protonated at NE2, neutral)
  AA[p].aa_name1 = 'H';            // <---- Note: HID = HIS (Protonated at ND1, neutral)
  AA[p].mass = 137.143;            // Note: HIP (Di-Protonated at ND1 & NE2, +)

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = true;
  AA[p].aa_is_nonpolar = false;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = false;

  AA[p].natoms = 17;
  AA[p].nheavyatoms = 10;
  AA[p].nchi = 2;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );
  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;

  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

  AA[p].chi_types = ( int * ) new int[AA[p].nchi];
  AA[p].chi_types[0] = 2;
  AA[p].chi_types[1] = 4;
  AA[p].natoms_EEF1 = 11;

  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " CG " );
  strcpy( AA[p].atom[6].atom_name, " ND1" ); // <--- Protonated here!
  strcpy( AA[p].atom[7].atom_name, " CD2" );
  strcpy( AA[p].atom[8].atom_name, " CE1" );
  strcpy( AA[p].atom[9].atom_name, " NE2" );
  strcpy( AA[p].atom[10].atom_name," H  " );
  strcpy( AA[p].atom[11].atom_name," HD1" ); //proton.
  strcpy( AA[p].atom[12].atom_name," HA " );
  strcpy( AA[p].atom[13].atom_name,"2HB " );
  strcpy( AA[p].atom[14].atom_name,"3HB " );
  strcpy( AA[p].atom[15].atom_name," HE1" );
  strcpy( AA[p].atom[16].atom_name," HD2" );

  strcpy( AA[p].atom[17].atom_name," OXT" );
  strcpy( AA[p].atom[18].atom_name,"2H  " );
  strcpy( AA[p].atom[19].atom_name,"3H  " );


  switch(opt)
  {
    //MODIFIED BY HSD
    case Rosseta:
    case Sybil:
      AA[p].atom[0].charge = -0.47; //    N
      AA[p].atom[1].charge = 0.07; //    CA
      AA[p].atom[2].charge = 0.51; //    C
      AA[p].atom[3].charge = -0.51; //    O
      AA[p].atom[4].charge = -0.09; //    CB
      AA[p].atom[5].charge = -0.05; //    CG
      AA[p].atom[6].charge = -0.36; //    ND1
      AA[p].atom[7].charge =  0.22; //    CD2
      AA[p].atom[8].charge = 0.25; //    CE1  (¿HE1?)
      AA[p].atom[9].charge = -0.70; //    NE2
      AA[p].atom[10].charge = 0.31; //    H
      AA[p].atom[11].charge = 0.32; //    HD1 NEW ADDED
      AA[p].atom[12].charge = 0.09; //    HA
      AA[p].atom[13].charge = 0.09; //    2HB
      AA[p].atom[14].charge = 0.09; //    3HB
      AA[p].atom[15].charge = 0.13; //    HE1
      AA[p].atom[16].charge = 0.10; //    HD2

      AA[p].atom[17].charge = -0.67; //   0XT
      AA[p].atom[18].charge = 0.33; //  1H
      AA[p].atom[19].charge = 0.33; //  2H
      break;
		//A REVISAR
    case ICM:
      AA[p].atom[0].charge = -0.356; //    N
      AA[p].atom[1].charge = 0.064; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = -0.05; //    CB
      AA[p].atom[5].charge = 0.15; //    CG
      AA[p].atom[6].charge = -0.2; //    ND1
      AA[p].atom[7].charge = 0.1; //    CD2
      AA[p].atom[8].charge = 0.275; //    CE1
      AA[p].atom[9].charge = -0.2; //    NE2
      AA[p].atom[10].charge = 0.176; //    H
      AA[p].atom[11].charge = 0.265; //    HD1 NEW ADDED
      AA[p].atom[12].charge = 0.02; //    HA
      AA[p].atom[13].charge = 0.065; //    2HB
      AA[p].atom[14].charge = 0.065; //    3HB
      AA[p].atom[15].charge = 0.15; //    HE1
      AA[p].atom[16].charge = 0.135; //    HD2

      AA[p].atom[17].charge = -0.67; //   0XT
      AA[p].atom[18].charge = 0.33; //  1H
      AA[p].atom[19].charge = 0.33; //  2H
      break;
    case EEF1:
      AA[p].atom[0].charge = -0.35; //    N
      AA[p].atom[1].charge = 0.10; //    CA
      AA[p].atom[2].charge = 0.55; //    C
      AA[p].atom[3].charge = -0.55; //    O
      AA[p].atom[4].charge = 0.00; //    CB
      AA[p].atom[5].charge = 0.10; //    CG
      AA[p].atom[6].charge = -0.40; //    ND1
      AA[p].atom[7].charge = 0.10; //    CD2
      AA[p].atom[8].charge = 0.30; //    CE1
      AA[p].atom[9].charge = -0.40; //    NE2
      AA[p].atom[10].charge = 0.25; //    H
      AA[p].atom[11].charge = 0.30; //    HD1
      AA[p].atom[12].charge = 0.00; //    HA
      AA[p].atom[13].charge = 0.00; //    2HB
      AA[p].atom[14].charge = 0.00; //    3HB
      AA[p].atom[15].charge = 0.00; //    HE1
      AA[p].atom[16].charge = 0.00; //    HD2
      break;
  }
  /*
   1.0 HSE-N:2H    0.33000 HC
   1.0 HSE-N:O    -0.51000 O
   1.0 HSE-N:HA    0.10000 HB
   1.0 HSE-N:3H    0.33000 HC
   1.0 HSE-N:CB   -0.08000 CT2
   1.0 HSE-N:HB1   0.09000 HA
   1.0 HSE-N:HB2   0.09000 HA
   1.0 HSE-N:ND1  -0.70000 NR2
   1.0 HSE-N:CG    0.22000 CPH1
   1.0 HSE-N:CE1   0.25000 CPH2
   1.0 HSE-N:HE1   0.13000 HR1
   1.0 HSE-N:NE2  -0.36000 NR1
   1.0 HSE-N:HE2   0.32000 H
   1.0 HSE-N:CD2  -0.05000 CPH1
   1.0 HSE-N:N    -0.30000 NH3
   1.0 HSE-N:HD2   0.09000 HR3
   1.0 HSE-N:1H    0.33000 HC
   1.0 HSE-N:C     0.51000 C
   1.0 HSE-N:CA    0.21000 CT1

   1.0 HSE-C:HA    0.09000 HB
   1.0 HSE-C:CB   -0.08000 CT2
   1.0 HSE-C:HB1   0.09000 HA
   1.0 HSE-C:HB2   0.09000 HA
   1.0 HSE-C:ND1  -0.70000 NR2
   1.0 HSE-C:CG    0.22000 CPH1
   1.0 HSE-C:CE1   0.25000 CPH2
   1.0 HSE-C:HE1   0.13000 HR1
   1.0 HSE-C:O    -0.67000 OC
   1.0 HSE-C:NE2  -0.36000 NR1
   1.0 HSE-C:OXT  -0.67000 OC
   1.0 HSE-C:HE2   0.32000 H
   1.0 HSE-C:CD2  -0.05000 CPH1
   1.0 HSE-C:N    -0.47000 NH1
   1.0 HSE-C:HD2   0.09000 HR3
   1.0 HSE-C:HN    0.31000 H
   1.0 HSE-C:C     0.34000 CC
   1.0 HSE-C:CA    0.07000 CT1
*/

  switch (opt) {
  	case EEF1:
	  AA[p].atom[0].fullatom_type =  12; // NH1    N
	  AA[p].atom[1].fullatom_type =   5; // CH1E   CA
	  AA[p].atom[2].fullatom_type =   4; // C      C
	  AA[p].atom[3].fullatom_type =  16; // O      O
	  AA[p].atom[4].fullatom_type =   6; // CH2E   CB
	  AA[p].atom[5].fullatom_type =  24; // CR     CG
	  AA[p].atom[6].fullatom_type =  12; // NH1    ND1
	  AA[p].atom[7].fullatom_type =   8; // CR1E   CD2
	  AA[p].atom[8].fullatom_type =   8; // CR1E   CE1
	  AA[p].atom[9].fullatom_type =  10; // NR     NE2
	  AA[p].atom[10].fullatom_type =  0; // H      H
	  AA[p].atom[11].fullatom_type =  0; // H      HD1 NEW ADDED NUT UNKNOWN
	  AA[p].atom[12].fullatom_type = -1; // Hapo   HA
	  AA[p].atom[13].fullatom_type = -1; // Hapo   2HB
	  AA[p].atom[14].fullatom_type = -1; // Hapo   3HB
	  AA[p].atom[15].fullatom_type = -1; // Hapo   HE1
	  AA[p].atom[16].fullatom_type = -1; // Hapo   HD2

	  AA[p].atom[17].fullatom_type = 17; // OC     OXT
	  AA[p].atom[18].fullatom_type =  1; // HC     2H
	  AA[p].atom[19].fullatom_type =  1; // HC     3H
	  break;
  	case Sybil:



  	  AA[p].atom[0].fullatom_type = 40; // Nbb    N
  	  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
  	  AA[p].atom[2].fullatom_type = 10; // CObb   C
  	  AA[p].atom[3].fullatom_type = 56; // OCbb   O
  	  AA[p].atom[4].fullatom_type = 8; // CH2    CB
  	  AA[p].atom[5].fullatom_type = 11; // aroC   CG
  	  AA[p].atom[6].fullatom_type = 44; // Nhis   ND1
  	  AA[p].atom[7].fullatom_type = 11; // aroC   CD2
  	  AA[p].atom[8].fullatom_type = 11; // aroC   CE1
  	  AA[p].atom[9].fullatom_type = 42; // Ntrp   NE2
	  AA[p].atom[10].fullatom_type = 66; // HNbb   H
	  AA[p].atom[11].fullatom_type = 66; // Hapo   HD1 NEW ADDED BUT UNKNOWN
	  AA[p].atom[12].fullatom_type = 66; // Hapo   HA
	  AA[p].atom[13].fullatom_type = 66; // Hapo  2HB
	  AA[p].atom[14].fullatom_type = 66; // Hapo  3HB
	  AA[p].atom[15].fullatom_type = 66; // Hapo   HE1
	  AA[p].atom[16].fullatom_type = 66; // Hapo   HD2

	  AA[p].atom[17].fullatom_type = 57; // OXT
	  AA[p].atom[18].fullatom_type = 66; // Hpol
	  AA[p].atom[19].fullatom_type = 66; // Hpol
	  break;
	default:
      AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 6; // aroC   CG
	  AA[p].atom[6].fullatom_type = 8; // Nhis   ND1
	  AA[p].atom[7].fullatom_type = 6; // aroC   CD2
	  AA[p].atom[8].fullatom_type = 6; // aroC   CE1
	  AA[p].atom[9].fullatom_type = 7; // Ntrp   NE2
	  AA[p].atom[10].fullatom_type = 25; // HNbb   H
	  AA[p].atom[11].fullatom_type = 23; // Hapo   HD1 NEW ADDED BUT UNKNOWN
	  AA[p].atom[12].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[13].fullatom_type = 23; // Hapo  2HB
	  AA[p].atom[14].fullatom_type = 23; // Hapo  3HB
	  AA[p].atom[15].fullatom_type = 23; // Hapo   HE1
	  AA[p].atom[16].fullatom_type = 23; // Hapo   HD2

	  AA[p].atom[17].fullatom_type = 20; // OXT
	  AA[p].atom[18].fullatom_type = 22; // Hpol
	  AA[p].atom[19].fullatom_type = 22; // Hpol
	  break;
  }

  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 10; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 12; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--CG
  AA[p].atom[4].bonded_neighbor[2] = 13; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 14; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 3; // CG
  AA[p].atom[5].bonded_neighbor[0] = 4; // CG--CB
  AA[p].atom[5].bonded_neighbor[1] = 6; // CG--ND1
  AA[p].atom[5].bonded_neighbor[2] = 7; // CG--CD2
  AA[p].atom[6].nbonded_neighbors = 3; // ND1
  AA[p].atom[6].bonded_neighbor[0] = 5; // ND1--CG
  AA[p].atom[6].bonded_neighbor[1] = 8; // ND1--CE1
  AA[p].atom[6].bonded_neighbor[2] = 11; // ND1--HD1 //NEW ADDED
  AA[p].atom[7].nbonded_neighbors = 3; // CD2
  AA[p].atom[7].bonded_neighbor[0] = 5; // CD2--CG
  AA[p].atom[7].bonded_neighbor[1] = 9; // CD2--NE2
  AA[p].atom[7].bonded_neighbor[2] = 16; // CD2--HD2
  AA[p].atom[8].nbonded_neighbors = 3; // CE1
  AA[p].atom[8].bonded_neighbor[0] = 6; // CE1--ND1
  AA[p].atom[8].bonded_neighbor[1] = 9; // CE1--NE2
  AA[p].atom[8].bonded_neighbor[2] = 15; // CE1--HE1
  AA[p].atom[9].nbonded_neighbors = 2; // NE2
  AA[p].atom[9].bonded_neighbor[0] = 7; // NE2--CD2
  AA[p].atom[9].bonded_neighbor[1] = 8; // NE2--CE1
  AA[p].atom[10].nbonded_neighbors = 1; // H
  AA[p].atom[10].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[11].nbonded_neighbors = 1; // HE2
  AA[p].atom[11].bonded_neighbor[0] = 6; // HD1--ND1 /NEW ADDED
  AA[p].atom[12].nbonded_neighbors = 1; // HA
  AA[p].atom[12].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[13].nbonded_neighbors = 1; //2HB
  AA[p].atom[13].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[14].nbonded_neighbors = 1; //3HB
  AA[p].atom[14].bonded_neighbor[0] = 4; //3HB--CB
  AA[p].atom[15].nbonded_neighbors = 1; // HE1
  AA[p].atom[15].bonded_neighbor[0] = 8; // HE1--CE1
  AA[p].atom[16].nbonded_neighbors = 1; // HD2
  AA[p].atom[16].bonded_neighbor[0] = 7; // HD2--CD2

  AA[p].atom[17].nbonded_neighbors = 1; // OXT
  AA[p].atom[17].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[18].nbonded_neighbors = 1; // 1H
  AA[p].atom[18].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[19].nbonded_neighbors = 1; // 2H
  AA[p].atom[19].bonded_neighbor[0] = 0; // H-N


  //A REVISAR TEMPLATES Y CHIS
  for ( k = 0; k < 11; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 11; k < 17; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template for building  HE2
  AA[p].atom[11].ta[0] = 9; //   NE2
  AA[p].atom[11].ta[1] = 7; //   CD2
  AA[p].atom[11].ta[2] = 8; //   CE1

  //bk   template for building  HA
  AA[p].atom[12].ta[0] = 1; //   CA
  AA[p].atom[12].ta[1] = 0; //   N
  AA[p].atom[12].ta[2] = 2; //   C

  //bk   template for building  2HB
  AA[p].atom[13].ta[0] = 4; //   CB
  AA[p].atom[13].ta[1] = 1; //   CA
  AA[p].atom[13].ta[2] = 5; //   CG

  //bk   template for building  3HB
  AA[p].atom[14].ta[0] = 4; //   CB
  AA[p].atom[14].ta[1] = 1; //   CA
  AA[p].atom[14].ta[2] = 5; //   CG

  //bk   template for building  HE1
  AA[p].atom[15].ta[0] = 8; //   CE1
  AA[p].atom[15].ta[1] = 6; //   ND1
  AA[p].atom[15].ta[2] = 9; //  NE2

  //bk   template for building  HD2
  AA[p].atom[16].ta[0] = 7; //   CD2
  AA[p].atom[16].ta[1] = 5; //   CG
  AA[p].atom[16].ta[2] = 9; //   NE2




  //bk   chi angles required to build atoms HIS
  //bk   chi angles needed for building  CG:HE1
  AA[p].chi_required[0] [5] = true;

  //bk   chi angles needed for building  ND1
  AA[p].chi_required[0] [6] = true;
  AA[p].chi_required[1] [6] = true;

  //bk   chi angles needed for building  CD2
  AA[p].chi_required[0] [7] = true;
  AA[p].chi_required[1] [7] = true;

  //bk   chi angles needed for building  CE1
  AA[p].chi_required[0] [8] = true;
  AA[p].chi_required[1] [8] = true;

  //bk   chi angles needed for building  NE2
  AA[p].chi_required[0] [9] = true;
  AA[p].chi_required[1] [9] = true;

  //bk   chi angles needed for building  HE2
  AA[p].chi_required[0] [11] = true;
  AA[p].chi_required[1] [11] = true;

  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [13] = true;

  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [14] = true;

  //bk   chi angles needed for building  HE1
  AA[p].chi_required[0] [15] = true;
  AA[p].chi_required[1] [15] = true;

  //bk   chi angles needed for building  HD2
  AA[p].chi_required[0] [16] = true;
  AA[p].chi_required[1] [16] = true;

  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   CG
  //bk   four atoms that define chi angle  2
  AA[p].chi_atoms[1] [0] = 1; //   CA
  AA[p].chi_atoms[1] [1] = 4; //   CB
  AA[p].chi_atoms[1] [2] = 5; //   CG:HE1
  AA[p].chi_atoms[1] [3] = 6; //   ND1



  AA[p].atom[0].icoor[0] = 0.0000; //   N
  AA[p].atom[0].icoor[1] = 12.0140; //   N
  AA[p].atom[0].icoor[2] = 22.5220; //   N
  AA[p].atom[1].icoor[0] = 0.0000; //   CA
  AA[p].atom[1].icoor[1] = 13.2300; //   CA
  AA[p].atom[1].icoor[2] = 23.3160; //   CA
  AA[p].atom[2].icoor[0] = 0.0000; //   C
  AA[p].atom[2].icoor[1] = 12.6360; //   C
  AA[p].atom[2].icoor[2] = 24.7310; //   C
  AA[p].atom[3].icoor[0] = 0.1540; //   O
  AA[p].atom[3].icoor[1] = 11.4260; //   O:HE1
  AA[p].atom[3].icoor[2] = 24.9130; //   O
  AA[p].atom[4].icoor[0] = -1.2080; //   CB
  AA[p].atom[4].icoor[1] = 14.1050; //   CB
  AA[p].atom[4].icoor[2] = 22.9660; //   CB
  AA[p].atom[5].icoor[0] = -1.2740; //   CG
  AA[p].atom[5].icoor[1] = 15.3800; //   CG
  AA[p].atom[5].icoor[2] = 23.7480; //   CG
  AA[p].atom[6].icoor[0] = -2.2840; //   ND1
  AA[p].atom[6].icoor[1] = 16.3040; //   ND1
  AA[p].atom[6].icoor[2] = 23.5800; //   ND1
  AA[p].atom[7].icoor[0] = -0.4560; //   CD2
  AA[p].atom[7].icoor[1] = 15.8850; //   CD2
  AA[p].atom[7].icoor[2] = 24.7010; //   CD2
  AA[p].atom[8].icoor[0] = -2.0840; //   CE1
  AA[p].atom[8].icoor[1] = 17.3230; //   CE1:HE1
  AA[p].atom[8].icoor[2] = 24.3980; //   CE1
  AA[p].atom[9].icoor[0] = -0.9820; //   NE2
  AA[p].atom[9].icoor[1] = 17.0930; //   NE2
  AA[p].atom[9].icoor[2] = 25.0880; //   NE2
  AA[p].atom[10].icoor[0] = -0.0410; //   H
  AA[p].atom[10].icoor[1] = 11.1270; //   H
  AA[p].atom[10].icoor[2] = 23.0030; //   H
  AA[p].atom[11].icoor[0] = 0.0; //   HD1 NEW ADDED BUT UNKNOWN
  AA[p].atom[11].icoor[1] = 0.0; //   HD1
  AA[p].atom[11].icoor[2] = 0.0; //   HD1
  AA[p].atom[12].icoor[0] = 0.9100; //   HA
  AA[p].atom[12].icoor[1] = 13.7960; //   HA
  AA[p].atom[12].icoor[2] = 23.1190; //   HA
  AA[p].atom[13].icoor[0] = -1.1760; //  2HB
  AA[p].atom[13].icoor[1] = 14.3870; //  2HB:HE1
  AA[p].atom[13].icoor[2] = 21.9130; //  2HB
  AA[p].atom[14].icoor[0] = -2.1330; //  3HB
  AA[p].atom[14].icoor[1] = 13.5660; //  3HB
  AA[p].atom[14].icoor[2] = 23.1690; //  3HB
  AA[p].atom[15].icoor[0] = -2.7780; //   HE1
  AA[p].atom[15].icoor[1] = 18.1630; //   HE1
  AA[p].atom[15].icoor[2] = 24.4180; //   HE1
  AA[p].atom[16].icoor[0] = 0.4620; //   HD2
  AA[p].atom[16].icoor[1] = 15.5170; //   HD2
  AA[p].atom[16].icoor[2] = 25.1600; //   HD2

  // atom number for backbone HN
  AA[p].HNpos=10;
  // atom number for backbone HA
  AA[p].HApos=12;

  // number of polar hydrogens
  AA[p].nH_polar_complete=2;
  // atom numbers for polar H
  AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
  AA[p].Hpos_polar_complete[0]=10;
  AA[p].Hpos_polar_complete[1]=11;

  // number of aromatic hydrogens
  AA[p].nH_aromatic_complete=0;

  // number of apolar hydrogens
  AA[p].nH_apolar_complete=5;
  //atom number for apolar hydrogens
  AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
  AA[p].Hpos_apolar_complete[0]=12;
  AA[p].Hpos_apolar_complete[1]=13;
  AA[p].Hpos_apolar_complete[2]=14;
  AA[p].Hpos_apolar_complete[3]=15;
  AA[p].Hpos_apolar_complete[4]=16;

  //A REVISAR
  // number of acceptors
  AA[p].nacceptors=2;
  //acceptor information
  AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
  (AA[p].aBase[1])[0]=6;(AA[p].aBase[1])[1]=5;
  (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;
  (AA[p].accpt_pos[1])[0]=1;(AA[p].accpt_pos[1])[1]=6;
  //OJO. HAY UNA RELACION ABASE2(6,8)

  //Number of Hydrogens connection
  AA[p].nH_hydrogen_connexions=7;
  //atoms hydrogens are connected too
  AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
  AA[p].Hbase[0][0]=10; AA[p].Hbase[0][1]=0;
  AA[p].Hbase[1][0]=11; AA[p].Hbase[1][1]=6;
  AA[p].Hbase[2][0]=12; AA[p].Hbase[2][1]=1;
  AA[p].Hbase[3][0]=13; AA[p].Hbase[3][1]=4;
  AA[p].Hbase[4][0]=14; AA[p].Hbase[4][1]=4;
  AA[p].Hbase[5][0]=15; AA[p].Hbase[5][1]=8;
  AA[p].Hbase[6][0]=16; AA[p].Hbase[6][1]=7;

  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=10;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=12;

  AA[p].atom[4].numHydrogens_atm=2;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=13;
  AA[p].atom[4].hydrogens_atm[1]=14;

  AA[p].atom[7].numHydrogens_atm=1;
  AA[p].atom[7].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[7].numHydrogens_atm);
  AA[p].atom[7].hydrogens_atm[0]=16;

  AA[p].atom[8].numHydrogens_atm=1;
  AA[p].atom[8].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[8].numHydrogens_atm);
  AA[p].atom[8].hydrogens_atm[0]=15;

  AA[p].atom[6].numHydrogens_atm=1;
  AA[p].atom[6].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[6].numHydrogens_atm);
  AA[p].atom[6].hydrogens_atm[0]=11;

  //FIN DE AMINOACIDO HISTIDINE NEUTRAL PROTON HD1
}


//NO SE USA Y NO LA HE MODIFICADO = HIE
/**
* Initializes Histidine NEURAL PROTON HE2
*
* @param opt: Rosseta or ICM
*/
void init_HIE(Convention opt)
{
  //AMINOACIDO HISTIDINE
  int k, p =HIE;


  strcpy( AA[p].aa_name3, "HIE" ); // <--- HIE = HIS nuestra = HSD (Protonated at NE2, neutral)
  AA[p].aa_name1 = 'H';            // Note: HID (Protonated at ND1, neutral)
  AA[p].mass = 137.143;            // Note: HIP (Di-Protonated at ND1 & NE2, +)

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = true;
  AA[p].aa_is_nonpolar = false;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = false;

  AA[p].natoms = 17;
  AA[p].nheavyatoms = 10;
  AA[p].nchi = 2;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );
  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;

  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

  AA[p].chi_types = ( int * ) new int[AA[p].nchi];
  AA[p].chi_types[0] = 2;
  AA[p].chi_types[1] = 4;
  AA[p].natoms_EEF1 = 11;


  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " CG " );
  strcpy( AA[p].atom[6].atom_name, " ND1" );
  strcpy( AA[p].atom[7].atom_name, " CD2" );
  strcpy( AA[p].atom[8].atom_name, " CE1" );
  strcpy( AA[p].atom[9].atom_name, " NE2" ); // <--- Protonated here!
  strcpy( AA[p].atom[10].atom_name," H  " );
  strcpy( AA[p].atom[11].atom_name," HE2" );// <--- The proton.
  strcpy( AA[p].atom[12].atom_name," HA " );
  strcpy( AA[p].atom[13].atom_name,"2HB " );
  strcpy( AA[p].atom[14].atom_name,"3HB " );
  strcpy( AA[p].atom[15].atom_name," HE1" );
  strcpy( AA[p].atom[16].atom_name," HD2" );

  strcpy( AA[p].atom[17].atom_name," OXT" );
  strcpy( AA[p].atom[18].atom_name,"2H  " );
  strcpy( AA[p].atom[19].atom_name,"3H  " );


  switch(opt)
  {
    case Rosseta:
    case Sybil:
      AA[p].atom[0].charge = -0.47; //    N
      AA[p].atom[1].charge = 0.07; //    CA
      AA[p].atom[2].charge = 0.51; //    C
      AA[p].atom[3].charge = -0.51; //    O
      AA[p].atom[4].charge = -0.08; //    CB
      AA[p].atom[5].charge = 0.22; //    CG
      AA[p].atom[6].charge = -0.70; //    ND1
      AA[p].atom[7].charge = -0.05; //    CD2
      AA[p].atom[8].charge = 0.25; //    CE1  (¿HE1?)
      AA[p].atom[9].charge = -0.36; //    NE2
      AA[p].atom[10].charge = 0.31; //    H
      AA[p].atom[11].charge = 0.32; //    HE2
      AA[p].atom[12].charge = 0.09; //    HA
      AA[p].atom[13].charge = 0.09; //    2HB
      AA[p].atom[14].charge = 0.09; //    3HB
      AA[p].atom[15].charge = 0.13; //    HE1
      AA[p].atom[16].charge = 0.09; //    HD2

      AA[p].atom[17].charge = -0.67; //   0XT
      AA[p].atom[18].charge = 0.33; //  1H
      AA[p].atom[19].charge = 0.33; //  2H
      break;

    case ICM:
      AA[p].atom[0].charge = -0.356; //    N
      AA[p].atom[1].charge = 0.064; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = -0.05; //    CB
      AA[p].atom[5].charge = 0.15; //    CG
      AA[p].atom[6].charge = -0.2; //    ND1
      AA[p].atom[7].charge = 0.1; //    CD2
      AA[p].atom[8].charge = 0.275; //    CE1
      AA[p].atom[9].charge = -0.2; //    NE2
      AA[p].atom[10].charge = 0.176; //    H
      AA[p].atom[11].charge = 0.265; //    HE2
      AA[p].atom[12].charge = 0.02; //    HA
      AA[p].atom[13].charge = 0.065; //    2HB
      AA[p].atom[14].charge = 0.065; //    3HB
      AA[p].atom[15].charge = 0.15; //    HE1
      AA[p].atom[16].charge = 0.135; //    HD2

      AA[p].atom[17].charge = -0.67; //   0XT
      AA[p].atom[18].charge = 0.33; //  1H
      AA[p].atom[19].charge = 0.33; //  2H
      break;
  case EEF1:
      AA[p].atom[0].charge = -0.35; //    N
      AA[p].atom[1].charge =  0.10; //    CA
      AA[p].atom[2].charge =  0.55; //    C
      AA[p].atom[3].charge = -0.55; //    O
      AA[p].atom[4].charge =  0.00; //    CB
      AA[p].atom[5].charge =  0.10; //    CG
      AA[p].atom[6].charge = -0.40; //    ND1
      AA[p].atom[7].charge =  0.10; //    CD2
      AA[p].atom[8].charge =  0.30; //    CE1
      AA[p].atom[9].charge = -0.40; //    NE2
      AA[p].atom[10].charge = 0.25; //    H
      AA[p].atom[11].charge = 0.30; //    HE2
      AA[p].atom[12].charge = 0.00; //    HA
      AA[p].atom[13].charge = 0.00; //    2HB
      AA[p].atom[14].charge = 0.00; //    3HB
      AA[p].atom[15].charge = 0.00; //    HE1
      AA[p].atom[16].charge = 0.00; //    HD2
      break;
  }
  /*
   1.0 HSE-N:2H    0.33000 HC
   1.0 HSE-N:O    -0.51000 O
   1.0 HSE-N:HA    0.10000 HB
   1.0 HSE-N:3H    0.33000 HC
   1.0 HSE-N:CB   -0.08000 CT2
   1.0 HSE-N:HB1   0.09000 HA
   1.0 HSE-N:HB2   0.09000 HA
   1.0 HSE-N:ND1  -0.70000 NR2
   1.0 HSE-N:CG    0.22000 CPH1
   1.0 HSE-N:CE1   0.25000 CPH2
   1.0 HSE-N:HE1   0.13000 HR1
   1.0 HSE-N:NE2  -0.36000 NR1
   1.0 HSE-N:HE2   0.32000 H
   1.0 HSE-N:CD2  -0.05000 CPH1
   1.0 HSE-N:N    -0.30000 NH3
   1.0 HSE-N:HD2   0.09000 HR3
   1.0 HSE-N:1H    0.33000 HC
   1.0 HSE-N:C     0.51000 C
   1.0 HSE-N:CA    0.21000 CT1

   1.0 HSE-C:HA    0.09000 HB
   1.0 HSE-C:CB   -0.08000 CT2
   1.0 HSE-C:HB1   0.09000 HA
   1.0 HSE-C:HB2   0.09000 HA
   1.0 HSE-C:ND1  -0.70000 NR2
   1.0 HSE-C:CG    0.22000 CPH1
   1.0 HSE-C:CE1   0.25000 CPH2
   1.0 HSE-C:HE1   0.13000 HR1
   1.0 HSE-C:O    -0.67000 OC
   1.0 HSE-C:NE2  -0.36000 NR1
   1.0 HSE-C:OXT  -0.67000 OC
   1.0 HSE-C:HE2   0.32000 H
   1.0 HSE-C:CD2  -0.05000 CPH1
   1.0 HSE-C:N    -0.47000 NH1
   1.0 HSE-C:HD2   0.09000 HR3
   1.0 HSE-C:HN    0.31000 H
   1.0 HSE-C:C     0.34000 CC
   1.0 HSE-C:CA    0.07000 CT1
*/

  switch (opt) {
  	case EEF1:
	  AA[p].atom[0].fullatom_type =  12; // NH1    N
	  AA[p].atom[1].fullatom_type =   5; // CH1E   CA
	  AA[p].atom[2].fullatom_type =   4; // C      C
	  AA[p].atom[3].fullatom_type =  16; // O      O
	  AA[p].atom[4].fullatom_type =   6; // CH2E   CB
	  AA[p].atom[5].fullatom_type =  24; // CR     CG
	  AA[p].atom[6].fullatom_type =  10; // NR     ND1
	  AA[p].atom[7].fullatom_type =   8; // CR1E   CD2
	  AA[p].atom[8].fullatom_type =   8; // CR1E   CE1
	  AA[p].atom[9].fullatom_type =  12; // NH1    NE2
	  AA[p].atom[10].fullatom_type =  0; // H      H
	  AA[p].atom[11].fullatom_type =  0; // H      HE2
	  AA[p].atom[12].fullatom_type = -1; // Hapo   HA
	  AA[p].atom[13].fullatom_type = -1; // Hapo   2HB
	  AA[p].atom[14].fullatom_type = -1; // Hapo   3HB
	  AA[p].atom[15].fullatom_type = -1; // Hapo   HE1
	  AA[p].atom[16].fullatom_type = -1; // Hapo   HD2

	  AA[p].atom[17].fullatom_type = 17; // OC     OXT
	  AA[p].atom[18].fullatom_type =  1; // HC     Hpol
	  AA[p].atom[19].fullatom_type =  1; // HC     Hpol
	  break;
  	case Sybil:

  	  AA[p].atom[0].fullatom_type = 40; // Nbb    N
  	  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
  	  AA[p].atom[2].fullatom_type = 10; // CObb   C
  	  AA[p].atom[3].fullatom_type = 56; // OCbb   O
  	  AA[p].atom[4].fullatom_type = 8; // CH2    CB
  	  AA[p].atom[5].fullatom_type = 11; // aroC   CG
  	  AA[p].atom[6].fullatom_type = 44; // Nhis   ND1
  	  AA[p].atom[7].fullatom_type = 11; // aroC   CD2
  	  AA[p].atom[8].fullatom_type = 11; // aroC   CE1
  	  AA[p].atom[9].fullatom_type = 42; // Ntrp   NE2
	  AA[p].atom[10].fullatom_type = 66; // HNbb   H
	  AA[p].atom[11].fullatom_type = 66; // Hpol   HE2
	  AA[p].atom[12].fullatom_type = 66; // Hapo   HA
	  AA[p].atom[13].fullatom_type = 66; // Hapo  2HB
	  AA[p].atom[14].fullatom_type = 66; // Hapo  3HB
	  AA[p].atom[15].fullatom_type = 66; // Hapo   HE1
	  AA[p].atom[16].fullatom_type = 66; // Hapo   HD2
	  AA[p].atom[17].fullatom_type = 57; // OXT
	  AA[p].atom[18].fullatom_type = 66; // Hpol
	  AA[p].atom[19].fullatom_type = 66; // Hpol
	  break;
	default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 6; // aroC   CG
	  AA[p].atom[6].fullatom_type = 8; // Nhis   ND1
	  AA[p].atom[7].fullatom_type = 6; // aroC   CD2
	  AA[p].atom[8].fullatom_type = 6; // aroC   CE1
	  AA[p].atom[9].fullatom_type = 7; // Ntrp   NE2
	  AA[p].atom[10].fullatom_type = 25; // HNbb   H
	  AA[p].atom[11].fullatom_type = 22; // Hpol   HE2
	  AA[p].atom[12].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[13].fullatom_type = 23; // Hapo  2HB
	  AA[p].atom[14].fullatom_type = 23; // Hapo  3HB
	  AA[p].atom[15].fullatom_type = 23; // Hapo   HE1
	  AA[p].atom[16].fullatom_type = 23; // Hapo   HD2

	  AA[p].atom[17].fullatom_type = 20; // OXT
	  AA[p].atom[18].fullatom_type = 22; // Hpol
	  AA[p].atom[19].fullatom_type = 22; // Hpol
	  break;
  }

  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 10; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 12; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--CG
  AA[p].atom[4].bonded_neighbor[2] = 13; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 14; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 3; // CG
  AA[p].atom[5].bonded_neighbor[0] = 4; // CG--CB
  AA[p].atom[5].bonded_neighbor[1] = 6; // CG--ND1
  AA[p].atom[5].bonded_neighbor[2] = 7; // CG--CD2
  AA[p].atom[6].nbonded_neighbors = 2; // ND1
  AA[p].atom[6].bonded_neighbor[0] = 5; // ND1--CG
  AA[p].atom[6].bonded_neighbor[1] = 8; // ND1--CE1
  AA[p].atom[7].nbonded_neighbors = 3; // CD2
  AA[p].atom[7].bonded_neighbor[0] = 5; // CD2--CG
  AA[p].atom[7].bonded_neighbor[1] = 9; // CD2--NE2
  AA[p].atom[7].bonded_neighbor[2] = 16; // CD2--HD2
  AA[p].atom[8].nbonded_neighbors = 3; // CE1
  AA[p].atom[8].bonded_neighbor[0] = 6; // CE1--ND1
  AA[p].atom[8].bonded_neighbor[1] = 9; // CE1--NE2
  AA[p].atom[8].bonded_neighbor[2] = 15; // CE1--HE1
  AA[p].atom[9].nbonded_neighbors = 3; // NE2
  AA[p].atom[9].bonded_neighbor[0] = 7; // NE2--CD2
  AA[p].atom[9].bonded_neighbor[1] = 8; // NE2--CE1
  AA[p].atom[9].bonded_neighbor[2] = 11; // NE2--HE2
  AA[p].atom[10].nbonded_neighbors = 1; // H
  AA[p].atom[10].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[11].nbonded_neighbors = 1; // HE2
  AA[p].atom[11].bonded_neighbor[0] = 9; // HE2--NE2
  AA[p].atom[12].nbonded_neighbors = 1; // HA
  AA[p].atom[12].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[13].nbonded_neighbors = 1; //2HB
  AA[p].atom[13].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[14].nbonded_neighbors = 1; //3HB
  AA[p].atom[14].bonded_neighbor[0] = 4; //3HB--CB
  AA[p].atom[15].nbonded_neighbors = 1; // HE1
  AA[p].atom[15].bonded_neighbor[0] = 8; // HE1--CE1
  AA[p].atom[16].nbonded_neighbors = 1; // HD2
  AA[p].atom[16].bonded_neighbor[0] = 7; // HD2--CD2

  AA[p].atom[17].nbonded_neighbors = 1; // OXT
  AA[p].atom[17].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[18].nbonded_neighbors = 1; // 1H
  AA[p].atom[18].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[19].nbonded_neighbors = 1; // 2H
  AA[p].atom[19].bonded_neighbor[0] = 0; // H-N


  for ( k = 0; k < 11; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 11; k < 17; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template for building  HE2
  AA[p].atom[11].ta[0] = 9; //   NE2
  AA[p].atom[11].ta[1] = 7; //   CD2
  AA[p].atom[11].ta[2] = 8; //   CE1

  //bk   template for building  HA
  AA[p].atom[12].ta[0] = 1; //   CA
  AA[p].atom[12].ta[1] = 0; //   N
  AA[p].atom[12].ta[2] = 2; //   C

  //bk   template for building  2HB
  AA[p].atom[13].ta[0] = 4; //   CB
  AA[p].atom[13].ta[1] = 1; //   CA
  AA[p].atom[13].ta[2] = 5; //   CG

  //bk   template for building  3HB
  AA[p].atom[14].ta[0] = 4; //   CB
  AA[p].atom[14].ta[1] = 1; //   CA
  AA[p].atom[14].ta[2] = 5; //   CG

  //bk   template for building  HE1
  AA[p].atom[15].ta[0] = 8; //   CE1
  AA[p].atom[15].ta[1] = 6; //   ND1
  AA[p].atom[15].ta[2] = 9; //  NE2

  //bk   template for building  HD2
  AA[p].atom[16].ta[0] = 7; //   CD2
  AA[p].atom[16].ta[1] = 5; //   CG
  AA[p].atom[16].ta[2] = 9; //   NE2




  //bk   chi angles required to build atoms HIS
  //bk   chi angles needed for building  CG:HE1
  AA[p].chi_required[0] [5] = true;

  //bk   chi angles needed for building  ND1
  AA[p].chi_required[0] [6] = true;
  AA[p].chi_required[1] [6] = true;

  //bk   chi angles needed for building  CD2
  AA[p].chi_required[0] [7] = true;
  AA[p].chi_required[1] [7] = true;

  //bk   chi angles needed for building  CE1
  AA[p].chi_required[0] [8] = true;
  AA[p].chi_required[1] [8] = true;

  //bk   chi angles needed for building  NE2
  AA[p].chi_required[0] [9] = true;
  AA[p].chi_required[1] [9] = true;

  //bk   chi angles needed for building  HE2
  AA[p].chi_required[0] [11] = true;
  AA[p].chi_required[1] [11] = true;

  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [13] = true;

  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [14] = true;

  //bk   chi angles needed for building  HE1
  AA[p].chi_required[0] [15] = true;
  AA[p].chi_required[1] [15] = true;

  //bk   chi angles needed for building  HD2
  AA[p].chi_required[0] [16] = true;
  AA[p].chi_required[1] [16] = true;

  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   CG
  //bk   four atoms that define chi angle  2
  AA[p].chi_atoms[1] [0] = 1; //   CA
  AA[p].chi_atoms[1] [1] = 4; //   CB
  AA[p].chi_atoms[1] [2] = 5; //   CG:HE1
  AA[p].chi_atoms[1] [3] = 6; //   ND1



  AA[p].atom[0].icoor[0] = 0.0000; //   N
  AA[p].atom[0].icoor[1] = 12.0140; //   N
  AA[p].atom[0].icoor[2] = 22.5220; //   N
  AA[p].atom[1].icoor[0] = 0.0000; //   CA
  AA[p].atom[1].icoor[1] = 13.2300; //   CA
  AA[p].atom[1].icoor[2] = 23.3160; //   CA
  AA[p].atom[2].icoor[0] = 0.0000; //   C
  AA[p].atom[2].icoor[1] = 12.6360; //   C
  AA[p].atom[2].icoor[2] = 24.7310; //   C
  AA[p].atom[3].icoor[0] = 0.1540; //   O
  AA[p].atom[3].icoor[1] = 11.4260; //   O:HE1
  AA[p].atom[3].icoor[2] = 24.9130; //   O
  AA[p].atom[4].icoor[0] = -1.2080; //   CB
  AA[p].atom[4].icoor[1] = 14.1050; //   CB
  AA[p].atom[4].icoor[2] = 22.9660; //   CB
  AA[p].atom[5].icoor[0] = -1.2740; //   CG
  AA[p].atom[5].icoor[1] = 15.3800; //   CG
  AA[p].atom[5].icoor[2] = 23.7480; //   CG
  AA[p].atom[6].icoor[0] = -2.2840; //   ND1
  AA[p].atom[6].icoor[1] = 16.3040; //   ND1
  AA[p].atom[6].icoor[2] = 23.5800; //   ND1
  AA[p].atom[7].icoor[0] = -0.4560; //   CD2
  AA[p].atom[7].icoor[1] = 15.8850; //   CD2
  AA[p].atom[7].icoor[2] = 24.7010; //   CD2
  AA[p].atom[8].icoor[0] = -2.0840; //   CE1
  AA[p].atom[8].icoor[1] = 17.3230; //   CE1:HE1
  AA[p].atom[8].icoor[2] = 24.3980; //   CE1
  AA[p].atom[9].icoor[0] = -0.9820; //   NE2
  AA[p].atom[9].icoor[1] = 17.0930; //   NE2
  AA[p].atom[9].icoor[2] = 25.0880; //   NE2
  AA[p].atom[10].icoor[0] = -0.0410; //   H
  AA[p].atom[10].icoor[1] = 11.1270; //   H
  AA[p].atom[10].icoor[2] = 23.0030; //   H
  AA[p].atom[11].icoor[0] = -0.5830; //   HE2
  AA[p].atom[11].icoor[1] = 17.7010; //   HE2
  AA[p].atom[11].icoor[2] = 25.7890; //   HE2
  AA[p].atom[12].icoor[0] = 0.9100; //   HA
  AA[p].atom[12].icoor[1] = 13.7960; //   HA
  AA[p].atom[12].icoor[2] = 23.1190; //   HA
  AA[p].atom[13].icoor[0] = -1.1760; //  2HB
  AA[p].atom[13].icoor[1] = 14.3870; //  2HB:HE1
  AA[p].atom[13].icoor[2] = 21.9130; //  2HB
  AA[p].atom[14].icoor[0] = -2.1330; //  3HB
  AA[p].atom[14].icoor[1] = 13.5660; //  3HB
  AA[p].atom[14].icoor[2] = 23.1690; //  3HB
  AA[p].atom[15].icoor[0] = -2.7780; //   HE1
  AA[p].atom[15].icoor[1] = 18.1630; //   HE1
  AA[p].atom[15].icoor[2] = 24.4180; //   HE1
  AA[p].atom[16].icoor[0] = 0.4620; //   HD2
  AA[p].atom[16].icoor[1] = 15.5170; //   HD2
  AA[p].atom[16].icoor[2] = 25.1600; //   HD2

  // atom number for backbone HN
  AA[p].HNpos=10;
  // atom number for backbone HA
  AA[p].HApos=12;

  // number of polar hydrogens
  AA[p].nH_polar_complete=2;
  // atom numbers for polar H
  AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
  AA[p].Hpos_polar_complete[0]=10;
  AA[p].Hpos_polar_complete[1]=11;

  // number of aromatic hydrogens
  AA[p].nH_aromatic_complete=0;

  // number of apolar hydrogens
  AA[p].nH_apolar_complete=5;
  //atom number for apolar hydrogens
  AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
  AA[p].Hpos_apolar_complete[0]=12;
  AA[p].Hpos_apolar_complete[1]=13;
  AA[p].Hpos_apolar_complete[2]=14;
  AA[p].Hpos_apolar_complete[3]=15;
  AA[p].Hpos_apolar_complete[4]=16;

  // number of acceptors
  AA[p].nacceptors=2;
  //acceptor information
  AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
  (AA[p].aBase[1])[0]=6;(AA[p].aBase[1])[1]=5;
  (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;
  (AA[p].accpt_pos[1])[0]=1;(AA[p].accpt_pos[1])[1]=6;
  //OJO. HAY UNA RELACION ABASE2(6,8)

  //Number of Hydrogens connection
  AA[p].nH_hydrogen_connexions=7;
  //atoms hydrogens are connected too
  AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
  AA[p].Hbase[0][0]=10; AA[p].Hbase[0][1]=0;
  AA[p].Hbase[1][0]=11; AA[p].Hbase[1][1]=9;
  AA[p].Hbase[2][0]=12; AA[p].Hbase[2][1]=1;
  AA[p].Hbase[3][0]=13; AA[p].Hbase[3][1]=4;
  AA[p].Hbase[4][0]=14; AA[p].Hbase[4][1]=4;
  AA[p].Hbase[5][0]=15; AA[p].Hbase[5][1]=8;
  AA[p].Hbase[6][0]=16; AA[p].Hbase[6][1]=7;

  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=10;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=12;

  AA[p].atom[4].numHydrogens_atm=2;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=13;
  AA[p].atom[4].hydrogens_atm[1]=14;

  AA[p].atom[7].numHydrogens_atm=1;
  AA[p].atom[7].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[7].numHydrogens_atm);
  AA[p].atom[7].hydrogens_atm[0]=16;

  AA[p].atom[8].numHydrogens_atm=1;
  AA[p].atom[8].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[8].numHydrogens_atm);
  AA[p].atom[8].hydrogens_atm[0]=15;

  AA[p].atom[9].numHydrogens_atm=1;
  AA[p].atom[9].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[9].numHydrogens_atm);
  AA[p].atom[9].hydrogens_atm[0]=11;

  //FIN DE AMINOACIDO HISTIDINE PROTON HE2
}

/**
* Initializes Lysine NEUTRAL
*
* @param opt: Rosseta or ICM
*/
void init_LYN(Convention opt)
{
  //AMINOACIDO LYSINA NEUTRAL
  int k, p =LYN;

  strcpy( AA[p].aa_name3, "LYN" );
  AA[p].aa_name1 = 'K';
  AA[p].mass = 129.184;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = true;
  AA[p].aa_is_nonpolar = false;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = true;

  AA[p].natoms = 21; // NEW ADDED
  AA[p].nheavyatoms = 9;
  AA[p].nchi = 4;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );
  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;


  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

AA[p].chi_types = ( int * ) new int[AA[p].nchi];
  AA[p].chi_types[0] = 2;
  AA[p].chi_types[1] = 2;
  AA[p].chi_types[2] = 2;
  AA[p].chi_types[3] = 2;


  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " CG " );
  strcpy( AA[p].atom[6].atom_name, " CD " );
  strcpy( AA[p].atom[7].atom_name, " CE " );
  strcpy( AA[p].atom[8].atom_name, " NZ " );
  strcpy( AA[p].atom[9].atom_name, " H  " );
  strcpy( AA[p].atom[10].atom_name, "1HZ " );
  strcpy( AA[p].atom[11].atom_name, "2HZ " );
  strcpy( AA[p].atom[12].atom_name, " HA " );
  strcpy( AA[p].atom[13].atom_name, "2HB " );
  strcpy( AA[p].atom[14].atom_name, "3HB " );
  strcpy( AA[p].atom[15].atom_name, "2HG " );
  strcpy( AA[p].atom[16].atom_name, "3HG " );
  strcpy( AA[p].atom[17].atom_name, "2HD " );
  strcpy( AA[p].atom[18].atom_name, "3HD " );
  strcpy( AA[p].atom[19].atom_name, "2HE " );
  strcpy( AA[p].atom[20].atom_name, "3HE " );

  strcpy( AA[p].atom[21].atom_name," OXT" );
  strcpy( AA[p].atom[22].atom_name,"2H  " );
  strcpy( AA[p].atom[23].atom_name,"3H  " );

  switch(opt)
  {
    //A REVISAR
    case Rosseta:
    case Sybil:
      AA[p].atom[0].charge = -0.47; //    N
      AA[p].atom[1].charge = 0.07; //    CA
      AA[p].atom[2].charge = 0.51; //    C
      AA[p].atom[3].charge = -0.51; //    O
      AA[p].atom[4].charge = -0.18; //    CB
      AA[p].atom[5].charge = -0.18; //    CG
      AA[p].atom[6].charge = -0.18; //    CD
      AA[p].atom[7].charge = 0.21; //    CE
      AA[p].atom[8].charge = -0.30; //    NZ
      AA[p].atom[9].charge = 0.31; //    H
      AA[p].atom[10].charge = 0.33; //    1HZ
      AA[p].atom[11].charge = 0.33; //    2HZ
      AA[p].atom[12].charge = 0.09; //    HA
      AA[p].atom[13].charge = 0.09; //    2HB
      AA[p].atom[14].charge = 0.09; //    3HB
      AA[p].atom[15].charge = 0.09; //    2HG
      AA[p].atom[16].charge = 0.09; //    3HG
      AA[p].atom[17].charge = 0.09; //    2HD
      AA[p].atom[18].charge = 0.09; //    3HD
      AA[p].atom[19].charge = 0.05; //    2HE
      AA[p].atom[20].charge = 0.05; //    3HE

      AA[p].atom[21].charge = -0.67; //   0XT
      AA[p].atom[22].charge = 0.33; //  1H
      AA[p].atom[23].charge = 0.33; //  2H
      break;

    //A REVISAR
    case ICM:
      AA[p].atom[0].charge = -0.356; //    N
      AA[p].atom[1].charge = 0.064; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = -0.03; //    CB
      AA[p].atom[5].charge = -0.025; //    CG
      AA[p].atom[6].charge = -0.115; //    CD
      AA[p].atom[7].charge = 0.05; //    CE
      AA[p].atom[8].charge = -0.32; //    NZ
      AA[p].atom[9].charge = 0.176; //    H
      AA[p].atom[10].charge = 0.32; //    1HZ
      AA[p].atom[11].charge = 0.32; //    2HZ
      AA[p].atom[12].charge = 0.02; //    HA
      AA[p].atom[13].charge = 0.015; //    2HB
      AA[p].atom[14].charge = 0.015; //    3HB
      AA[p].atom[15].charge = 0.02; //    2HG
      AA[p].atom[16].charge = 0.02; //    3HG
      AA[p].atom[17].charge = 0.1; //    2HD
      AA[p].atom[18].charge = 0.1; //    3HD
      AA[p].atom[19].charge = 0.12; //    2HE
      AA[p].atom[20].charge = 0.12; //    3HE

      AA[p].atom[21].charge = -0.67; //   0XT
      AA[p].atom[22].charge = 0.33; //  1H
      AA[p].atom[23].charge = 0.33; //  2H
      break;
	default:
	  break;
  }


  /*
   1.0 LYS-N:2H    0.33000 HC
   1.0 LYS-N:NZ   -0.30000 NH3
   1.0 LYS-N:HA    0.10000 HB
   1.0 LYS-N:3H    0.33000 HC
   1.0 LYS-N:1HZ   0.33000 HC
   1.0 LYS-N:CB   -0.18000 CT2
   1.0 LYS-N:2HZ   0.33000 HC
   1.0 LYS-N:HB1   0.09000 HA
   1.0 LYS-N:1HB   0.09000 HA
   1.0 LYS-N:CG   -0.18000 CT2
   1.0 LYS-N:HG1   0.09000 HA
   1.0 LYS-N:1HG   0.09000 HA
   1.0 LYS-N:3HZ   0.33000 HC
   1.0 LYS-N:CD   -0.18000 CT2
   1.0 LYS-N:C     0.51000 C
   1.0 LYS-N:HD1   0.09000 HA
   1.0 LYS-N:O    -0.51000 O
   1.0 LYS-N:1HD   0.09000 HA
   1.0 LYS-N:CE    0.21000 CT2
   1.0 LYS-N:N    -0.30000 NH3
   1.0 LYS-N:HE1   0.05000 HA
   1.0 LYS-N:1H    0.33000 HC
   1.0 LYS-N:1HE   0.05000 HA
   1.0 LYS-N:CA    0.21000 CT1

   1.0 LYS-C:NZ   -0.30000 NH3
   1.0 LYS-C:HA    0.09000 HB
   1.0 LYS-C:1HZ   0.33000 HC
   1.0 LYS-C:CB   -0.18000 CT2
   1.0 LYS-C:2HZ   0.33000 HC
   1.0 LYS-C:HB1   0.09000 HA
   1.0 LYS-C:1HB   0.09000 HA
   1.0 LYS-C:CG   -0.18000 CT2
   1.0 LYS-C:HG1   0.09000 HA
   1.0 LYS-C:1HG   0.09000 HA
   1.0 LYS-C:3HZ   0.33000 HC
   1.0 LYS-C:CD   -0.18000 CT2
   1.0 LYS-C:O    -0.67000 OC
   1.0 LYS-C:C     0.34000 CC
   1.0 LYS-C:HD1   0.09000 HA
   1.0 LYS-C:OXT  -0.67000 OC
   1.0 LYS-C:1HD   0.09000 HA
   1.0 LYS-C:CE    0.21000 CT2
   1.0 LYS-C:N    -0.47000 NH1
   1.0 LYS-C:HE1   0.05000 HA
   1.0 LYS-C:HN    0.31000 H
   1.0 LYS-C:1HE   0.05000 HA
   1.0 LYS-C:CA    0.07000 CT1
*/
  switch (opt) {
  case Sybil:

	  AA[p].atom[0].fullatom_type = 40; // Nbb    N
	  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 10; // CObb   C
	  AA[p].atom[3].fullatom_type = 56; // OCbb   O
	  AA[p].atom[4].fullatom_type = 8; // CH2    CB
	  AA[p].atom[5].fullatom_type = 8; // CH2    CG
	  AA[p].atom[6].fullatom_type = 8; // CH2    CD
	  AA[p].atom[7].fullatom_type = 5; // CH2    CE
	  AA[p].atom[8].fullatom_type = 35; // Nlys   NZ
	  AA[p].atom[9].fullatom_type = 66; // HNbb   H
	  AA[p].atom[10].fullatom_type = 66; // Hpol  1HZ
	  AA[p].atom[11].fullatom_type = 66; // Hpol  2HZ
	  AA[p].atom[12].fullatom_type = 66; // Hapo   HA
	  AA[p].atom[13].fullatom_type = 66; // Hapo  2HB
	  AA[p].atom[14].fullatom_type = 66; // Hapo  3HB
	  AA[p].atom[15].fullatom_type = 66; // Hapo  2HG
	  AA[p].atom[16].fullatom_type = 66; // Hapo  3HG
	  AA[p].atom[17].fullatom_type = 66; // Hapo  2HD
	  AA[p].atom[18].fullatom_type = 66; // Hapo  3HD
	  AA[p].atom[19].fullatom_type = 66; // Hapo  2HE
	  AA[p].atom[20].fullatom_type = 66; // Hapo  3HE

	  AA[p].atom[21].fullatom_type = 57; // OXT
	  AA[p].atom[22].fullatom_type = 66; // Hpol
	  AA[p].atom[23].fullatom_type = 66; // Hpol
	  break;
  default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 4; // CH2    CG
	  AA[p].atom[6].fullatom_type = 4; // CH2    CD
	  AA[p].atom[7].fullatom_type = 4; // CH2    CE
	  AA[p].atom[8].fullatom_type = 10; // Nlys   NZ
	  AA[p].atom[9].fullatom_type = 25; // HNbb   H
	  AA[p].atom[10].fullatom_type = 22; // Hpol  1HZ
	  AA[p].atom[11].fullatom_type = 22; // Hpol  2HZ
	  AA[p].atom[12].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[13].fullatom_type = 23; // Hapo  2HB
	  AA[p].atom[14].fullatom_type = 23; // Hapo  3HB
	  AA[p].atom[15].fullatom_type = 23; // Hapo  2HG
	  AA[p].atom[16].fullatom_type = 23; // Hapo  3HG
	  AA[p].atom[17].fullatom_type = 23; // Hapo  2HD
	  AA[p].atom[18].fullatom_type = 23; // Hapo  3HD
	  AA[p].atom[19].fullatom_type = 23; // Hapo  2HE
	  AA[p].atom[20].fullatom_type = 23; // Hapo  3HE

	  AA[p].atom[21].fullatom_type = 20; // OXT
	  AA[p].atom[22].fullatom_type = 22; // Hpol
	  AA[p].atom[23].fullatom_type = 22; // Hpol
	  break;
  }


  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 9; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 12; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--CG
  AA[p].atom[4].bonded_neighbor[2] = 13; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 14; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 4; // CG
  AA[p].atom[5].bonded_neighbor[0] = 4; // CG--CB
  AA[p].atom[5].bonded_neighbor[1] = 6; // CG--CD
  AA[p].atom[5].bonded_neighbor[2] = 15; // CG--2HG
  AA[p].atom[5].bonded_neighbor[3] = 16; // CG--3HG
  AA[p].atom[6].nbonded_neighbors = 4; // CD
  AA[p].atom[6].bonded_neighbor[0] = 5; // CD--CG
  AA[p].atom[6].bonded_neighbor[1] = 7; // CD--CE
  AA[p].atom[6].bonded_neighbor[2] = 17; // CD--2HD
  AA[p].atom[6].bonded_neighbor[3] = 18; // CD--3HD
  AA[p].atom[7].nbonded_neighbors = 4; // CE
  AA[p].atom[7].bonded_neighbor[0] = 6; // CE--CD
  AA[p].atom[7].bonded_neighbor[1] = 8; // CE--NZ
  AA[p].atom[7].bonded_neighbor[2] = 19; // CE--2HE
  AA[p].atom[7].bonded_neighbor[3] = 20; // CE--3HE
  AA[p].atom[8].nbonded_neighbors = 3; // NZ
  AA[p].atom[8].bonded_neighbor[0] = 7; // NZ--CE
  AA[p].atom[8].bonded_neighbor[1] = 10; // NZ--1HZ
  AA[p].atom[8].bonded_neighbor[2] = 11; // NZ--2HZ
  AA[p].atom[9].nbonded_neighbors = 1; // H
  AA[p].atom[9].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[10].nbonded_neighbors = 1; //1HZ
  AA[p].atom[10].bonded_neighbor[0] = 8; //1HZ--NZ
  AA[p].atom[11].nbonded_neighbors = 1; //2HZ
  AA[p].atom[11].bonded_neighbor[0] = 8; //2HZ--NZ
  AA[p].atom[12].nbonded_neighbors = 1; // HA
  AA[p].atom[12].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[13].nbonded_neighbors = 1; //2HB
  AA[p].atom[13].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[14].nbonded_neighbors = 1; //3HB
  AA[p].atom[14].bonded_neighbor[0] = 4; //3HB--CB
  AA[p].atom[15].nbonded_neighbors = 1; //2HG
  AA[p].atom[15].bonded_neighbor[0] = 5; //2HG--CG
  AA[p].atom[16].nbonded_neighbors = 1; //3HG
  AA[p].atom[16].bonded_neighbor[0] = 5; //3HG--CG
  AA[p].atom[17].nbonded_neighbors = 1; //2HD
  AA[p].atom[17].bonded_neighbor[0] = 6; //2HD--CD
  AA[p].atom[18].nbonded_neighbors = 1; //3HD
  AA[p].atom[18].bonded_neighbor[0] = 6; //3HD--CD
  AA[p].atom[19].nbonded_neighbors = 1; //2HE
  AA[p].atom[19].bonded_neighbor[0] = 7; //2HE--CE
  AA[p].atom[20].nbonded_neighbors = 1; //3HE
  AA[p].atom[20].bonded_neighbor[0] = 7; //3HE--CE

  AA[p].atom[22].nbonded_neighbors = 1; // OXT
  AA[p].atom[22].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[23].nbonded_neighbors = 1; // 1H
  AA[p].atom[23].bonded_neighbor[0] = 0; // H-N
// Bug 23/11/2009
//  AA[p].atom[24].nbonded_neighbors = 1; // 2H
//  AA[p].atom[24].bonded_neighbor[0] = 0; // H-N


  for ( k = 0; k < 10; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 10; k < 21; k++ )
    AA[p].atom[k].hastemplate = true;

	//A REVISAR TEMPLATES Y CHIS
  //bk   template for building 1HZ
  AA[p].atom[10].ta[0] = 8; //   NZ
  AA[p].atom[10].ta[1] = 7; //   CE
  AA[p].atom[10].ta[2] = 6; //   CD

  //bk   template for building 2HZ
  AA[p].atom[11].ta[0] = 8; //   NZ
  AA[p].atom[11].ta[1] = 7; //   CE
  AA[p].atom[11].ta[2] = 6; //   CD


  //bk   template for building  HA
  AA[p].atom[12].ta[0] = 1; //   CA
  AA[p].atom[12].ta[1] = 0; //   N
  AA[p].atom[12].ta[2] = 2; //   C

  //bk   template for building 2HB
  AA[p].atom[13].ta[0] = 4; //   CB
  AA[p].atom[13].ta[1] = 1; //   CA
  AA[p].atom[13].ta[2] = 5; //   CG

  //bk   template for building 3HB
  AA[p].atom[14].ta[0] = 4; //   CB
  AA[p].atom[14].ta[1] = 1; //   CA
  AA[p].atom[14].ta[2] = 5; //   CG

  //bk   template for building 2HG
  AA[p].atom[15].ta[0] = 5; //   CG
  AA[p].atom[15].ta[1] = 4; //   CB
  AA[p].atom[15].ta[2] = 6; //   CD

  //bk   template for building 3HG
  AA[p].atom[16].ta[0] = 5; //   CG
  AA[p].atom[16].ta[1] = 4; //   CB
  AA[p].atom[16].ta[2] = 6; //   CD

  //bk   template for building 2HD
  AA[p].atom[17].ta[0] = 6; //   CD
  AA[p].atom[17].ta[1] = 5; //   CG
  AA[p].atom[17].ta[2] = 7; //   CE

  //bk   template for building 3HD
  AA[p].atom[18].ta[0] = 6; //   CD
  AA[p].atom[18].ta[1] = 5; //   CG
  AA[p].atom[18].ta[2] = 7; //   CE

  //bk   template for building 2HE
  AA[p].atom[19].ta[0] = 7; //   CE
  AA[p].atom[19].ta[1] = 6; //   CD
  AA[p].atom[19].ta[2] = 8; //   NZ

  //bk   template for building 3HE
  AA[p].atom[20].ta[0] = 7; //   CE
  AA[p].atom[20].ta[1] = 6; //   CD
  AA[p].atom[20].ta[2] = 8; //   NZ


  //bk   chi angles required to build atoms LYS

  //bk   chi angles needed for building  N
  //bk   chi angles needed for building  CA
  //bk   chi angles needed for building  C
  //bk   chi angles needed for building  O
  //bk   chi angles needed for building  CB
  //bk   chi angles needed for building  CG
  AA[p].chi_required[0] [5] = true;


  //bk   chi angles needed for building  CD
  AA[p].chi_required[0] [6] = true;
  AA[p].chi_required[1] [6] = true;

  //bk   chi angles needed for building  CE
  AA[p].chi_required[0] [7] = true;
  AA[p].chi_required[1] [7] = true;
  AA[p].chi_required[2] [7] = true;

  //bk   chi angles needed for building  NZ
  AA[p].chi_required[0] [8] = true;
  AA[p].chi_required[1] [8] = true;
  AA[p].chi_required[2] [8] = true;
  AA[p].chi_required[3] [8] = true;

  //bk   chi angles needed for building  H

  //bk   chi angles needed for building 1HZ
  AA[p].chi_required[0] [10] = true;
  AA[p].chi_required[1] [10] = true;
  AA[p].chi_required[2] [10] = true;
  AA[p].chi_required[3] [10] = true;

  //bk   chi angles needed for building 2HZ
  AA[p].chi_required[0] [11] = true;
  AA[p].chi_required[1] [11] = true;
  AA[p].chi_required[2] [11] = true;
  AA[p].chi_required[3] [11] = true;


  //bk   chi angles needed for building  HA

  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [13] = true;

  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [14] = true;

  //bk   chi angles needed for building 2HG
  AA[p].chi_required[0] [15] = true;
  AA[p].chi_required[1] [15] = true;


  //bk   chi angles needed for building 3HG
  AA[p].chi_required[0] [16] = true;
  AA[p].chi_required[1] [16] = true;


  //bk   chi angles needed for building 2HD
  AA[p].chi_required[0] [17] = true;
  AA[p].chi_required[1] [17] = true;
  AA[p].chi_required[2] [17] = true;


  //bk   chi angles needed for building 3HD
  AA[p].chi_required[0] [18] = true;
  AA[p].chi_required[1] [18] = true;
  AA[p].chi_required[2] [18] = true;


  //bk   chi angles needed for building 2HE
  AA[p].chi_required[0] [19] = true;
  AA[p].chi_required[1] [19] = true;
  AA[p].chi_required[2] [19] = true;
  AA[p].chi_required[3] [19] = true;

  //bk   chi angles needed for building 3HE
  AA[p].chi_required[0] [20] = true;
  AA[p].chi_required[1] [20] = true;
  AA[p].chi_required[2] [20] = true;
  AA[p].chi_required[3] [20] = true;


  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   CG
  //bk   four atoms that define chi angle  2
  AA[p].chi_atoms[1] [0] = 1; //   CA
  AA[p].chi_atoms[1] [1] = 4; //   CB
  AA[p].chi_atoms[1] [2] = 5; //   CG
  AA[p].chi_atoms[1] [3] = 6; //   CD
  //bk   four atoms that define chi angle  3
  AA[p].chi_atoms[2] [0] = 4; //   CB
  AA[p].chi_atoms[2] [1] = 5; //   CG
  AA[p].chi_atoms[2] [2] = 6; //   CD
  AA[p].chi_atoms[2] [3] = 7; //   CE
  //bk   four atoms that define chi angle  4
  AA[p].chi_atoms[3] [0] = 5; //   CG
  AA[p].chi_atoms[3] [1] = 6; //   CD
  AA[p].chi_atoms[3] [2] = 7; //   CE
  AA[p].chi_atoms[3] [3] = 8; //   NZ


  AA[p].atom[0].icoor[0] = 0.0000; //   N
  AA[p].atom[0].icoor[1] = 13.9730; //   N
  AA[p].atom[0].icoor[2] = 29.3820; //   N
  AA[p].atom[1].icoor[0] = 0.0000; //   CA
  AA[p].atom[1].icoor[1] = 15.0130; //   CA
  AA[p].atom[1].icoor[2] = 30.4040; //   CA
  AA[p].atom[2].icoor[0] = 0.0000; //   C
  AA[p].atom[2].icoor[1] = 14.4110; //   C
  AA[p].atom[2].icoor[2] = 31.8030; //   C
  AA[p].atom[3].icoor[0] = 0.0030; //   O
  AA[p].atom[3].icoor[1] = 13.1890; //   O
  AA[p].atom[3].icoor[2] = 31.9540; //   O
  AA[p].atom[4].icoor[0] = -1.2070; //   CB
  AA[p].atom[4].icoor[1] = 15.9360; //   CB
  AA[p].atom[4].icoor[2] = 30.2290; //   CB
  AA[p].atom[5].icoor[0] = -1.2940; //   CG
  AA[p].atom[5].icoor[1] = 17.0590; //   CG
  AA[p].atom[5].icoor[2] = 31.2540; //   CG
  AA[p].atom[6].icoor[0] = -2.5170; //   CD
  AA[p].atom[6].icoor[1] = 17.9320; //   CD
  AA[p].atom[6].icoor[2] = 31.0160; //   CD
  AA[p].atom[7].icoor[0] = -2.6050; //   CE
  AA[p].atom[7].icoor[1] = 19.0540; //   CE
  AA[p].atom[7].icoor[2] = 32.0400; //   CE
  AA[p].atom[8].icoor[0] = -3.7970; //   NZ
  AA[p].atom[8].icoor[1] = 19.9170; //   NZ
  AA[p].atom[8].icoor[2] = 31.8190; //   NZ
  AA[p].atom[9].icoor[0] = -0.0400; //   H
  AA[p].atom[9].icoor[1] = 13.0090; //   H
  AA[p].atom[9].icoor[2] = 29.6800; //   H
  AA[p].atom[10].icoor[0] = -3.8180; //  1HZ
  AA[p].atom[10].icoor[1] = 20.6470; //  1HZ
  AA[p].atom[10].icoor[2] = 32.5170; //  1HZ
  AA[p].atom[11].icoor[0] = -3.7470; //  2HZ
  AA[p].atom[11].icoor[1] = 20.3320; //  2HZ
  AA[p].atom[11].icoor[2] = 30.8990; //  2HZ
  AA[p].atom[12].icoor[0] = 0.9090; //   HA
  AA[p].atom[12].icoor[1] = 15.6100; //   HA
  AA[p].atom[12].icoor[2] = 30.3230; //   HA
  AA[p].atom[13].icoor[0] = -1.1400; //  2HB
  AA[p].atom[13].icoor[1] = 16.3630; //  2HB
  AA[p].atom[13].icoor[2] = 29.2280; //  2HB
  AA[p].atom[14].icoor[0] = -2.0990; //  3HB
  AA[p].atom[14].icoor[1] = 15.3130; //  3HB
  AA[p].atom[14].icoor[2] = 30.2980; //  3HB
  AA[p].atom[15].icoor[0] = -1.3510; //  2HG
  AA[p].atom[15].icoor[1] = 16.6160; //  2HG
  AA[p].atom[15].icoor[2] = 32.2490; //  2HG
  AA[p].atom[16].icoor[0] = -0.3930; //  3HG
  AA[p].atom[16].icoor[1] = 17.6670; //  3HG
  AA[p].atom[16].icoor[2] = 31.1790; //  3HG
  AA[p].atom[17].icoor[0] = -2.4500; //  2HD
  AA[p].atom[17].icoor[1] = 18.3580; //  2HD
  AA[p].atom[17].icoor[2] = 30.0140; //  2HD
  AA[p].atom[18].icoor[0] = -3.4080; //  3HD
  AA[p].atom[18].icoor[1] = 17.3080; //  3HD
  AA[p].atom[18].icoor[2] = 31.0850; //  3HD
  AA[p].atom[19].icoor[0] = -2.6590; //  2HE
  AA[p].atom[19].icoor[1] = 18.6080; //  2HE
  AA[p].atom[19].icoor[2] = 33.0320; //  2HE
  AA[p].atom[20].icoor[0] = -1.7010; //  3HE
  AA[p].atom[20].icoor[1] = 19.6580; //  3HE
  AA[p].atom[20].icoor[2] = 31.9630; //  3HE


  // atom number for backbone HN
  AA[p].HNpos=9;
  // atom number for backbone HA
  AA[p].HApos=12;

  // number of polar hydrogens
  AA[p].nH_polar_complete=3;
  // atom numbers for polar H
  AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
  AA[p].Hpos_polar_complete[0]=9;
  AA[p].Hpos_polar_complete[1]=10;
  AA[p].Hpos_polar_complete[2]=11;

  // number of aromatic hydrogens
  AA[p].nH_aromatic_complete=0;

  // number of apolar hydrogens
  AA[p].nH_apolar_complete=9;
  //atom number for apolar hydrogens
  AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
  AA[p].Hpos_apolar_complete[0]=12;
  AA[p].Hpos_apolar_complete[1]=13;
  AA[p].Hpos_apolar_complete[2]=14;
  AA[p].Hpos_apolar_complete[3]=15;
  AA[p].Hpos_apolar_complete[4]=16;
  AA[p].Hpos_apolar_complete[5]=17;
  AA[p].Hpos_apolar_complete[6]=18;
  AA[p].Hpos_apolar_complete[7]=19;
  AA[p].Hpos_apolar_complete[8]=20;

	//A REVISAR
  // number of acceptors
  AA[p].nacceptors=1;
  //acceptor information
  AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
  (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;

  //Number of Hydrogens connection
  AA[p].nH_hydrogen_connexions=12;
  //atoms hydrogens are connected too
  AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
  AA[p].Hbase[0][0]=9; AA[p].Hbase[0][1]=0;
  AA[p].Hbase[1][0]=10; AA[p].Hbase[1][1]=8;
  AA[p].Hbase[2][0]=11; AA[p].Hbase[2][1]=8;
  AA[p].Hbase[3][0]=12; AA[p].Hbase[3][1]=1;
  AA[p].Hbase[4][0]=13; AA[p].Hbase[4][1]=4;
  AA[p].Hbase[5][0]=14; AA[p].Hbase[5][1]=4;
  AA[p].Hbase[6][0]=15; AA[p].Hbase[6][1]=5;
  AA[p].Hbase[7][0]=16; AA[p].Hbase[7][1]=5;
  AA[p].Hbase[8][0]=17; AA[p].Hbase[8][1]=6;
  AA[p].Hbase[9][0]=18; AA[p].Hbase[9][1]=6;
  AA[p].Hbase[10][0]=19; AA[p].Hbase[10][1]=7;
  AA[p].Hbase[11][0]=20; AA[p].Hbase[11][1]=7;


  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=9;
  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=12;
  AA[p].atom[4].numHydrogens_atm=2;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=13;
  AA[p].atom[4].hydrogens_atm[1]=14;
  AA[p].atom[5].numHydrogens_atm=2;
  AA[p].atom[5].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[5].numHydrogens_atm);
  AA[p].atom[5].hydrogens_atm[0]=15;
  AA[p].atom[5].hydrogens_atm[1]=16;
  AA[p].atom[6].numHydrogens_atm=2;
  AA[p].atom[6].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[6].numHydrogens_atm);
  AA[p].atom[6].hydrogens_atm[0]=17;
  AA[p].atom[6].hydrogens_atm[1]=18;
  AA[p].atom[7].numHydrogens_atm=2;
  AA[p].atom[7].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[7].numHydrogens_atm);
  AA[p].atom[7].hydrogens_atm[0]=19;
  AA[p].atom[7].hydrogens_atm[1]=20;
  AA[p].atom[8].numHydrogens_atm=2;
  AA[p].atom[8].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[8].numHydrogens_atm);
  AA[p].atom[8].hydrogens_atm[0]=10;
  AA[p].atom[8].hydrogens_atm[1]=11;

}

//FIN DE AMINOACIDO  LYSINA NEUTRAL

/**
* Initializes Tyrosine NEGATIVE
*
* @param opt: Rosseta or ICM
*/
void init_TYM(Convention opt)
{
  //AMINOACIDO TIROSINA
  int k, p =TYM;

  strcpy( AA[p].aa_name3, "TYM" );
  AA[p].aa_name1 = 'Y';

  AA[p].mass = 163.178;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = false;
  AA[p].aa_is_nonpolar = true;
  AA[p].aa_is_aromatic = true;
  AA[p].aa_is_charged = false;

  AA[p].natoms = 20; //NEW ADDED
  AA[p].nheavyatoms = 12;
  AA[p].nchi = 2;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );
  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;

  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

 AA[p].chi_types = ( int * ) new int[AA[p].nchi];
 AA[p].chi_types[0] = 2;
 AA[p].chi_types[1] = 4;


  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " CG " );
  strcpy( AA[p].atom[6].atom_name, " CD1" );
  strcpy( AA[p].atom[7].atom_name, " CD2" );
  strcpy( AA[p].atom[8].atom_name, " CE1" );
  strcpy( AA[p].atom[9].atom_name, " CE2" );
  strcpy( AA[p].atom[10].atom_name, " CZ " );
  strcpy( AA[p].atom[11].atom_name, " OH " );
  strcpy( AA[p].atom[12].atom_name, " H  " );
  strcpy( AA[p].atom[13].atom_name, " HD1" );
  strcpy( AA[p].atom[14].atom_name, " HE1" );
  strcpy( AA[p].atom[15].atom_name, " HE2" );
  strcpy( AA[p].atom[16].atom_name, " HD2" );
  strcpy( AA[p].atom[17].atom_name, " HA " );
  strcpy( AA[p].atom[18].atom_name, "2HB " );
  strcpy( AA[p].atom[19].atom_name, "3HB " );

  strcpy( AA[p].atom[20].atom_name," OXT" );
  strcpy( AA[p].atom[21].atom_name,"2H  " );
  strcpy( AA[p].atom[22].atom_name,"3H  " );

  switch(opt)
  {
    //A REVISAR
    case Rosseta:
    case Sybil:
      AA[p].atom[0].charge = -0.47; //    N
      AA[p].atom[1].charge = 0.07; //    CA
      AA[p].atom[2].charge = 0.51; //    C
      AA[p].atom[3].charge = -0.51; //    O
      AA[p].atom[4].charge = -0.18; //    CB
      AA[p].atom[5].charge = 0.00; //    CG
      AA[p].atom[6].charge = -0.115; //    CD1
      AA[p].atom[7].charge = -0.115; //    CD2
      AA[p].atom[8].charge = -0.115; //    CE1
      AA[p].atom[9].charge = -0.115; //    CE2
      AA[p].atom[10].charge = 0.11; //    CZ
      AA[p].atom[11].charge = -0.54; //    OH
      AA[p].atom[12].charge = 0.31; //    H
      AA[p].atom[13].charge = 0.115; //    HD1
      AA[p].atom[14].charge = 0.115; //    HE1
      AA[p].atom[15].charge = 0.115; //    HE2
      AA[p].atom[16].charge = 0.115; //    HD2
      AA[p].atom[17].charge = 0.09; //    HA
      AA[p].atom[18].charge = 0.09; //   2HB
      AA[p].atom[19].charge = 0.09; //   3HB

      AA[p].atom[20].charge = -0.67; //   0XT
      AA[p].atom[21].charge = 0.33; //  1H
      AA[p].atom[22].charge = 0.33; //  2H
    break;
    //A REVISAR
    case ICM:
      AA[p].atom[0].charge = -0.356; //    N
      AA[p].atom[1].charge = 0.064; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = -0.04; //    CB
      AA[p].atom[5].charge = 0.02; //    CG
      AA[p].atom[6].charge = -0.01; //    CD1
      AA[p].atom[7].charge = -0.01; //    CD2
      AA[p].atom[8].charge = -0.06; //    CE1
      AA[p].atom[9].charge = -0.06; //    CE2
      AA[p].atom[10].charge = 0.225; //    CZ
      AA[p].atom[11].charge = -0.33; //    OH
      AA[p].atom[12].charge = 0.176; //    H
      AA[p].atom[13].charge = 0.01; //    HD1
      AA[p].atom[14].charge = 0.03; //    HE1
      AA[p].atom[15].charge = 0.03; //    HE2
      AA[p].atom[16].charge = 0.01; //    HD2
      AA[p].atom[17].charge = 0.02; //    HA
      AA[p].atom[18].charge = 0.025; //   2HB
      AA[p].atom[19].charge = 0.025; //   3HB

      AA[p].atom[20].charge = -0.67; //   0XT
      AA[p].atom[21].charge = 0.33; //  1H
      AA[p].atom[22].charge = 0.33; //  2H

    break;

    default:
      AA[p].atom[0].charge = -0.356; //    N
      AA[p].atom[1].charge = 0.064; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = -0.04; //    CB
      AA[p].atom[5].charge = 0.02; //    CG
      AA[p].atom[6].charge = -0.01; //    CD1
      AA[p].atom[7].charge = -0.01; //    CD2
      AA[p].atom[8].charge = -0.06; //    CE1
      AA[p].atom[9].charge = -0.06; //    CE2
      AA[p].atom[10].charge = 0.225; //    CZ
      AA[p].atom[11].charge = -0.33; //    OH
      AA[p].atom[12].charge = 0.176; //    H
      AA[p].atom[13].charge = 0.01; //    HD1
      AA[p].atom[14].charge = 0.03; //    HE1
      AA[p].atom[15].charge = 0.03; //    HE2
      AA[p].atom[16].charge = 0.01; //    HD2
      AA[p].atom[17].charge = 0.02; //    HA
      AA[p].atom[18].charge = 0.025; //   2HB
      AA[p].atom[19].charge = 0.025; //   3HB

      AA[p].atom[20].charge = -0.67; //   0XT
      AA[p].atom[21].charge = 0.33; //  1H
      AA[p].atom[22].charge = 0.33; //  2H

    break;
  }

  /*
   1.0 TYR-N:2H    0.33000 HC
   1.0 TYR-N:HD2   0.11500 HP
   1.0 TYR-N:HA    0.10000 HB
   1.0 TYR-N:3H    0.33000 HC
   1.0 TYR-N:CE2  -0.11500 CA
   1.0 TYR-N:CB   -0.18000 CT2
   1.0 TYR-N:HE2   0.11500 HP
   1.0 TYR-N:HB1   0.09000 HA
   1.0 TYR-N:1HB   0.09000 HA
   1.0 TYR-N:CG    0.00000 CA
   1.0 TYR-N:CD1  -0.11500 CA
   1.0 TYR-N:HD1   0.11500 HP
   1.0 TYR-N:C     0.51000 C
   1.0 TYR-N:CE1  -0.11500 CA
   1.0 TYR-N:O    -0.51000 O
   1.0 TYR-N:HE1   0.11500 HP
   1.0 TYR-N:CZ    0.11000 CA
   1.0 TYR-N:OH   -0.54000 OH1
   1.0 TYR-N:N    -0.30000 NH3
   1.0 TYR-N:HH    0.43000 H
   1.0 TYR-N:1H    0.33000 HC
   1.0 TYR-N:CD2  -0.11500 CA
   1.0 TYR-N:CA    0.21000 CT1

   1.0 TYR-C:HD2   0.11500 HP
   1.0 TYR-C:HA    0.09000 HB
   1.0 TYR-C:CE2  -0.11500 CA
   1.0 TYR-C:CB   -0.18000 CT2
   1.0 TYR-C:HE2   0.11500 HP
   1.0 TYR-C:HB1   0.09000 HA
   1.0 TYR-C:1HB   0.09000 HA
   1.0 TYR-C:CG    0.00000 CA
   1.0 TYR-C:CD1  -0.11500 CA
   1.0 TYR-C:HD1   0.11500 HP
   1.0 TYR-C:C     0.34000 CC
   1.0 TYR-C:CE1  -0.11500 CA
   1.0 TYR-C:O    -0.67000 OC
   1.0 TYR-C:HE1   0.11500 HP
   1.0 TYR-C:OXT  -0.67000 OC
   1.0 TYR-C:CZ    0.11000 CA
   1.0 TYR-C:OH   -0.54000 OH1
   1.0 TYR-C:N    -0.47000 NH1
   1.0 TYR-C:HH    0.43000 H
   1.0 TYR-C:HN    0.31000 H
   1.0 TYR-C:CD2  -0.11500 CA
   1.0 TYR-C:CA    0.07000 CT1
*/
  switch (opt) {
  case Sybil:
	  AA[p].atom[0].fullatom_type = 40; // Nbb    N
      AA[p].atom[1].fullatom_type = 5; // CAbb   CA
      AA[p].atom[2].fullatom_type = 10; // CObb   C
	  AA[p].atom[3].fullatom_type = 56; // OCbb   O
	  AA[p].atom[4].fullatom_type = 8; // CH2    CB
	  AA[p].atom[5].fullatom_type = 12; // aroC   CG
	  AA[p].atom[6].fullatom_type = 12; // aroC   CD1
	  AA[p].atom[7].fullatom_type = 12; // aroC   CD2
	  AA[p].atom[8].fullatom_type = 12; // aroC   CE1
	  AA[p].atom[9].fullatom_type = 12; // aroC   CE2
	  AA[p].atom[10].fullatom_type = 12; // aroC   CZ
	  AA[p].atom[11].fullatom_type = 55; // OH     OH
 	  AA[p].atom[12].fullatom_type = 66; // HNbb   H
 	  AA[p].atom[13].fullatom_type = 66; // Haro   HD1
 	  AA[p].atom[14].fullatom_type = 66; // Haro   HE1
 	  AA[p].atom[15].fullatom_type = 66; // Haro   HE2
 	  AA[p].atom[16].fullatom_type = 66; // Haro   HD2
 	  AA[p].atom[17].fullatom_type = 66; // Hapo   HA
 	  AA[p].atom[18].fullatom_type = 66; // Hapo  2HB
 	  AA[p].atom[19].fullatom_type = 66; // Hapo  3HB

 	  AA[p].atom[20].fullatom_type = 57; // OXT
 	  AA[p].atom[21].fullatom_type = 66; // Hpol
 	  AA[p].atom[22].fullatom_type = 66; // Hpol
 	  break;
  default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 6; // aroC   CG
	  AA[p].atom[6].fullatom_type = 6; // aroC   CD1
	  AA[p].atom[7].fullatom_type = 6; // aroC   CD2
	  AA[p].atom[8].fullatom_type = 6; // aroC   CE1
	  AA[p].atom[9].fullatom_type = 6; // aroC   CE2
	  AA[p].atom[10].fullatom_type = 6; // aroC   CZ
	  AA[p].atom[11].fullatom_type = 13; // OH     OH
	  AA[p].atom[12].fullatom_type = 25; // HNbb   H
	  AA[p].atom[13].fullatom_type = 24; // Haro   HD1
	  AA[p].atom[14].fullatom_type = 24; // Haro   HE1
	  AA[p].atom[15].fullatom_type = 24; // Haro   HE2
	  AA[p].atom[16].fullatom_type = 24; // Haro   HD2
	  AA[p].atom[17].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[18].fullatom_type = 23; // Hapo  2HB
	  AA[p].atom[19].fullatom_type = 23; // Hapo  3HB

	  AA[p].atom[20].fullatom_type = 20; // OXT
	  AA[p].atom[21].fullatom_type = 22; // Hpol
	  AA[p].atom[22].fullatom_type = 22; // Hpol
	  break;
  }

  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 12; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 17; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--CG
  AA[p].atom[4].bonded_neighbor[2] = 18; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 19; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 3; // CG
  AA[p].atom[5].bonded_neighbor[0] = 4; // CG--CB
  AA[p].atom[5].bonded_neighbor[1] = 6; // CG--CD1
  AA[p].atom[5].bonded_neighbor[2] = 7; // CG--CD2
  AA[p].atom[6].nbonded_neighbors = 3; // CD1
  AA[p].atom[6].bonded_neighbor[0] = 5; // CD1--CG
  AA[p].atom[6].bonded_neighbor[1] = 8; // CD1--CE1
  AA[p].atom[6].bonded_neighbor[2] = 13; // CD1--HD1
  AA[p].atom[7].nbonded_neighbors = 3; // CD2
  AA[p].atom[7].bonded_neighbor[0] = 5; // CD2--CG
  AA[p].atom[7].bonded_neighbor[1] = 9; // CD2--CE2
  AA[p].atom[7].bonded_neighbor[2] = 16; // CD2--HD2
  AA[p].atom[8].nbonded_neighbors = 3; // CE1
  AA[p].atom[8].bonded_neighbor[0] = 6; // CE1--CD1
  AA[p].atom[8].bonded_neighbor[1] = 10; // CE1--CZ
  AA[p].atom[8].bonded_neighbor[2] = 14; // CE1--HE1
  AA[p].atom[9].nbonded_neighbors = 3; // CE2
  AA[p].atom[9].bonded_neighbor[0] = 7; // CE2--CD2
  AA[p].atom[9].bonded_neighbor[1] = 10; // CE2--CZ
  AA[p].atom[9].bonded_neighbor[2] = 15; // CE2--HE2
  AA[p].atom[10].nbonded_neighbors = 3; // CZ
  AA[p].atom[10].bonded_neighbor[0] = 8; // CZ--CE1
  AA[p].atom[10].bonded_neighbor[1] = 9; // CZ--CE2
  AA[p].atom[10].bonded_neighbor[2] = 11; // CZ--OH
  AA[p].atom[11].nbonded_neighbors = 1; // OH
  AA[p].atom[11].bonded_neighbor[0] = 10; // OH--CZ
  AA[p].atom[12].nbonded_neighbors = 1; // H
  AA[p].atom[12].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[13].nbonded_neighbors = 1; // HD1
  AA[p].atom[13].bonded_neighbor[0] = 6; // HD1--CD1
  AA[p].atom[14].nbonded_neighbors = 1; // HE1
  AA[p].atom[14].bonded_neighbor[0] = 8; // HE1--CE1
  AA[p].atom[15].nbonded_neighbors = 1; // HE2
  AA[p].atom[15].bonded_neighbor[0] = 9; // HE2--CE2
  AA[p].atom[16].nbonded_neighbors = 1; // HD2
  AA[p].atom[16].bonded_neighbor[0] = 7; // HD2--CD2
  AA[p].atom[17].nbonded_neighbors = 1; // HA
  AA[p].atom[17].bonded_neighbor[0] = 1; // HA--CA
  AA[p].atom[18].nbonded_neighbors = 1; //2HB
  AA[p].atom[18].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[19].nbonded_neighbors = 1; //3HB
  AA[p].atom[19].bonded_neighbor[0] = 4; //3HB--CB

  AA[p].atom[20].nbonded_neighbors = 1; // OXT
  AA[p].atom[20].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[21].nbonded_neighbors = 1; // 1H
  AA[p].atom[21].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[22].nbonded_neighbors = 1; // 2H
  AA[p].atom[22].bonded_neighbor[0] = 0; // H-N


  //A REVISAR TEMPLATES Y CHIS
  for ( k = 0; k < 13; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 6; k < 20; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template atoms used for placing atoms TYR

  //bk   template for building HD1
  AA[p].atom[13].ta[0] = 6; //   CD1
  AA[p].atom[13].ta[1] = 5; //   CG
  AA[p].atom[13].ta[2] = 8; //   CE1

  //bk   template for building  HE1
  AA[p].atom[14].ta[0] = 8; //   CE1
  AA[p].atom[14].ta[1] = 6; //   CD1
  AA[p].atom[14].ta[2] = 10; //   CZ

  //bk   template for building HE2
  AA[p].atom[15].ta[0] = 9; //   CE2
  AA[p].atom[15].ta[1] = 7; //   CD2
  AA[p].atom[15].ta[2] = 10; //   CZ

  //bk   template for building HD2
  AA[p].atom[16].ta[0] = 7; //   CD2
  AA[p].atom[16].ta[1] = 5; //   CG
  AA[p].atom[16].ta[2] = 9; //   CE2

  //bk   template for building HA
  AA[p].atom[17].ta[0] = 1; //   CA
  AA[p].atom[17].ta[1] = 0; //   N
  AA[p].atom[17].ta[2] = 2; //   C

  //bk   template for building 2HB
  AA[p].atom[18].ta[0] = 4; //   CB
  AA[p].atom[18].ta[1] = 1; //   CA
  AA[p].atom[18].ta[2] = 5; //   CG

  //bk   template for building 3HB
  AA[p].atom[19].ta[0] = 4; //   CB
  AA[p].atom[19].ta[1] = 1; //   CA
  AA[p].atom[19].ta[2] = 5; //   CG


  //bk   chi angles needed for building  CG
  AA[p].chi_required[0] [5] = true;

  //bk   chi angles needed for building  CD1
  AA[p].chi_required[0] [6] = true;
  AA[p].chi_required[1] [6] = true;

  //bk   chi angles needed for building  CD2
  AA[p].chi_required[0] [7] = true;
  AA[p].chi_required[1] [7] = true;

  //bk   chi angles needed for building  CE1
  AA[p].chi_required[0] [8] = true;
  AA[p].chi_required[1] [8] = true;

  //bk   chi angles needed for building  CE2
  AA[p].chi_required[0] [9] = true;
  AA[p].chi_required[1] [9] = true;

  //bk   chi angles needed for building  CZ
  AA[p].chi_required[0] [10] = true;
  AA[p].chi_required[1] [10] = true;

  //bk   chi angles needed for building  OH
  AA[p].chi_required[0] [11] = true;
  AA[p].chi_required[1] [11] = true;


  //bk   chi angles needed for building  HD1
  AA[p].chi_required[0] [13] = true;
  AA[p].chi_required[1] [13] = true;

  //bk   chi angles needed for building  HE1
  AA[p].chi_required[0] [14] = true;
  AA[p].chi_required[1] [14] = true;

  //bk   chi angles needed for building  HE2
  AA[p].chi_required[0] [15] = true;
  AA[p].chi_required[1] [15] = true;

  //bk   chi angles needed for building  HD2
  AA[p].chi_required[0] [16] = true;
  AA[p].chi_required[1] [16] = true;

  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [18] = true;

  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [19] = true;


  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   CG

  //bk   four atoms that define chi angle  2
  AA[p].chi_atoms[1] [0] = 1; //   CA
  AA[p].chi_atoms[1] [1] = 4; //   CB
  AA[p].chi_atoms[1] [2] = 5; //   CG
  AA[p].chi_atoms[1] [3] = 6; //   CD1

  //bk   four atoms that define chi angle  3
  //AA[p].chi_atoms[2] [0] = 9; //   CE2
  //AA[p].chi_atoms[2] [1] = 10; //   CZ
  //AA[p].chi_atoms[2] [2] = 11; //   OH
  //AA[p].chi_atoms[2] [3] = 13; //   HH


  //bk   template coordinates for the amino acid
  //jjh new template from Kevin Karplus

  AA[p].atom[0].icoor[0] = 0.0000; //   N
  AA[p].atom[0].icoor[1] = 17.7910; //   N
  AA[p].atom[0].icoor[2] = 69.1770; //   N
  AA[p].atom[1].icoor[0] = 0.0000; //   CA
  AA[p].atom[1].icoor[1] = 17.0770; //   CA
  AA[p].atom[1].icoor[2] = 70.4470; //   CA
  AA[p].atom[2].icoor[0] = 0.0000; //   C
  AA[p].atom[2].icoor[1] = 18.0440; //   C
  AA[p].atom[2].icoor[2] = 71.6230; //   C
  AA[p].atom[3].icoor[0] = -0.0010; //   O
  AA[p].atom[3].icoor[1] = 19.2610; //   O
  AA[p].atom[3].icoor[2] = 71.4380; //   O
  AA[p].atom[4].icoor[0] = 1.2090; //   CB
  AA[p].atom[4].icoor[1] = 16.1430; //   CB
  AA[p].atom[4].icoor[2] = 70.5360; //   CB
  AA[p].atom[5].icoor[0] = 1.2850; //   CG
  AA[p].atom[5].icoor[1] = 15.3600; //   CG
  AA[p].atom[5].icoor[2] = 71.8280; //   CG
  AA[p].atom[6].icoor[0] = 2.3230; //   CD1
  AA[p].atom[6].icoor[1] = 14.4700; //   CD1
  AA[p].atom[6].icoor[2] = 72.0620; //   CD1
  AA[p].atom[7].icoor[0] = 0.3190; //   CD2
  AA[p].atom[7].icoor[1] = 15.5150; //   CD2
  AA[p].atom[7].icoor[2] = 72.8110; //   CD2
  AA[p].atom[8].icoor[0] = 2.3970; //   CE1
  AA[p].atom[8].icoor[1] = 13.7520; //   CE1
  AA[p].atom[8].icoor[2] = 73.2400; //   CE1
  AA[p].atom[9].icoor[0] = 0.3830; //   CE2
  AA[p].atom[9].icoor[1] = 14.8030; //   CE2
  AA[p].atom[9].icoor[2] = 73.9930; //   CE2
  AA[p].atom[10].icoor[0] = 1.4240; //   CZ
  AA[p].atom[10].icoor[1] = 13.9220; //   CZ
  AA[p].atom[10].icoor[2] = 74.2040; //   CZ
  AA[p].atom[11].icoor[0] = 1.4930; //   OH
  AA[p].atom[11].icoor[1] = 13.2110; //   OH
  AA[p].atom[11].icoor[2] = 75.3800; //   OH
  AA[p].atom[12].icoor[0] = 0.0390; //   H
  AA[p].atom[12].icoor[1] = 18.8010; //   H
  AA[p].atom[12].icoor[2] = 69.1950; //   H
  AA[p].atom[13].icoor[0] = 3.0880; //   HD1
  AA[p].atom[13].icoor[1] = 14.3410; //   HD1
  AA[p].atom[13].icoor[2] = 71.2960; //   HD1
  AA[p].atom[14].icoor[0] = 3.2180; //   HE1
  AA[p].atom[14].icoor[1] = 13.0560; //   HE1
  AA[p].atom[14].icoor[2] = 73.4100; //   HE1
  AA[p].atom[15].icoor[0] = -0.3870; //   HE2
  AA[p].atom[15].icoor[1] = 14.9380; //   HE2
  AA[p].atom[15].icoor[2] = 74.7530; //   HE2
  AA[p].atom[16].icoor[0] = -0.5010; //   HD2
  AA[p].atom[16].icoor[1] = 16.2120; //   HD2
  AA[p].atom[16].icoor[2] = 72.6380; //   HD2
  AA[p].atom[17].icoor[0] = -0.9080; //   HA
  AA[p].atom[17].icoor[1] = 16.4800; //   HA
  AA[p].atom[17].icoor[2] = 70.5360; //   HA
  AA[p].atom[18].icoor[0] = 1.1450; //  2HB
  AA[p].atom[18].icoor[1] = 15.4510; //  2HB
  AA[p].atom[18].icoor[2] = 69.6950; //  2HB
  AA[p].atom[19].icoor[0] = 2.1010; //  3HB
  AA[p].atom[19].icoor[1] = 16.7600; //  3HB
  AA[p].atom[19].icoor[2] = 70.4320; //  3HB

  // atom number for backbone HN
  AA[p].HNpos=12;
  // atom number for backbone HA
  AA[p].HApos=17;

  // number of polar hydrogens
  AA[p].nH_polar_complete=1;
  //atom number for apolar hydrogens
  AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
  AA[p].Hpos_polar_complete[0]=12;


  // number of aromatic hydrogens
  AA[p].nH_aromatic_complete=4;
  AA[p].Hpos_aromatic_complete=(int*)malloc(sizeof(int)*AA[p].nH_aromatic_complete);
  AA[p].Hpos_aromatic_complete[0]=13;
  AA[p].Hpos_aromatic_complete[1]=14;
  AA[p].Hpos_aromatic_complete[2]=15;
  AA[p].Hpos_aromatic_complete[3]=16;

  // number of apolar hydrogens
  AA[p].nH_apolar_complete=3;
  //atom number for apolar hydrogens
  AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
  AA[p].Hpos_apolar_complete[0]=17;
  AA[p].Hpos_apolar_complete[1]=18;
  AA[p].Hpos_apolar_complete[2]=19;

  //A REVISAR
  // number of acceptors
  AA[p].nacceptors=2;
  //acceptor information
  AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
  (AA[p].aBase[1])[0]=11;(AA[p].aBase[1])[1]=13;
  (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;
  (AA[p].accpt_pos[1])[0]=1;(AA[p].accpt_pos[1])[1]=11;

  //Number of Hydrogens connection
  AA[p].nH_hydrogen_connexions=8;
  //atoms hydrogens are connected too
  AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
  AA[p].Hbase[0][0]=12; AA[p].Hbase[0][1]=0;
  AA[p].Hbase[1][0]=13; AA[p].Hbase[1][1]=6;
  AA[p].Hbase[2][0]=14; AA[p].Hbase[2][1]=8;
  AA[p].Hbase[3][0]=15; AA[p].Hbase[3][1]=9;
  AA[p].Hbase[4][0]=16; AA[p].Hbase[4][1]=7;
  AA[p].Hbase[5][0]=17; AA[p].Hbase[5][1]=1;
  AA[p].Hbase[6][0]=18; AA[p].Hbase[6][1]=4;
  AA[p].Hbase[7][0]=19; AA[p].Hbase[7][1]=4;

  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=12;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=17;

  AA[p].atom[4].numHydrogens_atm=2;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=18;
  AA[p].atom[4].hydrogens_atm[1]=19;

  AA[p].atom[6].numHydrogens_atm=1;
  AA[p].atom[6].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[6].numHydrogens_atm);
  AA[p].atom[6].hydrogens_atm[0]=15;

  AA[p].atom[7].numHydrogens_atm=1;
  AA[p].atom[7].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[7].numHydrogens_atm);
  AA[p].atom[7].hydrogens_atm[0]=16;

  AA[p].atom[8].numHydrogens_atm=1;
  AA[p].atom[8].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[8].numHydrogens_atm);
  AA[p].atom[8].hydrogens_atm[0]=14;

  AA[p].atom[9].numHydrogens_atm=1;
  AA[p].atom[9].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[9].numHydrogens_atm);
  AA[p].atom[9].hydrogens_atm[0]=15;


  //FIN DE AMINOACIDO TIROSINA NEGATIVE
}

// Mon added. (3/5/2007)
// Pay attention, ONLY NAMES MODIFIED FROM "MET" ONES !!!
/**
* Initializes SelenoMethionine
*
* @param opt: Rosseta or ICM
*/
void init_MSE(Convention opt)
{
  //AMINOACIDO SELENOMETHIONINE
  int k, p =MSE;

  strcpy( AA[p].aa_name3, "MSE" );
  AA[p].aa_name1 = 'M';

  AA[p].mass = 131.2;

  AA[p].is_protein = true;
  AA[p].is_RNA = false;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = false;
  AA[p].aa_is_nonpolar = true;
  AA[p].aa_is_aromatic = false;
  AA[p].aa_is_charged = false;

  AA[p].natoms = 17;
  AA[p].nheavyatoms = 8;
  AA[p].nchi = 3;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );
  AA[p].chi_required = ( bool * * ) malloc( 4 * sizeof( bool * * ) );
  for ( int j = 0; j < 4; j++ )
    AA[p].chi_required[j] = ( bool * ) malloc( AA[p].natoms * sizeof( bool * ) );

  for ( int j = 0; j < 4; j++ )
    for ( int i = 0; i < AA[p].natoms; i++ )
      AA[p].chi_required[j] [i] = false;

  AA[p].atom = ( t_aa_atom_p * ) malloc( (AA[p].natoms+3) * sizeof( t_aa_atom_p ) );

AA[p].chi_types = ( int * ) new int[AA[p].nchi];
  AA[p].chi_types[0] = 2;
  AA[p].chi_types[1] = 2;
  AA[p].chi_types[2] = 5;


  strcpy( AA[p].atom[0].atom_name, " N  " );
  strcpy( AA[p].atom[1].atom_name, " CA " );
  strcpy( AA[p].atom[2].atom_name, " C  " );
  strcpy( AA[p].atom[3].atom_name, " O  " );
  strcpy( AA[p].atom[4].atom_name, " CB " );
  strcpy( AA[p].atom[5].atom_name, " CG " );
  strcpy( AA[p].atom[6].atom_name, " SD " );
  strcpy( AA[p].atom[7].atom_name, " CE " );
  strcpy( AA[p].atom[8].atom_name, " H  " );
  strcpy( AA[p].atom[9].atom_name, " HA " );
  strcpy( AA[p].atom[10].atom_name, "2HB " );
  strcpy( AA[p].atom[11].atom_name, "3HB " );
  strcpy( AA[p].atom[12].atom_name, "2HG " );
  strcpy( AA[p].atom[13].atom_name, "3HG " );
  strcpy( AA[p].atom[14].atom_name, "1HE " );
  strcpy( AA[p].atom[15].atom_name, "2HE " );
  strcpy( AA[p].atom[16].atom_name, "3HE " );

  strcpy( AA[p].atom[17].atom_name," OXT" );
  strcpy( AA[p].atom[18].atom_name,"2H  " );
  strcpy( AA[p].atom[19].atom_name,"3H  " );


  switch(opt)
  {
    case Rosseta:
    case Sybil:
      AA[p].atom[0].charge = -0.47; //    N
      AA[p].atom[1].charge = 0.07; //    CA
      AA[p].atom[2].charge = 0.51; //    C
      AA[p].atom[3].charge = -0.51; //    O
      AA[p].atom[4].charge = -0.18; //    CB
      AA[p].atom[5].charge = -0.14; //    CG
      AA[p].atom[6].charge = -0.09; //    SD
      AA[p].atom[7].charge = -0.22; //    CE
      AA[p].atom[8].charge = 0.31; //    H
      AA[p].atom[9].charge = 0.09; //    HA
      AA[p].atom[10].charge = 0.09; //   2HB
      AA[p].atom[11].charge = 0.09; //   3HB
      AA[p].atom[12].charge = 0.09; //   2HG
      AA[p].atom[13].charge = 0.09; //   3HG
      AA[p].atom[14].charge = 0.09; //   1HE
      AA[p].atom[15].charge = 0.09; //   2HE
      AA[p].atom[16].charge = 0.09; //   3HE

      AA[p].atom[17].charge = -0.67; //   0XT
      AA[p].atom[18].charge = 0.33; //  1H
      AA[p].atom[19].charge = 0.33; //  2H
    break;

    case ICM:
      AA[p].atom[0].charge = -0.356; //    N
      AA[p].atom[1].charge = 0.064; //    CA
      AA[p].atom[2].charge = 0.45; //    C
      AA[p].atom[3].charge = -0.384; //    O
      AA[p].atom[4].charge = -0.01; //    CB
      AA[p].atom[5].charge = -0.12; //    CG
      AA[p].atom[6].charge = -0.035; //    SD
      AA[p].atom[7].charge = -0.19; //    CE
      AA[p].atom[8].charge = 0.176; //    H
      AA[p].atom[9].charge = 0.02; //    HA
      AA[p].atom[10].charge = 0.03; //   2HB
      AA[p].atom[11].charge = 0.03; //   3HB
      AA[p].atom[12].charge = 0.045; //   2HG
      AA[p].atom[13].charge = 0.045; //   3HG
      AA[p].atom[14].charge = 0.055; //   1HE
      AA[p].atom[15].charge = 0.055; //   2HE
      AA[p].atom[16].charge = 0.055; //   3HE

      AA[p].atom[17].charge = -0.67; //   0XT
      AA[p].atom[18].charge = 0.33; //  1H
      AA[p].atom[19].charge = 0.33; //  2H
      break;
	default:
	  break;

  }
  /*
   1.0 MET-N:2H    0.33000 HC
   1.0 MET-N:O    -0.51000 O
   1.0 MET-N:HA    0.10000 HB
   1.0 MET-N:3H    0.33000 HC
   1.0 MET-N:CB   -0.18000 CT2
   1.0 MET-N:HB1   0.09000 HA
   1.0 MET-N:1HB   0.09000 HA
   1.0 MET-N:CG   -0.14000 CT2
   1.0 MET-N:HG1   0.09000 HA
   1.0 MET-N:1HG   0.09000 HA
   1.0 MET-N:SD   -0.09000 S
   1.0 MET-N:CE   -0.22000 CT3
   1.0 MET-N:1HE   0.09000 HA
   1.0 MET-N:2HE   0.09000 HA
   1.0 MET-N:N    -0.30000 NH3
   1.0 MET-N:3HE   0.09000 HA
   1.0 MET-N:1H    0.33000 HC
   1.0 MET-N:C     0.51000 C
   1.0 MET-N:CA    0.21000 CT1

   1.0 MET-C:HA    0.09000 HB
   1.0 MET-C:CB   -0.18000 CT2
   1.0 MET-C:HB1   0.09000 HA
   1.0 MET-C:1HB   0.09000 HA
   1.0 MET-C:CG   -0.14000 CT2
   1.0 MET-C:HG1   0.09000 HA
   1.0 MET-C:1HG   0.09000 HA
   1.0 MET-C:SD   -0.09000 S
   1.0 MET-C:O    -0.67000 OC
   1.0 MET-C:CE   -0.22000 CT3
   1.0 MET-C:OXT  -0.67000 OC
   1.0 MET-C:1HE   0.09000 HA
   1.0 MET-C:2HE   0.09000 HA
   1.0 MET-C:N    -0.47000 NH1
   1.0 MET-C:3HE   0.09000 HA
   1.0 MET-C:HN    0.31000 H
   1.0 MET-C:C     0.34000 CC
   1.0 MET-C:CA    0.07000 CT1
*/
  switch (opt) {
  case Sybil:
 	  AA[p].atom[0].fullatom_type = 40; // Nbb    N
 	  AA[p].atom[1].fullatom_type = 5; // CAbb   CA
 	  AA[p].atom[2].fullatom_type = 10; // CObb   C
 	  AA[p].atom[3].fullatom_type = 56; // OCbb   O
 	  AA[p].atom[4].fullatom_type = 8; // CH2    CB
 	  AA[p].atom[5].fullatom_type = 8; // CH2    CG
 	  AA[p].atom[6].fullatom_type = 63; // Se??      SD
 	  AA[p].atom[7].fullatom_type = 6; // CH3    CE
 	  AA[p].atom[8].fullatom_type = 66; // HNbb   H
 	  AA[p].atom[9].fullatom_type = 66; // Hapo   HA
 	  AA[p].atom[10].fullatom_type = 66; // Hapo  2HB
 	  AA[p].atom[11].fullatom_type = 66; // Hapo  3HB
 	  AA[p].atom[12].fullatom_type = 66; // Hapo  2HG
 	  AA[p].atom[13].fullatom_type = 66; // Hapo  3HG
 	  AA[p].atom[14].fullatom_type = 66; // Hapo  1HE
 	  AA[p].atom[15].fullatom_type = 66; // Hapo  2HE
 	  AA[p].atom[16].fullatom_type = 66; // Hapo  3HE

 	  AA[p].atom[17].fullatom_type = 57; // OXT
 	  AA[p].atom[18].fullatom_type = 66; // Hpol
 	  AA[p].atom[19].fullatom_type = 66; // Hpol
 	  break;
  default:
	  AA[p].atom[0].fullatom_type = 17; // Nbb    N
	  AA[p].atom[1].fullatom_type = 18; // CAbb   CA
	  AA[p].atom[2].fullatom_type = 19; // CObb   C
	  AA[p].atom[3].fullatom_type = 20; // OCbb   O
	  AA[p].atom[4].fullatom_type = 4; // CH2    CB
	  AA[p].atom[5].fullatom_type = 4; // CH2    CG
	  AA[p].atom[6].fullatom_type = 16; // S      SD
	  AA[p].atom[7].fullatom_type = 5; // CH3    CE
	  AA[p].atom[8].fullatom_type = 25; // HNbb   H
	  AA[p].atom[9].fullatom_type = 23; // Hapo   HA
	  AA[p].atom[10].fullatom_type = 23; // Hapo  2HB
	  AA[p].atom[11].fullatom_type = 23; // Hapo  3HB
	  AA[p].atom[12].fullatom_type = 23; // Hapo  2HG
	  AA[p].atom[13].fullatom_type = 23; // Hapo  3HG
	  AA[p].atom[14].fullatom_type = 23; // Hapo  1HE
	  AA[p].atom[15].fullatom_type = 23; // Hapo  2HE
	  AA[p].atom[16].fullatom_type = 23; // Hapo  3HE

	  AA[p].atom[17].fullatom_type = 20; // OXT
	  AA[p].atom[18].fullatom_type = 22; // Hpol
	  AA[p].atom[19].fullatom_type = 22; // Hpol
	  break;
  }

  AA[p].atom[0].nbonded_neighbors = 2; // N
  AA[p].atom[0].bonded_neighbor[0] = 1; // N--CA
  AA[p].atom[0].bonded_neighbor[1] = 8; // N--H
  AA[p].atom[1].nbonded_neighbors = 4; // CA
  AA[p].atom[1].bonded_neighbor[0] = 0; // CA--N
  AA[p].atom[1].bonded_neighbor[1] = 2; // CA--C
  AA[p].atom[1].bonded_neighbor[2] = 4; // CA--CB
  AA[p].atom[1].bonded_neighbor[3] = 9; // CA--HA
  AA[p].atom[2].nbonded_neighbors = 2; // C
  AA[p].atom[2].bonded_neighbor[0] = 1; // C--CA
  AA[p].atom[2].bonded_neighbor[1] = 3; // C--O
  AA[p].atom[3].nbonded_neighbors = 1; // O
  AA[p].atom[3].bonded_neighbor[0] = 2; // O--C
  AA[p].atom[4].nbonded_neighbors = 4; // CB
  AA[p].atom[4].bonded_neighbor[0] = 1; // CB--CA
  AA[p].atom[4].bonded_neighbor[1] = 5; // CB--CG
  AA[p].atom[4].bonded_neighbor[2] = 10; // CB--2HB
  AA[p].atom[4].bonded_neighbor[3] = 11; // CB--3HB
  AA[p].atom[5].nbonded_neighbors = 4; // CG
  AA[p].atom[5].bonded_neighbor[0] = 4; // CG--CB
  AA[p].atom[5].bonded_neighbor[1] = 6; // CG--CD
  AA[p].atom[5].bonded_neighbor[2] = 12; // CG--2HG
  AA[p].atom[5].bonded_neighbor[3] = 13; // CG--3HG
  AA[p].atom[6].nbonded_neighbors = 2; // SD
  AA[p].atom[6].bonded_neighbor[0] = 5; // SD--CG
  AA[p].atom[6].bonded_neighbor[1] = 7; // SD--CE
  AA[p].atom[7].nbonded_neighbors = 4; // CE
  AA[p].atom[7].bonded_neighbor[0] = 6; // CE--SD
  AA[p].atom[7].bonded_neighbor[1] = 14; // CE--1HE
  AA[p].atom[7].bonded_neighbor[2] = 15; // CE--2HE
  AA[p].atom[7].bonded_neighbor[3] = 16; // CE--3HE
  AA[p].atom[8].nbonded_neighbors = 1; // H
  AA[p].atom[8].bonded_neighbor[0] = 0; // H--N
  AA[p].atom[9].nbonded_neighbors = 1; //HA
  AA[p].atom[9].bonded_neighbor[0] = 1; //HA--CA
  AA[p].atom[10].nbonded_neighbors = 1; //2HB
  AA[p].atom[10].bonded_neighbor[0] = 4; //2HB--CB
  AA[p].atom[11].nbonded_neighbors = 1; //3HB
  AA[p].atom[11].bonded_neighbor[0] = 4; //3HB--CB
  AA[p].atom[12].nbonded_neighbors = 1; // 2HG
  AA[p].atom[12].bonded_neighbor[0] = 5; //2HG--CG
  AA[p].atom[13].nbonded_neighbors = 1; //3HG
  AA[p].atom[13].bonded_neighbor[0] = 5; //3HG--CG
  AA[p].atom[14].nbonded_neighbors = 1; //1HE
  AA[p].atom[14].bonded_neighbor[0] = 7; //1HE--CE
  AA[p].atom[15].nbonded_neighbors = 1; //2HE
  AA[p].atom[15].bonded_neighbor[0] = 7; //2HE--CE
  AA[p].atom[16].nbonded_neighbors = 1; //3HE
  AA[p].atom[16].bonded_neighbor[0] = 7; //3HE--CG

  AA[p].atom[17].nbonded_neighbors = 1; // OXT
  AA[p].atom[17].bonded_neighbor[0] = 2; // OXT-C
  AA[p].atom[18].nbonded_neighbors = 1; // 1H
  AA[p].atom[18].bonded_neighbor[0] = 0; // H-N
  AA[p].atom[19].nbonded_neighbors = 1; // 2H
  AA[p].atom[19].bonded_neighbor[0] = 0; // H-N


  for ( k = 0; k < 9; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 9; k < 17; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template for building  HA
  AA[p].atom[9].ta[0] = 1; //   CA
  AA[p].atom[9].ta[1] = 0; //   N
  AA[p].atom[9].ta[2] = 2; //   C

  //bk   template for building   2HB
  AA[p].atom[10].ta[0] = 4; //   CB
  AA[p].atom[10].ta[1] = 1; //   CA
  AA[p].atom[10].ta[2] = 5; //   CG

  //bk   template for building  3HB
  AA[p].atom[11].ta[0] = 4; //   CB
  AA[p].atom[11].ta[1] = 1; //   CA
  AA[p].atom[11].ta[2] = 5; //   CG

  //bk   template for building  2HG
  AA[p].atom[12].ta[0] = 5; //   CG
  AA[p].atom[12].ta[1] = 4; //   CB
  AA[p].atom[12].ta[2] = 6; //   SD

  //bk   template for building  3HG
  AA[p].atom[13].ta[0] = 5; //   CG
  AA[p].atom[13].ta[1] = 4; //   CB
  AA[p].atom[13].ta[2] = 6; //   SD

  //bk   template for building  1HE
  AA[p].atom[14].ta[0] = 7; //   CE
  AA[p].atom[14].ta[1] = 6; //   SD
  AA[p].atom[14].ta[2] = 5; //   CG

  //bk   template for building  2HE
  AA[p].atom[15].ta[0] = 7; //   CE
  AA[p].atom[15].ta[1] = 6; //   SD
  AA[p].atom[15].ta[2] = 5; //   CG

  //bk   template for building  3HE
  AA[p].atom[16].ta[0] = 7; //   CE
  AA[p].atom[16].ta[1] = 6; //   SD
  AA[p].atom[16].ta[2] = 5; //   CG
  //bk   chi angles required to build atoms MET

  //bk   chi angles needed for building  CG
  AA[p].chi_required[0] [5] = true;

  //bk   chi angles needed for building  SD
  AA[p].chi_required[0] [6] = true;
  AA[p].chi_required[1] [6] = true;


  //bk   chi angles needed for building  CE
  AA[p].chi_required[0] [7] = true;
  AA[p].chi_required[1] [7] = true;
  AA[p].chi_required[2] [7] = true;

  //bk   chi angles needed for building 2HB
  AA[p].chi_required[0] [10] = true;

  //bk   chi angles needed for building 3HB
  AA[p].chi_required[0] [11] = true;


  //bk   chi angles needed for building 2HG
  AA[p].chi_required[0] [12] = true;
  AA[p].chi_required[1] [12] = true;


  //bk   chi angles needed for building 3HG
  AA[p].chi_required[0] [13] = true;
  AA[p].chi_required[1] [13] = true;


  //bk   chi angles needed for building 1HE
  AA[p].chi_required[0] [14] = true;
  AA[p].chi_required[1] [14] = true;
  AA[p].chi_required[2] [14] = true;


  //bk   chi angles needed for building 2HE
  AA[p].chi_required[0] [15] = true;
  AA[p].chi_required[1] [15] = true;
  AA[p].chi_required[2] [15] = true;


  //bk   chi angles needed for building 3HE
  AA[p].chi_required[0] [16] = true;
  AA[p].chi_required[1] [16] = true;
  AA[p].chi_required[2] [16] = true;





  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 0; //   N
  AA[p].chi_atoms[0] [1] = 1; //   CA
  AA[p].chi_atoms[0] [2] = 4; //   CB
  AA[p].chi_atoms[0] [3] = 5; //   CG
  //bk   four atoms that define chi angle  2
  AA[p].chi_atoms[1] [0] = 1; //   CA
  AA[p].chi_atoms[1] [1] = 4; //   CB
  AA[p].chi_atoms[1] [2] = 5; //   CG
  AA[p].chi_atoms[1] [3] = 6; //   SD
  //bk   four atoms that define chi angle  3
  AA[p].chi_atoms[2] [0] = 4; //   CB
  AA[p].chi_atoms[2] [1] = 5; //   CG
  AA[p].chi_atoms[2] [2] = 6; //   SD
  AA[p].chi_atoms[2] [3] = 7; //   CE



  AA[p].atom[0].icoor[0] = 0.0000; //   N
  AA[p].atom[0].icoor[1] = 15.5330; //   N
  AA[p].atom[0].icoor[2] = 36.4900; //   N
  AA[p].atom[1].icoor[0] = 0.0000; //   CA
  AA[p].atom[1].icoor[1] = 16.4650; //   CA
  AA[p].atom[1].icoor[2] = 37.6220; //   CA
  AA[p].atom[2].icoor[0] = 0.0000; //   C
  AA[p].atom[2].icoor[1] = 15.5850; //   C
  AA[p].atom[2].icoor[2] = 38.8610; //   C
  AA[p].atom[3].icoor[0] = 0.0630; //   O
  AA[p].atom[3].icoor[1] = 14.3610; //   O
  AA[p].atom[3].icoor[2] = 38.7570; //   O
  AA[p].atom[4].icoor[0] = -1.2080; //   CB
  AA[p].atom[4].icoor[1] = 17.3960; //   CB
  AA[p].atom[4].icoor[2] = 37.5390; //   CB
  AA[p].atom[5].icoor[0] = -1.3030; //   CG
  AA[p].atom[5].icoor[1] = 18.4070; //   CG
  AA[p].atom[5].icoor[2] = 38.6730; //   CG
  AA[p].atom[6].icoor[0] = -2.7450; //   SD
  AA[p].atom[6].icoor[1] = 19.4810; //   SD
  AA[p].atom[6].icoor[2] = 38.5280; //   SD
  AA[p].atom[7].icoor[0] = -2.5470; //   CE
  AA[p].atom[7].icoor[1] = 20.5110; //   CE
  AA[p].atom[7].icoor[2] = 39.9790; //   CE
  AA[p].atom[8].icoor[0] = -0.0400; //   H
  AA[p].atom[8].icoor[1] = 14.5430; //   H
  AA[p].atom[8].icoor[2] = 36.6850; //   H
  AA[p].atom[9].icoor[0] = 0.9070; //   HA
  AA[p].atom[9].icoor[1] = 17.0690; //   HA
  AA[p].atom[9].icoor[2] = 37.6060; //   HA
  AA[p].atom[10].icoor[0] = -1.1380; //   2HB
  AA[p].atom[10].icoor[1] = 17.9220; //   2HB
  AA[p].atom[10].icoor[2] = 36.5880; //   2HB
  AA[p].atom[11].icoor[0] = -2.0970; //   3HB
  AA[p].atom[11].icoor[1] = 16.7640; //   3HB
  AA[p].atom[11].icoor[2] = 37.5410; //   3HB
  AA[p].atom[12].icoor[0] = -1.3580; //   2HG
  AA[p].atom[12].icoor[1] = 17.8590; //   2HG
  AA[p].atom[12].icoor[2] = 39.6130; //   2HG
  AA[p].atom[13].icoor[0] = -0.4000; //   3HG
  AA[p].atom[13].icoor[1] = 19.0170; //   3HG
  AA[p].atom[13].icoor[2] = 38.6600; //   3HG
  AA[p].atom[14].icoor[0] = -3.3650; //   1HE
  AA[p].atom[14].icoor[1] = 21.2300; //   1HE
  AA[p].atom[14].icoor[2] = 40.0300; //   1HE
  AA[p].atom[15].icoor[0] = -2.5560; //   2HE
  AA[p].atom[15].icoor[1] = 19.8860; //   2HE
  AA[p].atom[15].icoor[2] = 40.8730; //   2HE
  AA[p].atom[16].icoor[0] = -1.5980; //   3HE
  AA[p].atom[16].icoor[1] = 21.0450; //   3HE
  AA[p].atom[16].icoor[2] = 39.9190; //   3HE

  // atom number for backbone HN
  AA[p].HNpos=8;
  // atom number for backbone HA
  AA[p].HApos=9;

  // number of polar hydrogens
  AA[p].nH_polar_complete=1;
  // atom numbers for polar H
  AA[p].Hpos_polar_complete=(int*)malloc(sizeof(int)*AA[p].nH_polar_complete);
  AA[p].Hpos_polar_complete[0]=8;

  // number of aromatic hydrogens
  AA[p].nH_aromatic_complete=0;

  // number of apolar hydrogens
  AA[p].nH_apolar_complete=8;
  //atom number for apolar hydrogens
  AA[p].Hpos_apolar_complete=(int*)malloc(sizeof(int)*AA[p].nH_apolar_complete);
  AA[p].Hpos_apolar_complete[0]=9;
  AA[p].Hpos_apolar_complete[1]=10;
  AA[p].Hpos_apolar_complete[2]=11;
  AA[p].Hpos_apolar_complete[3]=12;
  AA[p].Hpos_apolar_complete[4]=13;
  AA[p].Hpos_apolar_complete[5]=14;
  AA[p].Hpos_apolar_complete[6]=15;
  AA[p].Hpos_apolar_complete[7]=16;

  // number of acceptors
  AA[p].nacceptors=1;
  //acceptor information
  AA[p].aBase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  AA[p].accpt_pos=(Tpair*)malloc(sizeof(Tpair)*AA[p].nacceptors);
  (AA[p].aBase[0])[0]=3;(AA[p].aBase[0])[1]=2;
  (AA[p].accpt_pos[0])[0]=0;(AA[p].accpt_pos[0])[1]=3;

  //Number of Hydrogens connection
  AA[p].nH_hydrogen_connexions=9;
  //atoms hydrogens are connected too
  AA[p].Hbase=(Tpair*)malloc(sizeof(Tpair)*AA[p].nH_hydrogen_connexions);
  AA[p].Hbase[0][0]=8; AA[p].Hbase[0][1]=0;
  AA[p].Hbase[1][0]=9; AA[p].Hbase[1][1]=1;
  AA[p].Hbase[2][0]=10; AA[p].Hbase[2][1]=4;
  AA[p].Hbase[3][0]=11; AA[p].Hbase[3][1]=4;
  AA[p].Hbase[4][0]=12; AA[p].Hbase[4][1]=5;
  AA[p].Hbase[5][0]=13; AA[p].Hbase[5][1]=5;
  AA[p].Hbase[6][0]=14; AA[p].Hbase[6][1]=7;
  AA[p].Hbase[7][0]=15; AA[p].Hbase[7][1]=7;
  AA[p].Hbase[8][0]=16; AA[p].Hbase[8][1]=7;

  AA[p].atom[0].numHydrogens_atm=1;
  AA[p].atom[0].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[0].numHydrogens_atm);
  AA[p].atom[0].hydrogens_atm[0]=8;

  AA[p].atom[1].numHydrogens_atm=1;
  AA[p].atom[1].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[1].numHydrogens_atm);
  AA[p].atom[1].hydrogens_atm[0]=9;

  AA[p].atom[4].numHydrogens_atm=2;
  AA[p].atom[4].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[4].numHydrogens_atm);
  AA[p].atom[4].hydrogens_atm[0]=10;
  AA[p].atom[4].hydrogens_atm[1]=11;

  AA[p].atom[7].numHydrogens_atm=3;
  AA[p].atom[7].hydrogens_atm=(int*)malloc(sizeof(int)*AA[p].atom[7].numHydrogens_atm);
  AA[p].atom[7].hydrogens_atm[0]=14;
  AA[p].atom[7].hydrogens_atm[1]=15;
  AA[p].atom[7].hydrogens_atm[2]=16;


  //FIN DE AMINOACIDO SELENOMETHIONINE
}


void aa_iupac2pqr() {
//ALA
  strcpy( AA[ALA].atom[7].atom_name, " HB1" );
  strcpy( AA[ALA].atom[8].atom_name, " HB2" );
  strcpy( AA[ALA].atom[9].atom_name, " HB3" );
  strcpy( AA[ALA].atom[11].atom_name," H2 " );
  strcpy( AA[ALA].atom[12].atom_name," H3 " );
//CYS
  strcpy( AA[CYS].atom[8].atom_name, " HB2" );
  strcpy( AA[CYS].atom[9].atom_name, " HB3" );
  strcpy( AA[CYS].atom[12].atom_name," H2 " );
  strcpy( AA[CYS].atom[13].atom_name," H3 " );
//CYX
  strcpy( AA[CYX].atom[8].atom_name, " HB2" );
  strcpy( AA[CYX].atom[9].atom_name, " HB3" );
  strcpy( AA[CYX].atom[11].atom_name," H2 " );
  strcpy( AA[CYX].atom[12].atom_name," H3 " );
//CYM
  strcpy( AA[CYM].atom[8].atom_name, " HB2" );
  strcpy( AA[CYM].atom[9].atom_name, " HB3" );
  strcpy( AA[CYM].atom[11].atom_name," H2 " );
  strcpy( AA[CYM].atom[12].atom_name," H3 " );
//ASP
  strcpy( AA[ASP].atom[10].atom_name," HB2" );
  strcpy( AA[ASP].atom[11].atom_name," HB3" );
  strcpy( AA[ASP].atom[13].atom_name," H2 " );
  strcpy( AA[ASP].atom[14].atom_name," H3 " );
//ASH
  strcpy( AA[ASH].atom[10].atom_name, " HB2" );
  strcpy( AA[ASH].atom[11].atom_name, " HB3" );
  strcpy( AA[ASH].atom[14].atom_name," H2 " );
  strcpy( AA[ASH].atom[15].atom_name," H3 " );
  strcpy( AA[ASH].atom[12].atom_name, " HD2" ); //NEW ATOM
//GLU
  strcpy( AA[GLU].atom[11].atom_name, " HB2" );
  strcpy( AA[GLU].atom[12].atom_name, " HB3" );
  strcpy( AA[GLU].atom[13].atom_name, " HG2" );
  strcpy( AA[GLU].atom[14].atom_name, " HG3" );
  strcpy( AA[GLU].atom[16].atom_name," H2 " );
  strcpy( AA[GLU].atom[17].atom_name," H3 " );

//PHE
  strcpy( AA[PHE].atom[18].atom_name, " HB2" );
  strcpy( AA[PHE].atom[19].atom_name, " HB3" );
  strcpy( AA[PHE].atom[21].atom_name, " H2 " );
  strcpy( AA[PHE].atom[22].atom_name, " H3 " );
//GLY
  strcpy( AA[GLY].atom[5].atom_name, " HA2" );
  strcpy( AA[GLY].atom[6].atom_name, " HA3" );
   strcpy( AA[GLY].atom[8].atom_name," H2 " );
  strcpy( AA[GLY].atom[9].atom_name, " H3 " );
//GLH
 strcpy( AA[GLH].atom[17].atom_name,  " H2 " );
  strcpy( AA[GLH].atom[18].atom_name, " H3 " );
  strcpy( AA[GLH].atom[15].atom_name, " HE2" );//NEW ADDED
//HIS
  strcpy( AA[HIS].atom[13].atom_name," HB2" );
  strcpy( AA[HIS].atom[14].atom_name," HB3" );
   strcpy( AA[HIS].atom[18].atom_name," H2 " );
  strcpy( AA[HIS].atom[19].atom_name," H3 " );
//HIP
  strcpy( AA[HIP].atom[11].atom_name," HE2" ); // <--- The proton.
  strcpy( AA[HIP].atom[12].atom_name," HD1" ); // <--- The proton.
  strcpy( AA[HIP].atom[13].atom_name," HA " );
  strcpy( AA[HIP].atom[14].atom_name," HB2" );
  strcpy( AA[HIP].atom[15].atom_name," HB3" );
  strcpy( AA[HIP].atom[16].atom_name," HE1" );
  strcpy( AA[HIP].atom[17].atom_name," HD2" );
  strcpy( AA[HIS].atom[19].atom_name," H2 " );
  strcpy( AA[HIS].atom[20].atom_name," H3 " );
//HID
  strcpy( AA[HID].atom[13].atom_name," HB2" );
  strcpy( AA[HID].atom[14].atom_name," HB3" );
 strcpy( AA[HID].atom[18].atom_name," H2 " );
  strcpy( AA[HID].atom[19].atom_name," H3 " );
  strcpy( AA[HID].atom[11].atom_name," HD1" );//NEW ADDED
//HIE
  strcpy( AA[HIE].atom[13].atom_name," HB2" );
  strcpy( AA[HIE].atom[14].atom_name," HB3" );
	strcpy( AA[HIE].atom[18].atom_name," H2 " );
  strcpy( AA[HIE].atom[19].atom_name," H3 " );
//ILE
  strcpy( AA[ILE].atom[11].atom_name, "HG21" );
  strcpy( AA[ILE].atom[12].atom_name, "HG22" );
  strcpy( AA[ILE].atom[13].atom_name, "HG23" );
  strcpy( AA[ILE].atom[14].atom_name, "HG12" );
  strcpy( AA[ILE].atom[15].atom_name, "HG13" );
  strcpy( AA[ILE].atom[16].atom_name, "HD11" );
  strcpy( AA[ILE].atom[17].atom_name, "HD12" );
  strcpy( AA[ILE].atom[18].atom_name, "HD13" );
  strcpy( AA[ILE].atom[20].atom_name," H2 " );
  strcpy( AA[ILE].atom[21].atom_name," H3 " );
//LYS


  strcpy( AA[LYS].atom[10].atom_name, " HZ1" );
  strcpy( AA[LYS].atom[11].atom_name, " HZ2" );
  strcpy( AA[LYS].atom[12].atom_name, " HZ3" );
  strcpy( AA[LYS].atom[13].atom_name, " HA " );
  strcpy( AA[LYS].atom[14].atom_name, " HB2" );
  strcpy( AA[LYS].atom[15].atom_name, " HB3" );
  strcpy( AA[LYS].atom[16].atom_name, " HG2" );
  strcpy( AA[LYS].atom[17].atom_name, " HG3" );
  strcpy( AA[LYS].atom[18].atom_name, " HD2" );
  strcpy( AA[LYS].atom[19].atom_name, " HD3" );
  strcpy( AA[LYS].atom[20].atom_name, " HE2" );
  strcpy( AA[LYS].atom[21].atom_name, " HE3" );
  strcpy( AA[LYS].atom[23].atom_name, " H2  " );
  strcpy( AA[LYS].atom[24].atom_name, " H3 " );

//LYN

  strcpy( AA[LYN].atom[10].atom_name, " HZ1" );
  strcpy( AA[LYN].atom[11].atom_name, " HZ2 " );
  strcpy( AA[LYN].atom[12].atom_name, " HA " );
  strcpy( AA[LYN].atom[13].atom_name, " HB2" );
  strcpy( AA[LYN].atom[14].atom_name, " HB3" );
  strcpy( AA[LYN].atom[15].atom_name, " HG2" );
  strcpy( AA[LYN].atom[16].atom_name, " HG3" );
  strcpy( AA[LYN].atom[17].atom_name, " HD2" );
  strcpy( AA[LYN].atom[18].atom_name, " HD3" );
  strcpy( AA[LYN].atom[19].atom_name, " HE2" );
  strcpy( AA[LYN].atom[20].atom_name, " HE3 " );
  strcpy( AA[LYN].atom[22].atom_name, " H2 " );
  strcpy( AA[LYN].atom[23].atom_name, " H3 " );


//LEU
  strcpy( AA[LEU].atom[10].atom_name, " HB2" );
  strcpy( AA[LEU].atom[11].atom_name, " HB3" );
  strcpy( AA[LEU].atom[13].atom_name, "HD11" );
  strcpy( AA[LEU].atom[14].atom_name, "HD12" );
  strcpy( AA[LEU].atom[15].atom_name, "HD13" );
  strcpy( AA[LEU].atom[16].atom_name, "HD21" );
  strcpy( AA[LEU].atom[17].atom_name, "HD22" );
  strcpy( AA[LEU].atom[18].atom_name, "HD23" );
  strcpy( AA[LEU].atom[20].atom_name," H2 " );
  strcpy( AA[LEU].atom[21].atom_name," H3 " );
//MET
 	strcpy( AA[MET].atom[10].atom_name, " HB2" );
  strcpy( AA[MET].atom[11].atom_name, " HB3" );
  strcpy( AA[MET].atom[12].atom_name, " HG2" );
  strcpy( AA[MET].atom[13].atom_name, " HG3" );
  strcpy( AA[MET].atom[14].atom_name, " HE1" );
  strcpy( AA[MET].atom[15].atom_name, " HE2" );
  strcpy( AA[MET].atom[16].atom_name, " HE3" );
  strcpy( AA[MET].atom[18].atom_name," H2 " );
  strcpy( AA[MET].atom[19].atom_name," H3 " );
//ASN
  strcpy( AA[ASN].atom[9].atom_name,  "HD21" );
  strcpy( AA[ASN].atom[10].atom_name, "HD22" );
  strcpy( AA[ASN].atom[12].atom_name, " HB2" );
  strcpy( AA[ASN].atom[13].atom_name, " HB3" );
  strcpy( AA[ASN].atom[15].atom_name," H2 " );
  strcpy( AA[ASN].atom[16].atom_name," H3 " );
//PRO
 strcpy( AA[PRO].atom[7].atom_name, " HD2" );
  strcpy( AA[PRO].atom[8].atom_name, " HD3" );
  strcpy( AA[PRO].atom[9].atom_name, " HG2" );
  strcpy( AA[PRO].atom[10].atom_name, " HG3" );
  strcpy( AA[PRO].atom[11].atom_name, " HB2" );
  strcpy( AA[PRO].atom[12].atom_name, " HB3" );
   strcpy( AA[PRO].atom[15].atom_name," H2 " );
  strcpy( AA[PRO].atom[16].atom_name," H  " );
//GLN
  strcpy( AA[GLN].atom[10].atom_name, "HE21" );
  strcpy( AA[GLN].atom[11].atom_name, "HE22" );
  strcpy( AA[GLN].atom[12].atom_name, " HA " );
  strcpy( AA[GLN].atom[13].atom_name, " HB2" );
  strcpy( AA[GLN].atom[14].atom_name, " HB3" );
  strcpy( AA[GLN].atom[15].atom_name, " HG2" );
  strcpy( AA[GLN].atom[16].atom_name, " HG3" );
  strcpy( AA[GLN].atom[18].atom_name," H2 " );
  strcpy( AA[GLN].atom[19].atom_name," H3 " );
//ARG
 strcpy( AA[ARG].atom[12].atom_name, "HH11" );
  strcpy( AA[ARG].atom[13].atom_name, "HH12" );
  strcpy( AA[ARG].atom[14].atom_name, "HH21" );
  strcpy( AA[ARG].atom[15].atom_name, "HH22" );
  strcpy( AA[ARG].atom[18].atom_name, " HB2" );
  strcpy( AA[ARG].atom[19].atom_name, " HB3" );
  strcpy( AA[ARG].atom[20].atom_name, " HG2" );
  strcpy( AA[ARG].atom[21].atom_name, " HG3" );
  strcpy( AA[ARG].atom[22].atom_name, " HD2" );
  strcpy( AA[ARG].atom[23].atom_name, " HD3" );
   strcpy( AA[ARG].atom[25].atom_name," H2 " );
  strcpy( AA[ARG].atom[26].atom_name," H3 " );
//SER
   strcpy( AA[SER].atom[9].atom_name,  " HB2" );
  strcpy( AA[SER].atom[10].atom_name, " HB3" );
   strcpy( AA[SER].atom[12].atom_name," H2 " );
  strcpy( AA[SER].atom[13].atom_name," H3 " );
//THR
  strcpy( AA[THR].atom[11].atom_name,"HG21" );
  strcpy( AA[THR].atom[12].atom_name,"HG22" );
  strcpy( AA[THR].atom[13].atom_name,"HG23" );
    strcpy( AA[THR].atom[15].atom_name," H2 " );
  strcpy( AA[THR].atom[16].atom_name," H3 " );
//VAL
 strcpy( AA[VAL].atom[10].atom_name, "HG11" );
  strcpy( AA[VAL].atom[11].atom_name, "HG12" );
  strcpy( AA[VAL].atom[12].atom_name, "HG13" );
  strcpy( AA[VAL].atom[13].atom_name, "HG21" );
  strcpy( AA[VAL].atom[14].atom_name, "HG22" );
  strcpy( AA[VAL].atom[15].atom_name, "HG23" );
   strcpy( AA[VAL].atom[17].atom_name," H2 " );
  strcpy( AA[VAL].atom[18].atom_name," H3 " );

//TRP
  strcpy( AA[TRP].atom[22].atom_name, " HB2" );
  strcpy( AA[TRP].atom[23].atom_name, " HB3" );
    strcpy( AA[TRP].atom[25].atom_name," H2 " );
  strcpy( AA[TRP].atom[26].atom_name," H3 " );

//TYR
 strcpy( AA[TYR].atom[19].atom_name, " HB2" );
 strcpy( AA[TYR].atom[20].atom_name, " HB3" );
 strcpy( AA[TYR].atom[22].atom_name," H2 " );
 strcpy( AA[TYR].atom[23].atom_name," H3 " );

//TYM
  strcpy( AA[TYM].atom[13].atom_name, " HD1" );
  strcpy( AA[TYM].atom[14].atom_name, " HE1" );
  strcpy( AA[TYM].atom[15].atom_name, " HE2" );
  strcpy( AA[TYM].atom[16].atom_name, " HD2" );
  strcpy( AA[TYM].atom[17].atom_name, " HA " );
  strcpy( AA[TYM].atom[18].atom_name, " HB2" );
  strcpy( AA[TYM].atom[19].atom_name, " HB2" );

}


void aa_iupac2pdb() {

  //fprintf( stderr, "ALA  \n" );

 //strcpy( AA[ALA].atom[5].atom_name,  " HN " );
 strcpy( AA[ALA].atom[11].atom_name,"1H  " );
 strcpy( AA[ALA].atom[12].atom_name,"2H  " );



//fprintf( stderr, "ARG  \n" );

 //strcpy( AA[ARG].atom[11].atom_name, " HN  " );
 strcpy( AA[ARG].atom[18].atom_name, "1HB " );
 strcpy( AA[ARG].atom[19].atom_name, "2HB " );
 strcpy( AA[ARG].atom[20].atom_name, "1HG " );
 strcpy( AA[ARG].atom[21].atom_name, "2HG " );
 strcpy( AA[ARG].atom[22].atom_name, "1HD " );
 strcpy( AA[ARG].atom[23].atom_name, "2HD " );
strcpy( AA[ARG].atom[25].atom_name,"1H  " );
  strcpy( AA[ARG].atom[26].atom_name,"2H  " );

 //fprintf( stderr, "ASN  \n" );
 //strcpy( AA[ASN].atom[8].atom_name,  " HN " );
 strcpy( AA[ASN].atom[12].atom_name, "1HB " );
 strcpy( AA[ASN].atom[13].atom_name, "2HB " );
 strcpy( AA[ASN].atom[15].atom_name,"1H  " );
  strcpy( AA[ASN].atom[16].atom_name,"2H  " );

 //fprintf( stderr, "ASP \n" );
 //strcpy( AA[ASP].atom[8].atom_name,  " HN " );
 strcpy( AA[ASP].atom[10].atom_name, "1HB " );
 strcpy( AA[ASP].atom[11].atom_name, "2HB " );
  strcpy( AA[ASP].atom[13].atom_name,"1H  " );
  strcpy( AA[ASP].atom[14].atom_name,"2H  " );

 //fprintf( stderr, "CYS  \n" );
 //strcpy( AA[CYS].atom[6].atom_name, " HN " );
 strcpy( AA[CYS].atom[8].atom_name, "1HB " );
 strcpy( AA[CYS].atom[9].atom_name, "2HB " );
  strcpy( AA[CYS].atom[12].atom_name,"1H  " );
  strcpy( AA[CYS].atom[13].atom_name,"2H  " );

 //fprintf( stderr, "GLN  \n" );
 //strcpy( AA[GLN].atom[9].atom_name,  " HN " );
 strcpy( AA[GLN].atom[13].atom_name, "1HB " );
 strcpy( AA[GLN].atom[14].atom_name, "2HB " );
 strcpy( AA[GLN].atom[15].atom_name, "1HG " );
 strcpy( AA[GLN].atom[16].atom_name, "2HG " );
 strcpy( AA[GLN].atom[18].atom_name,"1H  " );
  strcpy( AA[GLN].atom[19].atom_name,"2H  " );

 //fprintf( stderr, "GLUTI  \n" );
 //strcpy( AA[GLU].atom[9].atom_name,  " HN " );
 strcpy( AA[GLU].atom[11].atom_name, "1HB " );
 strcpy( AA[GLU].atom[12].atom_name, "2HB " );
 strcpy( AA[GLU].atom[13].atom_name, "1HG " );
 strcpy( AA[GLU].atom[14].atom_name, "2HG " );
 strcpy( AA[GLU].atom[16].atom_name,"1H  " );
  strcpy( AA[GLU].atom[17].atom_name,"2H  " );

 //fprintf( stderr, "GLY S \n" );
 //strcpy( AA[GLY].atom[4].atom_name, " HN " );
 strcpy( AA[GLY].atom[5].atom_name, "1HA " ); // ??
 strcpy( AA[GLY].atom[6].atom_name, "2HA " ); // ??
strcpy( AA[GLY].atom[8].atom_name,"1H  " );
  strcpy( AA[GLY].atom[9].atom_name,"2H  " );

 //fprintf( stderr, "HIS  \n" );
 //strcpy( AA[HIS].atom[10].atom_name, " HN " );
 strcpy( AA[HIS].atom[13].atom_name, "1HB " );
 strcpy( AA[HIS].atom[14].atom_name, "2HB " );
 strcpy( AA[HIS].atom[18].atom_name,"1H  " );
  strcpy( AA[HIS].atom[19].atom_name,"2H  " );

 //fprintf( stderr, "ILE  \n" );
 //strcpy( AA[ILE].atom[8].atom_name,  " HN " );
 strcpy( AA[ILE].atom[14].atom_name, "1HG1" );
 strcpy( AA[ILE].atom[15].atom_name, "2HG1" );
  strcpy( AA[ILE].atom[20].atom_name,"1H  " );
  strcpy( AA[ILE].atom[21].atom_name,"2H  " );

 //fprintf( stderr, "LEU  \n" );
 //strcpy( AA[LEU].atom[8].atom_name,  " HN " );
 strcpy( AA[LEU].atom[10].atom_name, "1HB " );
 strcpy( AA[LEU].atom[11].atom_name, "2HB " );
  strcpy( AA[LEU].atom[20].atom_name,"1H  " );
  strcpy( AA[LEU].atom[21].atom_name,"2H  " );

 //fprintf( stderr, "LYS  \n" );
 //strcpy( AA[LYS].atom[9].atom_name,  " HN " );
 strcpy( AA[LYS].atom[14].atom_name, "1HB " );
 strcpy( AA[LYS].atom[15].atom_name, "2HB " );
 strcpy( AA[LYS].atom[16].atom_name, "1HG " );
 strcpy( AA[LYS].atom[17].atom_name, "2HG " );
 strcpy( AA[LYS].atom[18].atom_name, "1HD " );
 strcpy( AA[LYS].atom[19].atom_name, "2HD " );
 strcpy( AA[LYS].atom[20].atom_name, "1HE " );
 strcpy( AA[LYS].atom[21].atom_name, "2HE " );
  strcpy( AA[LYS].atom[23].atom_name,"1H  " );
  strcpy( AA[LYS].atom[24].atom_name,"2H  " );

 //fprintf( stderr, "MET  \n" );
 //strcpy( AA[MET].atom[8].atom_name,  " HN " );
 strcpy( AA[MET].atom[10].atom_name, "1HB " );
 strcpy( AA[MET].atom[11].atom_name, "2HB " );
 strcpy( AA[MET].atom[12].atom_name, "1HG " );
 strcpy( AA[MET].atom[13].atom_name, "2HG " );
  strcpy( AA[MET].atom[18].atom_name,"1H  " );
  strcpy( AA[MET].atom[19].atom_name,"2H  " );

 //fprintf( stderr, "PHE  \n" );

 //strcpy( AA[PHE].atom[11].atom_name, " HN " );
 strcpy( AA[PHE].atom[18].atom_name, "1HB " );
 strcpy( AA[PHE].atom[19].atom_name, "2HB " );
  strcpy( AA[PHE].atom[21].atom_name,"1H  " );
  strcpy( AA[PHE].atom[22].atom_name,"2H  " );

 //fprintf( stderr, "PRO  \n" );
 strcpy( AA[PRO].atom[7].atom_name,  "1HD " );
 strcpy( AA[PRO].atom[8].atom_name,  "2HD " );
 strcpy( AA[PRO].atom[9].atom_name,  "1HG " );
 strcpy( AA[PRO].atom[10].atom_name, "2HG " );
 strcpy( AA[PRO].atom[11].atom_name, "1HB " );
 strcpy( AA[PRO].atom[12].atom_name, "2HB " );
  strcpy( AA[PRO].atom[15].atom_name,"1H  " );
  strcpy( AA[PRO].atom[16].atom_name,"2H  " );

 //fprintf( stderr, "SER \n" );
 //strcpy( AA[SER].atom[6].atom_name,  " HN " );
 strcpy( AA[SER].atom[9].atom_name,  "1HB " );
 strcpy( AA[SER].atom[10].atom_name, "2HB " );
strcpy( AA[SER].atom[12].atom_name,"1H  " );
  strcpy( AA[SER].atom[13].atom_name,"2H  " );

 //fprintf( stderr, "THR  \n" );
 //strcpy( AA[THR].atom[7].atom_name,  " HN " );
 strcpy( AA[THR].atom[15].atom_name,"1H  " );
  strcpy( AA[THR].atom[16].atom_name,"2H  " );

 //strcpy( AA[TRP].atom[14].atom_name, " HN " );
 strcpy( AA[TRP].atom[22].atom_name, "1HB " );
 strcpy( AA[TRP].atom[23].atom_name, "2HB " );
  strcpy( AA[TRP].atom[25].atom_name,"1H  " );
  strcpy( AA[TRP].atom[26].atom_name,"2H  " );

 //strcpy( AA[TYR].atom[12].atom_name, " HN " );
 strcpy( AA[TYR].atom[19].atom_name, "1HB " );
 strcpy( AA[TYR].atom[20].atom_name, "2HB " );
strcpy( AA[TYR].atom[22].atom_name,"1H  " );
  strcpy( AA[TYR].atom[23].atom_name,"2H  " );

 //strcpy( AA[VAL].atom[7].atom_name,  " HN " );
strcpy( AA[VAL].atom[17].atom_name,"1H  " );
  strcpy( AA[VAL].atom[18].atom_name,"2H  " );

	//ASH
	 strcpy( AA[ASH].atom[10].atom_name, "1HB " );
  strcpy( AA[ASH].atom[11].atom_name, "2HB " );
  strcpy( AA[ASH].atom[12].atom_name, "1HD " ); //NEW ATOM
  strcpy( AA[ASH].atom[14].atom_name,"1H  " );
  strcpy( AA[ASH].atom[15].atom_name,"2H  " );
  //CYX
  strcpy( AA[CYX].atom[8].atom_name, "1HB " );
  strcpy( AA[CYX].atom[9].atom_name, "2HB " );
  strcpy( AA[CYX].atom[11].atom_name,"1H  " );
  strcpy( AA[CYX].atom[12].atom_name,"2H  " );
  //CYM
  strcpy( AA[CYM].atom[8].atom_name, "1HB " );
  strcpy( AA[CYM].atom[9].atom_name, "2HB " );
  strcpy( AA[CYM].atom[11].atom_name,"1H  " );
  strcpy( AA[CYM].atom[12].atom_name,"2H  " );
  //GLH
   strcpy( AA[GLH].atom[11].atom_name, "1HB " );
  strcpy( AA[GLH].atom[12].atom_name, "2HB " );
  strcpy( AA[GLH].atom[13].atom_name, "1HG " );
  strcpy( AA[GLH].atom[14].atom_name, "2HG " );
  strcpy( AA[GLH].atom[15].atom_name, "1HE " );//NEW ADDED
  strcpy( AA[GLH].atom[17].atom_name,"1H  " );
  strcpy( AA[GLH].atom[18].atom_name,"2H  " );

  //HIP
  strcpy( AA[HIP].atom[11].atom_name," HE2" );// <--- The proton.????????
  strcpy( AA[HIP].atom[13].atom_name,"1HB " );
  strcpy( AA[HIP].atom[14].atom_name,"2HB " );
  strcpy( AA[HIP].atom[15].atom_name," HE1" ); // <--- The proton.????????
  strcpy( AA[HIP].atom[16].atom_name," HD2" );
  strcpy( AA[HIP].atom[17].atom_name," HD1" ); // NEW ADDED
	strcpy( AA[HIP].atom[19].atom_name,"1H  " );
  strcpy( AA[HIP].atom[20].atom_name,"2H  " );
  //HID
  strcpy( AA[HID].atom[11].atom_name," HD1" );//NEW ADDED
  strcpy( AA[HID].atom[12].atom_name," HA " );
  strcpy( AA[HID].atom[13].atom_name,"1HB " );
  strcpy( AA[HID].atom[14].atom_name,"2HB " );
  strcpy( AA[HID].atom[15].atom_name," HE1" ); // <--- The proton.????????
  strcpy( AA[HID].atom[16].atom_name," HD2" );
  strcpy( AA[HID].atom[18].atom_name,"1H  " );
  strcpy( AA[HID].atom[19].atom_name,"2H  " );
  //HIE
  strcpy( AA[HIE].atom[11].atom_name," HE2" );// <--- The proton.????????
  strcpy( AA[HIE].atom[13].atom_name,"1HB " );
  strcpy( AA[HIE].atom[14].atom_name,"2HB " );
  strcpy( AA[HIE].atom[15].atom_name," HE1" ); // <--- The proton.????????
  strcpy( AA[HIE].atom[16].atom_name," HD2" );
	strcpy( AA[HIE].atom[18].atom_name,"1H  " );
  strcpy( AA[HIE].atom[19].atom_name,"2H  " );
  //LYN
  strcpy( AA[LYN].atom[10].atom_name, "1HZ " );
  strcpy( AA[LYN].atom[11].atom_name, "2HZ " );
  strcpy( AA[LYN].atom[12].atom_name, " HA " );
  strcpy( AA[LYN].atom[13].atom_name, "1HB " );
  strcpy( AA[LYN].atom[14].atom_name, "2HB " );
  strcpy( AA[LYN].atom[15].atom_name, "1HG " );
  strcpy( AA[LYN].atom[16].atom_name, "2HG " );
  strcpy( AA[LYN].atom[17].atom_name, "1HD " );
  strcpy( AA[LYN].atom[18].atom_name, "2HD " );
  strcpy( AA[LYN].atom[19].atom_name, "1HE " );
  strcpy( AA[LYN].atom[20].atom_name, "2HE " );
  strcpy( AA[LYN].atom[22].atom_name,"1H  " );
  strcpy( AA[LYN].atom[23].atom_name,"2H  " );
  //TYM
  strcpy( AA[TYM].atom[13].atom_name, " HD1" );
  strcpy( AA[TYM].atom[14].atom_name, " HE1" );
  strcpy( AA[TYM].atom[15].atom_name, " HE2" );
  strcpy( AA[TYM].atom[16].atom_name, " HD2" );
  strcpy( AA[TYM].atom[17].atom_name, " HA " );
  strcpy( AA[TYM].atom[18].atom_name, "1HB " );
  strcpy( AA[TYM].atom[19].atom_name, "2HB " );
  strcpy( AA[TYM].atom[21].atom_name,"1H  " );
  strcpy( AA[TYM].atom[22].atom_name,"2H  " );

}


void aa_pdb2iupac() {

 //strcpy( AA[ALA].atom[5].atom_name,  " H  " );
 strcpy( AA[ALA].atom[11].atom_name,"2H  " );
 strcpy( AA[ALA].atom[12].atom_name,"3H  " );

 //CYS
 //strcpy( AA[CYS].atom[6].atom_name, " H  " );
 strcpy( AA[CYS].atom[8].atom_name, "2HB " );
 strcpy( AA[CYS].atom[9].atom_name, "3HB " );
  strcpy( AA[CYS].atom[12].atom_name,"2H  " );
  strcpy( AA[CYS].atom[13].atom_name,"3H  " );

 //ASP
  //strcpy( AA[ASP].atom[8].atom_name,  " H  " );
 	strcpy( AA[ASP].atom[10].atom_name, "2HB " );
 	strcpy( AA[ASP].atom[11].atom_name, "3HB " );
  strcpy( AA[ASP].atom[13].atom_name,"2H  " );
  strcpy( AA[ASP].atom[14].atom_name,"3H  " );
 //GLU
 //strcpy( AA[GLU].atom[9].atom_name, " H  " );
 strcpy( AA[GLU].atom[11].atom_name, "2HB " );
 strcpy( AA[GLU].atom[12].atom_name, "3HB " );
 strcpy( AA[GLU].atom[13].atom_name, "2HG " );
 strcpy( AA[GLU].atom[14].atom_name, "3HG " );
  strcpy( AA[GLU].atom[16].atom_name,"2H  " );
  strcpy( AA[GLU].atom[17].atom_name,"3H  " );
 //PHE
 //strcpy( AA[PHE].atom[11].atom_name, " H  " );
 strcpy( AA[PHE].atom[17].atom_name, " HA " );
 strcpy( AA[PHE].atom[18].atom_name, "2HB " );
 strcpy( AA[PHE].atom[19].atom_name, "3HB " );
  strcpy( AA[PHE].atom[21].atom_name,"2H  " );
  strcpy( AA[PHE].atom[22].atom_name,"3H  " );
 //GLY
 //strcpy( AA[GLY].atom[4].atom_name, " H  " );
 strcpy( AA[GLY].atom[5].atom_name, "2HA " );
 strcpy( AA[GLY].atom[6].atom_name, "3HA " );
  strcpy( AA[GLY].atom[8].atom_name,"2H  " );
  strcpy( AA[GLY].atom[9].atom_name,"3H  " );
 //HIS
 //strcpy( AA[HIS].atom[10].atom_name, " H  " );
 strcpy( AA[HIS].atom[13].atom_name, "2HB " );
 strcpy( AA[HIS].atom[14].atom_name, "3HB " );
  strcpy( AA[HIS].atom[18].atom_name,"2H  " );
  strcpy( AA[HIS].atom[19].atom_name,"3H  " );
 //ILE
 //strcpy( AA[ILE].atom[8].atom_name,  " H  " );
 strcpy( AA[ILE].atom[14].atom_name, "2HG1" );
 strcpy( AA[ILE].atom[15].atom_name, "3HG1" );
  strcpy( AA[ILE].atom[20].atom_name,"2H  " );
  strcpy( AA[ILE].atom[21].atom_name,"3H  " );
 //LYS
 //strcpy( AA[LYS].atom[9].atom_name, " H  " );
 strcpy( AA[LYS].atom[14].atom_name, "2HB " );
 strcpy( AA[LYS].atom[15].atom_name, "3HB " );
 strcpy( AA[LYS].atom[16].atom_name, "2HG " );
 strcpy( AA[LYS].atom[17].atom_name, "3HG " );
 strcpy( AA[LYS].atom[18].atom_name, "2HD " );
 strcpy( AA[LYS].atom[19].atom_name, "3HD " );
 strcpy( AA[LYS].atom[20].atom_name, "2HE " );
 strcpy( AA[LYS].atom[21].atom_name, "3HE " );
  strcpy( AA[LYS].atom[23].atom_name,"2H  " );
  strcpy( AA[LYS].atom[24].atom_name,"3H  " );
 //LEU
 //strcpy( AA[LEU].atom[8].atom_name, " H  " );
 strcpy( AA[LEU].atom[10].atom_name, "2HB " );
 strcpy( AA[LEU].atom[11].atom_name, "3HB " );
  strcpy( AA[LEU].atom[20].atom_name,"2H  " );
  strcpy( AA[LEU].atom[21].atom_name,"3H  " );
  //MET
   //strcpy( AA[MET].atom[8].atom_name, " H  " );
 strcpy( AA[MET].atom[10].atom_name, "2HB " );
 strcpy( AA[MET].atom[11].atom_name, "3HB " );
 strcpy( AA[MET].atom[12].atom_name, "2HG " );
 strcpy( AA[MET].atom[13].atom_name, "3HG " );
  strcpy( AA[MET].atom[18].atom_name,"2H  " );
  strcpy( AA[MET].atom[19].atom_name,"3H  " );
  //ASN
 //fprintf( stderr, "ASN  \n" );
  //strcpy( AA[ASN].atom[8].atom_name,  " H  " );
  strcpy( AA[ASN].atom[12].atom_name, "2HB " );
  strcpy( AA[ASN].atom[13].atom_name, "3HB " );
  strcpy( AA[ASN].atom[15].atom_name,"2H  " );
  strcpy( AA[ASN].atom[16].atom_name,"3H  " );
  //PRO
 //fprintf( stderr, "PRO  \n" );
 //strcpy( AA[PRO].atom[7].atom_name, "2HD " );
 strcpy( AA[PRO].atom[8].atom_name, "3HD " );
 strcpy( AA[PRO].atom[9].atom_name, "2HG " );
 strcpy( AA[PRO].atom[10].atom_name, "3HG " );
 strcpy( AA[PRO].atom[11].atom_name, "2HB " );
 strcpy( AA[PRO].atom[12].atom_name, "3HB " );
  strcpy( AA[PRO].atom[15].atom_name,"2H  " );
  strcpy( AA[PRO].atom[16].atom_name,"3H  " );
  //GLN
// fprintf( stderr, "GLN  \n" );
 //strcpy( AA[GLN].atom[9].atom_name,  " H  " );
  strcpy( AA[GLN].atom[13].atom_name, "2HB " );
 strcpy( AA[GLN].atom[14].atom_name, "3HB " );
 strcpy( AA[GLN].atom[15].atom_name, "2HG " );
  strcpy( AA[GLN].atom[16].atom_name, "3HG " );
  strcpy( AA[GLN].atom[18].atom_name,"2H  " );
  strcpy( AA[GLN].atom[19].atom_name,"3H  " );
  //ARG
 //fprintf( stderr, "ARG  \n" );
  //strcpy( AA[ARG].atom[11].atom_name, " H  " );
  strcpy( AA[ARG].atom[18].atom_name, "2HB " );
  strcpy( AA[ARG].atom[19].atom_name, "3HB " );
  strcpy( AA[ARG].atom[20].atom_name, "2HG " );
  strcpy( AA[ARG].atom[21].atom_name, "3HG " );
  strcpy( AA[ARG].atom[22].atom_name, "2HD " );
  strcpy( AA[ARG].atom[23].atom_name, "3HD " );
  strcpy( AA[ARG].atom[25].atom_name,"2H  " );
  strcpy( AA[ARG].atom[26].atom_name,"3H  " );
  //SER
 //fprintf( stderr, "SER  \n" );
 //strcpy( AA[SER].atom[6].atom_name, " H  " );
 strcpy( AA[SER].atom[9].atom_name, "2HB " );
 strcpy( AA[SER].atom[10].atom_name, "3HB " );
  strcpy( AA[SER].atom[12].atom_name,"2H  " );
  strcpy( AA[SER].atom[13].atom_name,"3H  " );
  //THR
 //fprintf( stderr, "THR  \n" );
  strcpy( AA[THR].atom[8].atom_name, "1HG " );
  strcpy( AA[THR].atom[9].atom_name, " HA " );
  strcpy( AA[THR].atom[10].atom_name," HB " );
  strcpy( AA[THR].atom[11].atom_name,"1HG2" );
  strcpy( AA[THR].atom[12].atom_name,"2HG2" );
  strcpy( AA[THR].atom[13].atom_name,"3HG2" );
  strcpy( AA[THR].atom[15].atom_name,"2H  " );
  strcpy( AA[THR].atom[16].atom_name,"3H  " );

  //VAL
 //fprintf( stderr, "VAL  \n" );
 //strcpy( AA[VAL].atom[7].atom_name, " H  " );
  strcpy( AA[VAL].atom[17].atom_name,"2H  " );
  strcpy( AA[VAL].atom[18].atom_name,"3H  " );
  //TRP
 //strcpy( AA[TRP].atom[14].atom_name, " H  " );
 //fprintf( stderr, "TRP  \n" );
 strcpy( AA[TRP].atom[22].atom_name, "2HB " );
 strcpy( AA[TRP].atom[23].atom_name, "3HB " );
  strcpy( AA[TRP].atom[25].atom_name,"2H  " );
  strcpy( AA[TRP].atom[26].atom_name,"3H  " );
  //TYR
 //strcpy( AA[TYR].atom[12].atom_name, " H  " );
 //fprintf( stderr, "TYR  \n" );
 strcpy( AA[TYR].atom[19].atom_name, "2HB " );
 strcpy( AA[TYR].atom[20].atom_name, "3HB " );
  strcpy( AA[TYR].atom[22].atom_name,"2H  " );
  strcpy( AA[TYR].atom[23].atom_name,"3H  " );
  //ASH
 //fprintf( stderr, "ASH  \n" );
	 strcpy( AA[ASH].atom[10].atom_name, "2HB " );
  strcpy( AA[ASH].atom[11].atom_name, "3HB " );
  strcpy( AA[ASH].atom[12].atom_name, "2HD " ); //NEW ATOM
  strcpy( AA[ASH].atom[14].atom_name,"2H  " );
  strcpy( AA[ASH].atom[15].atom_name,"3H  " );
  //CYX
 //fprintf( stderr, "CYX  \n" );
  strcpy( AA[CYX].atom[8].atom_name, "2HB " );
  strcpy( AA[CYX].atom[9].atom_name, "3HB " );
  strcpy( AA[CYX].atom[11].atom_name,"2H  " );
  strcpy( AA[CYX].atom[12].atom_name,"3H  " );
  //CYM
 //fprintf( stderr, "CYM  \n" );
  strcpy( AA[CYM].atom[8].atom_name, "2HB " );
  strcpy( AA[CYM].atom[9].atom_name, "3HB " );
  strcpy( AA[CYM].atom[11].atom_name,"2H  " );
  strcpy( AA[CYM].atom[12].atom_name,"3H  " );
  //GLH
 //fprintf( stderr, "GLH  \n" );
   strcpy( AA[GLH].atom[11].atom_name, "2HB " );
  strcpy( AA[GLH].atom[12].atom_name, "3HB " );
  strcpy( AA[GLH].atom[13].atom_name, "2HG " );
  strcpy( AA[GLH].atom[14].atom_name, "3HG " );
  strcpy( AA[GLH].atom[15].atom_name, "2HE " );//NEW ADDED
  strcpy( AA[GLH].atom[17].atom_name,"2H  " );
  strcpy( AA[GLH].atom[18].atom_name,"3H  " );
  //HIP
 //fprintf( stderr, "HIP  \n" );
  strcpy( AA[HIP].atom[10].atom_name," H  " );
  strcpy( AA[HIP].atom[11].atom_name," HE2" ); // <--- The proton.
  strcpy( AA[HIP].atom[12].atom_name," HD1" ); // <--- The proton.
  strcpy( AA[HIP].atom[13].atom_name," HA " );
  strcpy( AA[HIP].atom[14].atom_name,"2HB " );
  strcpy( AA[HIP].atom[15].atom_name,"3HB " );
  strcpy( AA[HIP].atom[16].atom_name," HE1" );
  strcpy( AA[HIP].atom[17].atom_name," HD2" );
	strcpy( AA[HIP].atom[19].atom_name,"2H  " );
  strcpy( AA[HIP].atom[20].atom_name,"3H  " );
  //HID
 //fprintf( stderr, "HID  \n" );
  strcpy( AA[HID].atom[18].atom_name,"2H  " );
  strcpy( AA[HID].atom[19].atom_name,"3H  " );
  strcpy( AA[HID].atom[11].atom_name," HD1" );//NEW ADDED
  strcpy( AA[HID].atom[12].atom_name," HA " );
  strcpy( AA[HID].atom[13].atom_name,"2HB " );
  strcpy( AA[HID].atom[14].atom_name,"3HB " );
  strcpy( AA[HID].atom[15].atom_name," HE1" ); // <--- The proton.????????
  strcpy( AA[HID].atom[16].atom_name," HD2" );


  //HIE
 //fprintf( stderr, "HIE  \n" );
  strcpy( AA[HIE].atom[11].atom_name," HE2" );// <--- The proton.????????
  strcpy( AA[HIE].atom[13].atom_name,"2HB " );
  strcpy( AA[HIE].atom[14].atom_name,"3HB " );
  strcpy( AA[HIE].atom[15].atom_name," HE1" ); // <--- The proton.????????
  strcpy( AA[HIE].atom[16].atom_name," HD2" );
	strcpy( AA[HIE].atom[18].atom_name,"2H  " );
  strcpy( AA[HIE].atom[19].atom_name,"3H  " );
  //LYN
  //fprintf( stderr, "LYN  \n" );
  strcpy( AA[LYN].atom[10].atom_name, "1HZ " );
  strcpy( AA[LYN].atom[11].atom_name, "2HZ " );
  strcpy( AA[LYN].atom[12].atom_name, " HA " );
  strcpy( AA[LYN].atom[13].atom_name, "2HB " );
  strcpy( AA[LYN].atom[14].atom_name, "3HB " );
  strcpy( AA[LYN].atom[15].atom_name, "2HG " );
  strcpy( AA[LYN].atom[16].atom_name, "3HG " );
  strcpy( AA[LYN].atom[17].atom_name, "2HD " );
  strcpy( AA[LYN].atom[18].atom_name, "3HD " );
  strcpy( AA[LYN].atom[19].atom_name, "2HE " );
  strcpy( AA[LYN].atom[20].atom_name, "3HE " );
  strcpy( AA[LYN].atom[22].atom_name,"2H  " );
  strcpy( AA[LYN].atom[23].atom_name,"3H  " );
  //TYM
  //fprintf( stderr, "TYM  \n" );

  strcpy( AA[TYM].atom[13].atom_name, " HD1" );
  strcpy( AA[TYM].atom[14].atom_name, " HE1" );
  strcpy( AA[TYM].atom[15].atom_name, " HE2" );
  strcpy( AA[TYM].atom[16].atom_name, " HD2" );
  strcpy( AA[TYM].atom[17].atom_name, " HA " );
  strcpy( AA[TYM].atom[18].atom_name, "2HB " );
  strcpy( AA[TYM].atom[19].atom_name, "3HB " );
  strcpy( AA[TYM].atom[21].atom_name,"2H  " );
  strcpy( AA[TYM].atom[22].atom_name,"3H  " );

}




void nucleotids(Convention opt)
{

  //GUA
//	  int k, p = GUA;
  int k, p = DGUA;

//  strcpy( AA[p].aa_name3, "  G" );
  strcpy( AA[p].aa_name3, " DG" );
//  AA[p].aa_name1 = 'X';
  AA[p].aa_name1 = 'G';
  AA[p].mass = 0;

  AA[p].is_protein = false;
  AA[p].is_RNA = false;
  AA[p].is_DNA = true;

  AA[p].aa_is_polar = true;
  AA[p].aa_is_nonpolar = false;
  AA[p].aa_is_aromatic = true;
  AA[p].aa_is_charged = true;

  AA[p].natoms = 33;
  AA[p].nheavyatoms = 22;
  AA[p].nchi = 0;

  AA[p].chi_atoms = NULL;
  AA[p].atom = ( t_aa_atom_p * ) malloc( AA[p].natoms * sizeof( t_aa_atom_p ) );


  strcpy( AA[p].atom[0].atom_name, " P  " );
//  strcpy( AA[p].atom[1].atom_name, " O1P" );
//  strcpy( AA[p].atom[2].atom_name, " O2P" );

  // New PDB format for phosphate oxygens
  strcpy( AA[p].atom[1].atom_name, " OP1" );
  strcpy( AA[p].atom[2].atom_name, " OP2" );


//	strcpy( AA[p].atom[3].atom_name, " O5*" );
//	strcpy( AA[p].atom[4].atom_name, " C5*" );
//	strcpy( AA[p].atom[5].atom_name, " C4*" );
//	strcpy( AA[p].atom[6].atom_name, " O4*" );
//	strcpy( AA[p].atom[7].atom_name, " C3*" );
//	strcpy( AA[p].atom[8].atom_name, " O3*" );
//	strcpy( AA[p].atom[9].atom_name, " C2*" );
//	strcpy( AA[p].atom[10].atom_name, " C1*" );

  // New PDB format for sugar...
  strcpy( AA[p].atom[3].atom_name, " O5'" );
  strcpy( AA[p].atom[4].atom_name, " C5'" );
  strcpy( AA[p].atom[5].atom_name, " C4'" );
  strcpy( AA[p].atom[6].atom_name, " O4'" );
  strcpy( AA[p].atom[7].atom_name, " C3'" );
  strcpy( AA[p].atom[8].atom_name, " O3'" );
  strcpy( AA[p].atom[9].atom_name, " C2'" );
  strcpy( AA[p].atom[10].atom_name, " C1'" );

  strcpy( AA[p].atom[11].atom_name, " N1 " );
  strcpy( AA[p].atom[12].atom_name, " C2 " );
  strcpy( AA[p].atom[13].atom_name, " N2 " );
  strcpy( AA[p].atom[14].atom_name, " N3 " );
  strcpy( AA[p].atom[15].atom_name, " C4 " );
  strcpy( AA[p].atom[16].atom_name, " C5 " );
  strcpy( AA[p].atom[17].atom_name, " C6 " );
  strcpy( AA[p].atom[18].atom_name, " O6 " );
  strcpy( AA[p].atom[19].atom_name, " N7 " );
  strcpy( AA[p].atom[20].atom_name, " C8 " );
  strcpy( AA[p].atom[21].atom_name, " N9 " );


//  strcpy( AA[p].atom[22].atom_name, "1H5*" );
//  strcpy( AA[p].atom[23].atom_name, "2H5*" );
//  strcpy( AA[p].atom[24].atom_name, " H4*" );
//  strcpy( AA[p].atom[25].atom_name, " H3*" );
//  strcpy( AA[p].atom[26].atom_name, "1H2*" );
//  strcpy( AA[p].atom[27].atom_name, "2H2*" );
//  strcpy( AA[p].atom[28].atom_name, " H1*" );

  // New PDB format for sugar...
  strcpy( AA[p].atom[22].atom_name, "1H5'" );
  strcpy( AA[p].atom[23].atom_name, "2H5'" );
  strcpy( AA[p].atom[24].atom_name, " H4'" );
  strcpy( AA[p].atom[25].atom_name, " H3'" );
  strcpy( AA[p].atom[26].atom_name, "1H2'" );
  strcpy( AA[p].atom[27].atom_name, "2H2'" );
  strcpy( AA[p].atom[28].atom_name, " H1'" );

  strcpy( AA[p].atom[29].atom_name, " H1 " );
  strcpy( AA[p].atom[30].atom_name, "1H2 " );
  strcpy( AA[p].atom[31].atom_name, "2H2 " );
  strcpy( AA[p].atom[32].atom_name, " H8 " );


  AA[p].atom[0].charge = 1.500; // P
  AA[p].atom[1].charge = -.780; // O1P
  AA[p].atom[2].charge = -.780; // O2P
  AA[p].atom[3].charge = -.570; // O5*
  AA[p].atom[4].charge = -.080; // C5*
  AA[p].atom[5].charge =.160; // C4*
  AA[p].atom[6].charge = -.500; // O4*
  AA[p].atom[7].charge =.010; // C3*
  AA[p].atom[8].charge = -.570; // O3*
  AA[p].atom[9].charge = -.180; // C2* !!!
  AA[p].atom[10].charge =.160; // C1*
  AA[p].atom[11].charge = -.340; // N1
  AA[p].atom[12].charge =.750; // C2
  AA[p].atom[13].charge = -.680; // N2
  AA[p].atom[14].charge = -.740; // N3
  AA[p].atom[15].charge =.260; // C4
  AA[p].atom[16].charge =.000; // C5
  AA[p].atom[17].charge =.540; // C6
  AA[p].atom[18].charge = -.510; // O6
  AA[p].atom[19].charge = -.600; // N7
  AA[p].atom[20].charge =.250; // C8
  AA[p].atom[21].charge = -.020; // N9
  AA[p].atom[22].charge =.090; // 1H5*
  AA[p].atom[23].charge =.090; // 2H5*
  AA[p].atom[24].charge =.090; // H4*
  AA[p].atom[25].charge =.090; // H3*
  AA[p].atom[26].charge =.090; // 1H2*
  AA[p].atom[27].charge =.090; // 2H2*
  AA[p].atom[28].charge =.090; // H1*
  AA[p].atom[29].charge =.260; // H1
  AA[p].atom[30].charge =.320; // 1H2
  AA[p].atom[31].charge =.350; // 2H2
  AA[p].atom[32].charge =.160; // H8

  switch (opt)
  {
  //NEW ADDED. MODIFIED BY GLUP
  case Sybil:
	  AA[p].atom[0].fullatom_type = 12; // Phos   P
	  AA[p].atom[1].fullatom_type = 10; // OOC    O1P
	  AA[p].atom[2].fullatom_type = 10; // OOC    O2P
	  AA[p].atom[3].fullatom_type = 9; // ONH2   O5*
	  AA[p].atom[4].fullatom_type = 2; // CH2    C5*
	  AA[p].atom[5].fullatom_type = 2; // CH1    C4*
	  AA[p].atom[6].fullatom_type = 9; // OH     O4*
	  AA[p].atom[7].fullatom_type = 2; // CH1    C3*
	  AA[p].atom[8].fullatom_type = 9; // ONH2   O3*
	  AA[p].atom[9].fullatom_type = 2; // CH2    C2*
	  AA[p].atom[10].fullatom_type = 2; // CH1    C1*
	  AA[p].atom[11].fullatom_type = 5; // Ntrp   N1
	  AA[p].atom[12].fullatom_type = 3; // aroC   C2
	  AA[p].atom[13].fullatom_type = 5; // NH2O   N2
	  AA[p].atom[14].fullatom_type = 5; // Nhis   N3
	  AA[p].atom[15].fullatom_type = 3; // aroC   C4
	  AA[p].atom[16].fullatom_type = 3; // aroC   C5
	  AA[p].atom[17].fullatom_type = 1; // CObb   C6
	  AA[p].atom[18].fullatom_type = 8; // OCbb   O6
	  AA[p].atom[19].fullatom_type = 5; // Nhis   N7
	  AA[p].atom[20].fullatom_type = 3; // aroC   C8
	  AA[p].atom[21].fullatom_type = 5; // NHis   N9
	  AA[p].atom[22].fullatom_type = 14; // Hapo  1H5*
	  AA[p].atom[23].fullatom_type = 14; // Hapo  2H5*
	  AA[p].atom[24].fullatom_type = 14; // Hapo   H4*
	  AA[p].atom[25].fullatom_type = 14; // Hapo   H3*
	  AA[p].atom[26].fullatom_type = 14; // Hapo  1H2*
	  AA[p].atom[27].fullatom_type = 14; // Hapo  2H2*
	  AA[p].atom[28].fullatom_type = 14; // Hapo   H1*
	  AA[p].atom[29].fullatom_type = 13; // Hpol   H1
	  AA[p].atom[30].fullatom_type = 13; // Hpol  1H2
	  AA[p].atom[31].fullatom_type = 13; // Hpol  2H2
	  AA[p].atom[32].fullatom_type = 15; // Haro   H8
	  break;

  default:
	  AA[p].atom[0].fullatom_type = 21; // Phos   P
	  AA[p].atom[1].fullatom_type = 15; // OOC    O1P
	  AA[p].atom[2].fullatom_type = 15; // OOC    O2P
	  AA[p].atom[3].fullatom_type = 14; // ONH2   O5*
	  AA[p].atom[4].fullatom_type = 4; // CH2    C5*
	  AA[p].atom[5].fullatom_type = 3; // CH1    C4*
	  AA[p].atom[6].fullatom_type = 13; // OH     O4*
	  AA[p].atom[7].fullatom_type = 3; // CH1    C3*
	  AA[p].atom[8].fullatom_type = 14; // ONH2   O3*
	  AA[p].atom[9].fullatom_type = 4; // CH2    C2*
	  AA[p].atom[10].fullatom_type = 3; // CH1    C1*
	  AA[p].atom[11].fullatom_type = 7; // Ntrp   N1
	  AA[p].atom[12].fullatom_type = 6; // aroC   C2
	  AA[p].atom[13].fullatom_type = 9; // NH2O   N2
	  AA[p].atom[14].fullatom_type = 8; // Nhis   N3
	  AA[p].atom[15].fullatom_type = 6; // aroC   C4
	  AA[p].atom[16].fullatom_type = 6; // aroC   C5
	  AA[p].atom[17].fullatom_type = 19; // CObb   C6
	  AA[p].atom[18].fullatom_type = 20; // OCbb   O6
	  AA[p].atom[19].fullatom_type = 8; // Nhis   N7
	  AA[p].atom[20].fullatom_type = 6; // aroC   C8
	  AA[p].atom[21].fullatom_type = 8; // NHis   N9
	  AA[p].atom[22].fullatom_type = 23; // Hapo  1H5*
	  AA[p].atom[23].fullatom_type = 23; // Hapo  2H5*
	  AA[p].atom[24].fullatom_type = 23; // Hapo   H4*
	  AA[p].atom[25].fullatom_type = 23; // Hapo   H3*
	  AA[p].atom[26].fullatom_type = 23; // Hapo  1H2*
	  AA[p].atom[27].fullatom_type = 23; // Hapo  2H2*
	  AA[p].atom[28].fullatom_type = 23; // Hapo   H1*
	  AA[p].atom[29].fullatom_type = 22; // Hpol   H1
	  AA[p].atom[30].fullatom_type = 22; // Hpol  1H2
	  AA[p].atom[31].fullatom_type = 22; // Hpol  2H2
	  AA[p].atom[32].fullatom_type = 24; // Haro   H8
	  break;
  }

  //jjh intra residue bonding

  AA[p].atom[0].nbonded_neighbors = 3; // P
  AA[p].atom[0].bonded_neighbor[0] = 1; // P--O1P
  AA[p].atom[0].bonded_neighbor[1] = 2; // P--O2P
  AA[p].atom[0].bonded_neighbor[2] = 3; // P--O5*
  AA[p].atom[1].nbonded_neighbors = 1; // O1P
  AA[p].atom[1].bonded_neighbor[0] = 0; // O1P--P
  AA[p].atom[2].nbonded_neighbors = 1; // O2P
  AA[p].atom[2].bonded_neighbor[0] = 0; // O2P--P
  AA[p].atom[3].nbonded_neighbors = 2; // O5*
  AA[p].atom[3].bonded_neighbor[0] = 0; // O5*--P
  AA[p].atom[3].bonded_neighbor[1] = 4; // O5*--C5*
  AA[p].atom[4].nbonded_neighbors = 4; // C5*
  AA[p].atom[4].bonded_neighbor[0] = 3; // C5*--O5*
  AA[p].atom[4].bonded_neighbor[1] = 5; // C5*--C4*
  AA[p].atom[4].bonded_neighbor[2] = 22; // C5*--1H5*
  AA[p].atom[4].bonded_neighbor[3] = 23; // C5*--2H5*
  AA[p].atom[5].nbonded_neighbors = 4; // C4*
  AA[p].atom[5].bonded_neighbor[0] = 4; // C4*--C5*
  AA[p].atom[5].bonded_neighbor[1] = 6; // C4*--O4*
  AA[p].atom[5].bonded_neighbor[2] = 7; // C4*--C3*
  AA[p].atom[5].bonded_neighbor[3] = 24; // C4*--H4*
  AA[p].atom[6].nbonded_neighbors = 2; // O4*
  AA[p].atom[6].bonded_neighbor[0] = 5; // O4*--C4*
  AA[p].atom[6].bonded_neighbor[1] = 10; // O4*--C1*
  AA[p].atom[7].nbonded_neighbors = 4; // C3*
  AA[p].atom[7].bonded_neighbor[0] = 5; // C3*--C4*
  AA[p].atom[7].bonded_neighbor[1] = 8; // C3*--O3*
  AA[p].atom[7].bonded_neighbor[2] = 9; // C3*--C2*
  AA[p].atom[7].bonded_neighbor[3] = 25; // C3*--H3*
  AA[p].atom[8].nbonded_neighbors = 1; // O3*
  AA[p].atom[8].bonded_neighbor[0] = 7; // O3*--C3*
  AA[p].atom[9].nbonded_neighbors = 4; // C2*
  AA[p].atom[9].bonded_neighbor[0] = 7; // C2*--C3*
  AA[p].atom[9].bonded_neighbor[1] = 10; // C2*--C1*
  AA[p].atom[9].bonded_neighbor[2] = 26; // C2*--1H2*
  AA[p].atom[9].bonded_neighbor[3] = 27; // C2*--2H2*
  AA[p].atom[10].nbonded_neighbors = 4; // C1*
  AA[p].atom[10].bonded_neighbor[0] = 6; // C1*--O4*
  AA[p].atom[10].bonded_neighbor[1] = 9; // C1*--C2*
  AA[p].atom[10].bonded_neighbor[2] = 21; // C1*--N9
  AA[p].atom[10].bonded_neighbor[3] = 28; // C1*--H1*
  AA[p].atom[11].nbonded_neighbors = 3; // N1
  AA[p].atom[11].bonded_neighbor[0] = 12; // N1--C2
  AA[p].atom[11].bonded_neighbor[1] = 17; // N1--C6
  AA[p].atom[11].bonded_neighbor[2] = 29; // N1--H1
  AA[p].atom[12].nbonded_neighbors = 3; // C2
  AA[p].atom[12].bonded_neighbor[0] = 11; // C2--N1
  AA[p].atom[12].bonded_neighbor[1] = 13; // C2--N2
  AA[p].atom[12].bonded_neighbor[2] = 14; // C2--N3
  AA[p].atom[13].nbonded_neighbors = 3; // N2
  AA[p].atom[13].bonded_neighbor[0] = 12; // N2--C2
  AA[p].atom[13].bonded_neighbor[1] = 30; // N2--1H2
  AA[p].atom[13].bonded_neighbor[2] = 31; // N2--2H2
  AA[p].atom[14].nbonded_neighbors = 2; // N3
  AA[p].atom[14].bonded_neighbor[0] = 12; // N3--C2
  AA[p].atom[14].bonded_neighbor[1] = 15; // N3--C4
  AA[p].atom[15].nbonded_neighbors = 3; // C4
  AA[p].atom[15].bonded_neighbor[0] = 14; // C4--N3
  AA[p].atom[15].bonded_neighbor[1] = 16; // C4--C5
  AA[p].atom[15].bonded_neighbor[2] = 21; // C4--N9
  AA[p].atom[16].nbonded_neighbors = 3; // C5
  AA[p].atom[16].bonded_neighbor[0] = 15; // C5--C4
  AA[p].atom[16].bonded_neighbor[1] = 17; // C5--C6
  AA[p].atom[16].bonded_neighbor[2] = 19; // C5--N7
  AA[p].atom[17].nbonded_neighbors = 3; // C6
  AA[p].atom[17].bonded_neighbor[0] = 11; // C6--N1
  AA[p].atom[17].bonded_neighbor[1] = 16; // C6--C5
  AA[p].atom[17].bonded_neighbor[2] = 18; // C6--O6
  AA[p].atom[18].nbonded_neighbors = 1; // O6
  AA[p].atom[18].bonded_neighbor[0] = 17; // O6--C6
  AA[p].atom[19].nbonded_neighbors = 2; // N7
  AA[p].atom[19].bonded_neighbor[0] = 16; // N7--C5
  AA[p].atom[19].bonded_neighbor[1] = 20; // N7--C8
  AA[p].atom[20].nbonded_neighbors = 3; // C8
  AA[p].atom[20].bonded_neighbor[0] = 19; // C8--N7
  AA[p].atom[20].bonded_neighbor[1] = 21; // C8--N9
  AA[p].atom[20].bonded_neighbor[2] = 32; // C8--H8
  AA[p].atom[21].nbonded_neighbors = 3; // N9
  AA[p].atom[21].bonded_neighbor[0] = 10; // N9--C1*
  AA[p].atom[21].bonded_neighbor[1] = 15; // N9--C4
  AA[p].atom[21].bonded_neighbor[2] = 20; // N9--C8
  AA[p].atom[22].nbonded_neighbors = 1; //1H5*
  AA[p].atom[22].bonded_neighbor[0] = 4; //1H5*--C5*
  AA[p].atom[23].nbonded_neighbors = 1; //2H5*
  AA[p].atom[23].bonded_neighbor[0] = 4; //2H5*--C5*
  AA[p].atom[24].nbonded_neighbors = 1; // H4*
  AA[p].atom[24].bonded_neighbor[0] = 5; // H4*--C4*
  AA[p].atom[25].nbonded_neighbors = 1; // H3*
  AA[p].atom[25].bonded_neighbor[0] = 7; // H3*--C3*
  AA[p].atom[26].nbonded_neighbors = 1; //1H2*
  AA[p].atom[26].bonded_neighbor[0] = 9; //1H2*--C2*
  AA[p].atom[27].nbonded_neighbors = 1; //2H2*
  AA[p].atom[27].bonded_neighbor[0] = 9; //2H2*--C2*
  AA[p].atom[28].nbonded_neighbors = 1; // H1*
  AA[p].atom[28].bonded_neighbor[0] = 10; // H1*--C1*
  AA[p].atom[29].nbonded_neighbors = 1; // H1
  AA[p].atom[29].bonded_neighbor[0] = 11; // H1--N1
  AA[p].atom[30].nbonded_neighbors = 1; //1H2
  AA[p].atom[30].bonded_neighbor[0] = 13; //1H2--N2
  AA[p].atom[31].nbonded_neighbors = 1; //2H2
  AA[p].atom[31].bonded_neighbor[0] = 13; //2H2--N2
  AA[p].atom[32].nbonded_neighbors = 1; // H8
  AA[p].atom[32].bonded_neighbor[0] = 20; // H8--C8

  for ( k = 0; k < 3; k++ )
    AA[p].atom[k].hastemplate = true;

  for ( k = 3; k < 22; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 22; k < 33; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template atoms used for placing atoms GUA

  //bk   template for building P
  AA[p].atom[0].ta[0] = 3; //   O5*
  AA[p].atom[0].ta[1] = 4; //   C5*
  AA[p].atom[0].ta[2] = 5; //   C4*

  //bk   template for building O1P
  AA[p].atom[1].ta[0] = 0; //   P
  AA[p].atom[1].ta[1] = 4; //   C5*
  AA[p].atom[1].ta[2] = 3; //   O5*

  //bk   template for building  O2P
  AA[p].atom[2].ta[0] = 0; //   P
  AA[p].atom[2].ta[1] = 1; //   O1P
  AA[p].atom[2].ta[2] = 3; //   O5*

  //bk   template for building 1H5*
  AA[p].atom[22].ta[0] = 4; //   C5*
  AA[p].atom[22].ta[1] = 3; //   O5*
  AA[p].atom[22].ta[2] = 5; //   C4*

  //bk   template for building 2H5*
  AA[p].atom[23].ta[0] = 4; //   C5*
  AA[p].atom[23].ta[1] = 3; //   O5*
  AA[p].atom[23].ta[2] = 5; //   C4*

  //bk   template for building H4*
  AA[p].atom[24].ta[0] = 5; //   C4*
  AA[p].atom[24].ta[1] = 6; //   O4*
  AA[p].atom[24].ta[2] = 4; //   C5*

  //bk   template for building H3*
  AA[p].atom[25].ta[0] = 7; //   C3*
  AA[p].atom[25].ta[1] = 9; //   C2*
  AA[p].atom[25].ta[2] = 5; //   C4*

  //bk   template for building 1H2*
  AA[p].atom[26].ta[0] = 9; //   C2*
  AA[p].atom[26].ta[1] = 7; //   C3*
  AA[p].atom[26].ta[2] = 10; //   C1*

  //bk   template for building 2H2*
  AA[p].atom[27].ta[0] = 9; //   C2*
  AA[p].atom[27].ta[1] = 7; //   C3*
  AA[p].atom[27].ta[2] = 10; //   C1*

  //bk   template for building H1*
  AA[p].atom[28].ta[0] = 10; //   C1*
  AA[p].atom[28].ta[1] = 9; //   C2*
  AA[p].atom[28].ta[2] = 6; //   O4*

  //bk   template for building H1
  AA[p].atom[29].ta[0] = 11; //   N1
  AA[p].atom[29].ta[1] = 12; //   C2
  AA[p].atom[29].ta[2] = 17; //   C6

  //bk   template for building 1H2
  AA[p].atom[30].ta[0] = 13; //   N2
  AA[p].atom[30].ta[1] = 12; //   C2
  AA[p].atom[30].ta[2] = 14; //   N3

  //bk   template for building 2H2
  AA[p].atom[31].ta[0] = 13; //   N2
  AA[p].atom[31].ta[1] = 12; //   C2
  AA[p].atom[31].ta[2] = 14; //   N3

  //bk   template for building H8
  AA[p].atom[32].ta[0] = 20; //   C8
  AA[p].atom[32].ta[1] = 19; //   N7
  AA[p].atom[32].ta[2] = 21; //   N9


  //bk   template coordinates for the amino acid
  //jjh new template from Kevin Karplus

  AA[p].atom[0].icoor[0] = 23.2380; //   P
  AA[p].atom[0].icoor[1] = 27.4740; //   P
  AA[p].atom[0].icoor[2] = -0.3190; //   P
  AA[p].atom[1].icoor[0] = 24.1140; //   O1P
  AA[p].atom[1].icoor[1] = 28.4410; //   O1P
  AA[p].atom[1].icoor[2] = -1.0010; //   O1P
  AA[p].atom[2].icoor[0] = 23.8070; //   O2P
  AA[p].atom[2].icoor[1] = 26.2020; //   O2P
  AA[p].atom[2].icoor[2] = 0.1760; //   O2P
  AA[p].atom[3].icoor[0] = 22.5050; //   O5*
  AA[p].atom[3].icoor[1] = 28.2470; //   O5*
  AA[p].atom[3].icoor[2] = 0.8590; //   O5*
  AA[p].atom[4].icoor[0] = 21.7980; //   C5*
  AA[p].atom[4].icoor[1] = 29.4430; //   C5*
  AA[p].atom[4].icoor[2] = 0.5720; //   C5*
  AA[p].atom[5].icoor[0] = 20.8380; //   C4*
  AA[p].atom[5].icoor[1] = 29.7890; //   C4*
  AA[p].atom[5].icoor[2] = 1.6850; //   C4*
  AA[p].atom[6].icoor[0] = 19.6900; //   O4*
  AA[p].atom[6].icoor[1] = 28.9030; //   O4*
  AA[p].atom[6].icoor[2] = 1.7070; //   O4*
  AA[p].atom[7].icoor[0] = 21.3970; //   C3*
  AA[p].atom[7].icoor[1] = 29.8070; //   C3*
  AA[p].atom[7].icoor[2] = 3.1100; //   C3*
  AA[p].atom[8].icoor[0] = 20.9400; //   O3*
  AA[p].atom[8].icoor[1] = 31.0170; //   O3*
  AA[p].atom[8].icoor[2] = 3.7240; //   O3*
  AA[p].atom[9].icoor[0] = 20.7290; //   C2*
  AA[p].atom[9].icoor[1] = 28.6020; //   C2*
  AA[p].atom[9].icoor[2] = 3.7520; //   C2*
  AA[p].atom[10].icoor[0] = 19.3850; //   C1*
  AA[p].atom[10].icoor[1] = 28.6210; //   C1*
  AA[p].atom[10].icoor[2] = 3.0560; //   C1*
  AA[p].atom[11].icoor[0] = 14.8590; //   N1
  AA[p].atom[11].icoor[1] = 26.5120; //   N1
  AA[p].atom[11].icoor[2] = 3.4000; //   N1
  AA[p].atom[12].icoor[0] = 15.1620; //   C2
  AA[p].atom[12].icoor[1] = 27.8120; //   C2
  AA[p].atom[12].icoor[2] = 3.4420; //   C2
  AA[p].atom[13].icoor[0] = 14.1390; //   N2
  AA[p].atom[13].icoor[1] = 28.6350; //   N2
  AA[p].atom[13].icoor[2] = 3.5600; //   N2
  AA[p].atom[14].icoor[0] = 16.3680; //   N3
  AA[p].atom[14].icoor[1] = 28.2550; //   N3
  AA[p].atom[14].icoor[2] = 3.3740; //   N3
  AA[p].atom[15].icoor[0] = 17.2210; //   C4
  AA[p].atom[15].icoor[1] = 27.2500; //   C4
  AA[p].atom[15].icoor[2] = 3.2160; //   C4
  AA[p].atom[16].icoor[0] = 16.9970; //   C5
  AA[p].atom[16].icoor[1] = 25.9070; //   C5
  AA[p].atom[16].icoor[2] = 3.1720; //   C5
  AA[p].atom[17].icoor[0] = 15.6950; //   C6
  AA[p].atom[17].icoor[1] = 25.4530; //   C6
  AA[p].atom[17].icoor[2] = 3.2810; //   C6
  AA[p].atom[18].icoor[0] = 15.2710; //   O6
  AA[p].atom[18].icoor[1] = 24.2950; //   O6
  AA[p].atom[18].icoor[2] = 3.2680; //   O6
  AA[p].atom[19].icoor[0] = 18.1780; //   N7
  AA[p].atom[19].icoor[1] = 25.1930; //   N7
  AA[p].atom[19].icoor[2] = 2.9990; //   N7
  AA[p].atom[20].icoor[0] = 19.1060; //   C8
  AA[p].atom[20].icoor[1] = 26.1210; //   C8
  AA[p].atom[20].icoor[2] = 2.9290; //   C8
  AA[p].atom[21].icoor[0] = 18.5990; //   N9
  AA[p].atom[21].icoor[1] = 27.3810; //   N9
  AA[p].atom[21].icoor[2] = 3.0780; //   N9
  AA[p].atom[22].icoor[0] = 21.2320; //  1H5*
  AA[p].atom[22].icoor[1] = 29.3070; //  1H5*
  AA[p].atom[22].icoor[2] = -0.3740; //  1H5*
  AA[p].atom[23].icoor[0] = 22.5280; //  2H5*
  AA[p].atom[23].icoor[1] = 30.2700; //  2H5*
  AA[p].atom[23].icoor[2] = 0.4450; //  2H5*
  AA[p].atom[24].icoor[0] = 20.4600; //   H4*
  AA[p].atom[24].icoor[1] = 30.8250; //   H4*
  AA[p].atom[24].icoor[2] = 1.5460; //   H4*
  AA[p].atom[25].icoor[0] = 22.5030; //   H3*
  AA[p].atom[25].icoor[1] = 29.7000; //   H3*
  AA[p].atom[25].icoor[2] = 3.0900; //   H3*
  AA[p].atom[26].icoor[0] = 21.2870; //  1H2*
  AA[p].atom[26].icoor[1] = 27.6650; //  1H2*
  AA[p].atom[26].icoor[2] = 3.5460; //  1H2*
  AA[p].atom[27].icoor[0] = 20.6250; //  2H2*
  AA[p].atom[27].icoor[1] = 28.7260; //  2H2*
  AA[p].atom[27].icoor[2] = 4.8510; //  2H2*
  AA[p].atom[28].icoor[0] = 18.7280; //   H1*
  AA[p].atom[28].icoor[1] = 29.4180; //   H1*
  AA[p].atom[28].icoor[2] = 3.4610; //   H1*
  AA[p].atom[29].icoor[0] = 13.8880; //   H1
  AA[p].atom[29].icoor[1] = 26.2710; //   H1
  AA[p].atom[29].icoor[2] = 3.5040; //   H1
  AA[p].atom[30].icoor[0] = 13.1860; //  1H2
  AA[p].atom[30].icoor[1] = 28.3260; //  1H2
  AA[p].atom[30].icoor[2] = 3.6200; //  1H2
  AA[p].atom[31].icoor[0] = 14.3470; //  2H2
  AA[p].atom[31].icoor[1] = 29.6020; //  2H2
  AA[p].atom[31].icoor[2] = 3.6190; //  2H2
  AA[p].atom[32].icoor[0] = 20.1740; //   H8
  AA[p].atom[32].icoor[1] = 25.9530; //   H8
  AA[p].atom[32].icoor[2] = 2.7630; //   H8


  //FIN DE GUA


  //ADE
//  p = ADE;
  p = DADE;

//  strcpy( AA[p].aa_name3, "  A" );
  strcpy( AA[p].aa_name3, " DA" );
//  AA[p].aa_name1 = 'X';
  AA[p].aa_name1 = 'A';

  AA[p].mass = 0;

  AA[p].is_protein = false;
  AA[p].is_RNA = false;
  AA[p].is_DNA = true;

  AA[p].aa_is_polar = true;
  AA[p].aa_is_nonpolar = false;
  AA[p].aa_is_aromatic = true;
  AA[p].aa_is_charged = true;

  AA[p].natoms = 32;
  AA[p].nheavyatoms = 21;
  AA[p].nchi = 0;

  AA[p].chi_atoms = NULL;
  AA[p].atom = ( t_aa_atom_p * ) malloc( AA[p].natoms * sizeof( t_aa_atom_p ) );


  strcpy( AA[p].atom[0].atom_name, " P  " );
//  strcpy( AA[p].atom[1].atom_name, " O1P" );
//  strcpy( AA[p].atom[2].atom_name, " O2P" );

  strcpy( AA[p].atom[1].atom_name, " OP1" );
  strcpy( AA[p].atom[2].atom_name, " OP2" );

//  strcpy( AA[p].atom[3].atom_name, " O5*" );
//  strcpy( AA[p].atom[4].atom_name, " C5*" );
//  strcpy( AA[p].atom[5].atom_name, " C4*" );
//  strcpy( AA[p].atom[6].atom_name, " O4*" );
//  strcpy( AA[p].atom[7].atom_name, " C3*" );
//  strcpy( AA[p].atom[8].atom_name, " O3*" );
//  strcpy( AA[p].atom[9].atom_name, " C2*" );
//  strcpy( AA[p].atom[10].atom_name, " C1*" );

  strcpy( AA[p].atom[3].atom_name, " O5'" );
  strcpy( AA[p].atom[4].atom_name, " C5'" );
  strcpy( AA[p].atom[5].atom_name, " C4'" );
  strcpy( AA[p].atom[6].atom_name, " O4'" );
  strcpy( AA[p].atom[7].atom_name, " C3'" );
  strcpy( AA[p].atom[8].atom_name, " O3'" );
  strcpy( AA[p].atom[9].atom_name, " C2'" );
  strcpy( AA[p].atom[10].atom_name, " C1'" );

  strcpy( AA[p].atom[11].atom_name, " N1 " );
  strcpy( AA[p].atom[12].atom_name, " C2 " );
  strcpy( AA[p].atom[13].atom_name, " N3 " );
  strcpy( AA[p].atom[14].atom_name, " C4 " );
  strcpy( AA[p].atom[15].atom_name, " C5 " );
  strcpy( AA[p].atom[16].atom_name, " C6 " );
  strcpy( AA[p].atom[17].atom_name, " N6 " );
  strcpy( AA[p].atom[18].atom_name, " N7 " );
  strcpy( AA[p].atom[19].atom_name, " C8 " );
  strcpy( AA[p].atom[20].atom_name, " N9 " );

//  strcpy( AA[p].atom[21].atom_name, "1H5*" );
//  strcpy( AA[p].atom[22].atom_name, "2H5*" );
//  strcpy( AA[p].atom[23].atom_name, " H4*" );
//  strcpy( AA[p].atom[24].atom_name, " H3*" );
//  strcpy( AA[p].atom[25].atom_name, "1H2*" );
//  strcpy( AA[p].atom[26].atom_name, "2H2*" );
//  strcpy( AA[p].atom[27].atom_name, " H1*" );

  strcpy( AA[p].atom[21].atom_name, "1H5'" );
  strcpy( AA[p].atom[22].atom_name, "2H5'" );
  strcpy( AA[p].atom[23].atom_name, " H4'" );
  strcpy( AA[p].atom[24].atom_name, " H3'" );
  strcpy( AA[p].atom[25].atom_name, "1H2'" );
  strcpy( AA[p].atom[26].atom_name, "2H2'" );
  strcpy( AA[p].atom[27].atom_name, " H1'" );

  strcpy( AA[p].atom[28].atom_name, " H2 " );
  strcpy( AA[p].atom[29].atom_name, "1H6 " );
  strcpy( AA[p].atom[30].atom_name, "2H6 " );
  strcpy( AA[p].atom[31].atom_name, " H8 " );


  AA[p].atom[0].charge = 1.500; // P
  AA[p].atom[1].charge = -.780; // O1P
  AA[p].atom[2].charge = -.780; // O2P
  AA[p].atom[3].charge = -.570; // O5*
  AA[p].atom[4].charge = -.080; // C5*
  AA[p].atom[5].charge =.160; // C4*
  AA[p].atom[6].charge = -.500; // O4*
  AA[p].atom[7].charge =.010; // C3*
  AA[p].atom[8].charge = -.570; // O3*
  AA[p].atom[9].charge = -.180; // C2*
  AA[p].atom[10].charge =.160; // C1*
  AA[p].atom[11].charge = -.740; // N1
  AA[p].atom[12].charge =.500; // C2
  AA[p].atom[13].charge = -.750; // N3
  AA[p].atom[14].charge =.430; // C4
  AA[p].atom[15].charge =.280; // C5
  AA[p].atom[16].charge =.460; // C6
  AA[p].atom[17].charge = -.770; // N6
  AA[p].atom[18].charge = -.710; // N7
  AA[p].atom[19].charge =.340; // C8
  AA[p].atom[20].charge = -.050; // N9
  AA[p].atom[21].charge =.090; // 1H5*
  AA[p].atom[22].charge =.090; // 2H5*
  AA[p].atom[23].charge =.090; // H4*
  AA[p].atom[24].charge =.090; // H3*
  AA[p].atom[25].charge =.090; // 1H2*
  AA[p].atom[26].charge =.090; // 2H2*
  AA[p].atom[27].charge =.090; // H1*
  AA[p].atom[28].charge =.130; // H2
  AA[p].atom[29].charge =.380; // 1H6
  AA[p].atom[30].charge =.380; // 2H6
  AA[p].atom[31].charge =.120; // H8

  switch(opt){
  case Sybil:
  	  AA[p].atom[0].fullatom_type = 12; // Phos   P
  	  AA[p].atom[1].fullatom_type = 10; // OOC    O1P
  	  AA[p].atom[2].fullatom_type = 10; // OOC    O2P
  	  AA[p].atom[3].fullatom_type = 9; // ONH2   O5*
  	  AA[p].atom[4].fullatom_type = 2; // CH2    C5*
  	  AA[p].atom[5].fullatom_type = 2; // CH1    C4*
  	  AA[p].atom[6].fullatom_type = 9; // OH     O4*
  	  AA[p].atom[7].fullatom_type = 2; // CH1    C3*
  	  AA[p].atom[8].fullatom_type = 9; // ONH2   O3*
  	  AA[p].atom[9].fullatom_type = 2; // CH2    C2*
  	  AA[p].atom[10].fullatom_type = 2; // CH1    C1*
  	  AA[p].atom[11].fullatom_type = 5; // Nhis   N1
  	  AA[p].atom[12].fullatom_type = 3; // aroC   C2
  	  AA[p].atom[13].fullatom_type = 5; // Nhis   N3
  	  AA[p].atom[14].fullatom_type = 3; // aroC   C4
  	  AA[p].atom[15].fullatom_type = 3; // aroC   C5
  	  AA[p].atom[16].fullatom_type = 3; // aroC   C6
  	  AA[p].atom[17].fullatom_type = 5; // NH2O   N6
  	  AA[p].atom[18].fullatom_type = 5; // Nhis   N7
  	  AA[p].atom[19].fullatom_type = 3; // aroC   C8
  	  AA[p].atom[20].fullatom_type = 5; // Nhis   N9
  	  AA[p].atom[21].fullatom_type = 14; // Hapo  1H5*
  	  AA[p].atom[22].fullatom_type = 14; // Hapo  2H5*
  	  AA[p].atom[23].fullatom_type = 14; // Hapo   H4*
  	  AA[p].atom[24].fullatom_type = 14; // Hapo   H3*
  	  AA[p].atom[25].fullatom_type = 14; // Hapo  1H2*
  	  AA[p].atom[26].fullatom_type = 14; // Hapo  2H2*
  	  AA[p].atom[27].fullatom_type = 14; // Hapo   H1*
  	  AA[p].atom[28].fullatom_type = 15; // Haro   H2
  	  AA[p].atom[29].fullatom_type = 13; // Hpol  1H6
  	  AA[p].atom[30].fullatom_type = 13; // Hpol  2H6
  	  AA[p].atom[31].fullatom_type = 15; // Haro   H8
  	  break;
  default:
	  AA[p].atom[0].fullatom_type = 21; // Phos   P
	  AA[p].atom[1].fullatom_type = 15; // OOC    O1P
	  AA[p].atom[2].fullatom_type = 15; // OOC    O2P
	  AA[p].atom[3].fullatom_type = 14; // ONH2   O5*
	  AA[p].atom[4].fullatom_type = 4; // CH2    C5*
	  AA[p].atom[5].fullatom_type = 3; // CH1    C4*
	  AA[p].atom[6].fullatom_type = 13; // OH     O4*
	  AA[p].atom[7].fullatom_type = 3; // CH1    C3*
	  AA[p].atom[8].fullatom_type = 14; // ONH2   O3*
	  AA[p].atom[9].fullatom_type = 4; // CH2    C2*
	  AA[p].atom[10].fullatom_type = 3; // CH1    C1*
	  AA[p].atom[11].fullatom_type = 8; // Nhis   N1
	  AA[p].atom[12].fullatom_type = 6; // aroC   C2
	  AA[p].atom[13].fullatom_type = 8; // Nhis   N3
	  AA[p].atom[14].fullatom_type = 6; // aroC   C4
	  AA[p].atom[15].fullatom_type = 6; // aroC   C5
	  AA[p].atom[16].fullatom_type = 6; // aroC   C6
	  AA[p].atom[17].fullatom_type = 9; // NH2O   N6
	  AA[p].atom[18].fullatom_type = 8; // Nhis   N7
	  AA[p].atom[19].fullatom_type = 6; // aroC   C8
	  AA[p].atom[20].fullatom_type = 8; // Nhis   N9
	  AA[p].atom[21].fullatom_type = 23; // Hapo  1H5*
	  AA[p].atom[22].fullatom_type = 23; // Hapo  2H5*
	  AA[p].atom[23].fullatom_type = 23; // Hapo   H4*
	  AA[p].atom[24].fullatom_type = 23; // Hapo   H3*
	  AA[p].atom[25].fullatom_type = 23; // Hapo  1H2*
	  AA[p].atom[26].fullatom_type = 23; // Hapo  2H2*
	  AA[p].atom[27].fullatom_type = 23; // Hapo   H1*
	  AA[p].atom[28].fullatom_type = 24; // Haro   H2
	  AA[p].atom[29].fullatom_type = 22; // Hpol  1H6
	  AA[p].atom[30].fullatom_type = 22; // Hpol  2H6
	  AA[p].atom[31].fullatom_type = 24; // Haro   H8
	  break;
  }

  //jjh intra residue bonding

  AA[p].atom[0].nbonded_neighbors = 3; // P
  AA[p].atom[0].bonded_neighbor[0] = 1; // P--O1P
  AA[p].atom[0].bonded_neighbor[1] = 2; // P--O2P
  AA[p].atom[0].bonded_neighbor[2] = 3; // P--O5*
  AA[p].atom[1].nbonded_neighbors = 1; // O1P
  AA[p].atom[1].bonded_neighbor[0] = 0; // O1P--P
  AA[p].atom[2].nbonded_neighbors = 1; // O2P
  AA[p].atom[2].bonded_neighbor[0] = 0; // O2P--P
  AA[p].atom[3].nbonded_neighbors = 2; // O5*
  AA[p].atom[3].bonded_neighbor[0] = 0; // O5*--P
  AA[p].atom[3].bonded_neighbor[1] = 4; // O5*--C5*
  AA[p].atom[4].nbonded_neighbors = 4; // C5*
  AA[p].atom[4].bonded_neighbor[0] = 3; // C5*--O5*
  AA[p].atom[4].bonded_neighbor[1] = 5; // C5*--C4*
  AA[p].atom[4].bonded_neighbor[2] = 21; // C5*--1H5*
  AA[p].atom[4].bonded_neighbor[3] = 22; // C5*--2H5*
  AA[p].atom[5].nbonded_neighbors = 4; // C4*
  AA[p].atom[5].bonded_neighbor[0] = 4; // C4*--C5*
  AA[p].atom[5].bonded_neighbor[1] = 6; // C4*--O4*
  AA[p].atom[5].bonded_neighbor[2] = 7; // C4*--C3*
  AA[p].atom[5].bonded_neighbor[3] = 23; // C4*--H4*
  AA[p].atom[6].nbonded_neighbors = 2; // O4*
  AA[p].atom[6].bonded_neighbor[0] = 5; // O4*--C4*
  AA[p].atom[6].bonded_neighbor[1] = 10; // O4*--C1*
  AA[p].atom[7].nbonded_neighbors = 4; // C3*
  AA[p].atom[7].bonded_neighbor[0] = 5; // C3*--C4*
  AA[p].atom[7].bonded_neighbor[1] = 8; // C3*--O3*
  AA[p].atom[7].bonded_neighbor[2] = 9; // C3*--C2*
  AA[p].atom[7].bonded_neighbor[3] = 24; // C3*--H3*
  AA[p].atom[8].nbonded_neighbors = 1; // O3*
  AA[p].atom[8].bonded_neighbor[0] = 7; // O3*--C3*
  AA[p].atom[9].nbonded_neighbors = 4; // C2*
  AA[p].atom[9].bonded_neighbor[0] = 7; // C2*--C3*
  AA[p].atom[9].bonded_neighbor[1] = 10; // C2*--C1*
  AA[p].atom[9].bonded_neighbor[2] = 25; // C2*--1H2*
  AA[p].atom[9].bonded_neighbor[3] = 26; // C2*--2H2*
  AA[p].atom[10].nbonded_neighbors = 4; // C1*
  AA[p].atom[10].bonded_neighbor[0] = 6; // C1*--O4*
  AA[p].atom[10].bonded_neighbor[1] = 9; // C1*--C2*
  AA[p].atom[10].bonded_neighbor[2] = 20; // C1*--N9
  AA[p].atom[10].bonded_neighbor[3] = 27; // C1*--H1*
  AA[p].atom[11].nbonded_neighbors = 2; // N1
  AA[p].atom[11].bonded_neighbor[0] = 12; // N1--C2
  AA[p].atom[11].bonded_neighbor[1] = 16; // N1--C6
  AA[p].atom[12].nbonded_neighbors = 3; // C2
  AA[p].atom[12].bonded_neighbor[0] = 11; // C2--N1
  AA[p].atom[12].bonded_neighbor[1] = 13; // C2--N3
  AA[p].atom[12].bonded_neighbor[2] = 28; // C2--H2
  AA[p].atom[13].nbonded_neighbors = 2; // N3
  AA[p].atom[13].bonded_neighbor[0] = 12; // N3--C2
  AA[p].atom[13].bonded_neighbor[1] = 14; // N3--C4
  AA[p].atom[14].nbonded_neighbors = 3; // C4
  AA[p].atom[14].bonded_neighbor[0] = 13; // C4--N3
  AA[p].atom[14].bonded_neighbor[1] = 15; // C4--C5
  AA[p].atom[14].bonded_neighbor[2] = 20; // C4--N9
  AA[p].atom[15].nbonded_neighbors = 3; // C5
  AA[p].atom[15].bonded_neighbor[0] = 14; // C5--C4
  AA[p].atom[15].bonded_neighbor[1] = 16; // C5--C6
  AA[p].atom[15].bonded_neighbor[2] = 18; // C5--N7
  AA[p].atom[16].nbonded_neighbors = 3; // C6
  AA[p].atom[16].bonded_neighbor[0] = 11; // C6--N1
  AA[p].atom[16].bonded_neighbor[1] = 15; // C6--C5
  AA[p].atom[16].bonded_neighbor[2] = 17; // C6--N6
  AA[p].atom[17].nbonded_neighbors = 3; // N6
  AA[p].atom[17].bonded_neighbor[0] = 16; // N6--C6
  AA[p].atom[17].bonded_neighbor[1] = 29; // N6--1H6
  AA[p].atom[17].bonded_neighbor[2] = 30; // N6--2H6
  AA[p].atom[18].nbonded_neighbors = 2; // N7
  AA[p].atom[18].bonded_neighbor[0] = 15; // N7--C5
  AA[p].atom[18].bonded_neighbor[1] = 19; // N7--C8
  AA[p].atom[19].nbonded_neighbors = 3; // C8
  AA[p].atom[19].bonded_neighbor[0] = 18; // C8--N7
  AA[p].atom[19].bonded_neighbor[1] = 20; // C8--N9
  AA[p].atom[19].bonded_neighbor[2] = 31; // C8--H8
  AA[p].atom[20].nbonded_neighbors = 3; // N9
  AA[p].atom[20].bonded_neighbor[0] = 10; // N9--C1*
  AA[p].atom[20].bonded_neighbor[1] = 14; // N9--C4
  AA[p].atom[20].bonded_neighbor[2] = 19; // N9--C8
  AA[p].atom[21].nbonded_neighbors = 1; //1H5*
  AA[p].atom[21].bonded_neighbor[0] = 4; //1H5*--C5*
  AA[p].atom[22].nbonded_neighbors = 1; //2H5*
  AA[p].atom[22].bonded_neighbor[0] = 4; //2H5*--C5*
  AA[p].atom[23].nbonded_neighbors = 1; // H4*
  AA[p].atom[23].bonded_neighbor[0] = 5; // H4*--C4*
  AA[p].atom[24].nbonded_neighbors = 1; // H3*
  AA[p].atom[24].bonded_neighbor[0] = 7; // H3*--C3*
  AA[p].atom[25].nbonded_neighbors = 1; //1H2*
  AA[p].atom[25].bonded_neighbor[0] = 9; //1H2*--C2*
  AA[p].atom[26].nbonded_neighbors = 1; //2H2*
  AA[p].atom[26].bonded_neighbor[0] = 9; //2H2*--C2*
  AA[p].atom[27].nbonded_neighbors = 1; // H1*
  AA[p].atom[27].bonded_neighbor[0] = 10; // H1*--C1*
  AA[p].atom[28].nbonded_neighbors = 1; // H2
  AA[p].atom[28].bonded_neighbor[0] = 12; // H2--C2
  AA[p].atom[29].nbonded_neighbors = 1; //1H6
  AA[p].atom[29].bonded_neighbor[0] = 17; //1H6--N6
  AA[p].atom[30].nbonded_neighbors = 1; //2H6
  AA[p].atom[30].bonded_neighbor[0] = 17; //2H6--N6
  AA[p].atom[31].nbonded_neighbors = 1; // H8
  AA[p].atom[31].bonded_neighbor[0] = 19; // H8--C8

  for ( k = 0; k < 3; k++ )
    AA[p].atom[k].hastemplate = true;

  for ( k = 3; k < 21; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 21; k < 32; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template atoms used for placing atoms ADE

  //bk   template for building P
  AA[p].atom[0].ta[0] = 3; //   O5*
  AA[p].atom[0].ta[1] = 4; //   C5*
  AA[p].atom[0].ta[2] = 5; //   C4*

  //bk   template for building O1P
  AA[p].atom[1].ta[0] = 0; //   P
  AA[p].atom[1].ta[1] = 4; //   C5*
  AA[p].atom[1].ta[2] = 3; //   O5*

  //bk   template for building  O2P
  AA[p].atom[2].ta[0] = 0; //   P
  AA[p].atom[2].ta[1] = 1; //   O1P
  AA[p].atom[2].ta[2] = 3; //   O5*

  //bk   template for building 1H5*
  AA[p].atom[21].ta[0] = 4; //   C5*
  AA[p].atom[21].ta[1] = 3; //   O5*
  AA[p].atom[21].ta[2] = 5; //   C4*

  //bk   template for building 2H5*
  AA[p].atom[22].ta[0] = 4; //   C5*
  AA[p].atom[22].ta[1] = 3; //   O5*
  AA[p].atom[22].ta[2] = 5; //   C4*

  //bk   template for building H4*
  AA[p].atom[23].ta[0] = 5; //   C4*
  AA[p].atom[23].ta[1] = 6; //   O4*
  AA[p].atom[23].ta[2] = 4; //   C5*

  //bk   template for building H3*
  AA[p].atom[24].ta[0] = 7; //   C3*
  AA[p].atom[24].ta[1] = 9; //   C2*
  AA[p].atom[24].ta[2] = 5; //   C4*

  //bk   template for building 1H2*
  AA[p].atom[25].ta[0] = 9; //   C2*
  AA[p].atom[25].ta[1] = 7; //   C3*
  AA[p].atom[25].ta[2] = 10; //   C1*

  //bk   template for building 2H2*
  AA[p].atom[26].ta[0] = 9; //   C2*
  AA[p].atom[26].ta[1] = 7; //   C3*
  AA[p].atom[26].ta[2] = 10; //   C1*

  //bk   template for building H1*
  AA[p].atom[27].ta[0] = 10; //   C1*
  AA[p].atom[27].ta[1] = 9; //   C2*
  AA[p].atom[27].ta[2] = 6; //   O4*

  //bk   template for building H2
  AA[p].atom[28].ta[0] = 12; //   C2
  AA[p].atom[28].ta[1] = 11; //   N1
  AA[p].atom[28].ta[2] = 13; //   N3

  //bk   template for building 1H6
  AA[p].atom[29].ta[0] = 17; //   N6
  AA[p].atom[29].ta[1] = 16; //   C6
  AA[p].atom[29].ta[2] = 11; //   N1

  //bk   template for building 2H6
  AA[p].atom[30].ta[0] = 17; //   N6
  AA[p].atom[30].ta[1] = 16; //   C6
  AA[p].atom[30].ta[2] = 11; //   N1

  //bk   template for building H8
  AA[p].atom[31].ta[0] = 19; //   C8
  AA[p].atom[31].ta[1] = 18; //   N7
  AA[p].atom[31].ta[2] = 20; //   N9


  //bk   template coordinates for the amino acid
  //jjh new template from Kevin Karplus

  AA[p].atom[0].icoor[0] = 11.0150; //   P
  AA[p].atom[0].icoor[1] = 33.6920; //   P
  AA[p].atom[0].icoor[2] = 12.8880; //   P
  AA[p].atom[1].icoor[0] = 10.1500; //   O1P
  AA[p].atom[1].icoor[1] = 34.8960; //   O1P
  AA[p].atom[1].icoor[2] = 12.8560; //   O1P
  AA[p].atom[2].icoor[0] = 12.2240; //   O2P
  AA[p].atom[2].icoor[1] = 33.6790; //   O2P
  AA[p].atom[2].icoor[2] = 13.7470; //   O2P
  AA[p].atom[3].icoor[0] = 10.1110; //   O5*
  AA[p].atom[3].icoor[1] = 32.4520; //   O5*
  AA[p].atom[3].icoor[2] = 13.3010; //   O5*
  AA[p].atom[4].icoor[0] = 9.0260; //   C5*
  AA[p].atom[4].icoor[1] = 32.0420; //   C5*
  AA[p].atom[4].icoor[2] = 12.4730; //   C5*
  AA[p].atom[5].icoor[0] = 8.3730; //   C4*
  AA[p].atom[5].icoor[1] = 30.8110; //   C4*
  AA[p].atom[5].icoor[2] = 13.0510; //   C4*
  AA[p].atom[6].icoor[0] = 9.2390; //   O4*
  AA[p].atom[6].icoor[1] = 29.6680; //   O4*
  AA[p].atom[6].icoor[2] = 12.8620; //   O4*
  AA[p].atom[7].icoor[0] = 8.1170; //   C3*
  AA[p].atom[7].icoor[1] = 30.9090; //   C3*
  AA[p].atom[7].icoor[2] = 14.5530; //   C3*
  AA[p].atom[8].icoor[0] = 6.8570; //   O3*
  AA[p].atom[8].icoor[1] = 30.3240; //   O3*
  AA[p].atom[8].icoor[2] = 14.8830; //   O3*
  AA[p].atom[9].icoor[0] = 9.2600; //   C2*
  AA[p].atom[9].icoor[1] = 30.1180; //   C2*
  AA[p].atom[9].icoor[2] = 15.1610; //   C2*
  AA[p].atom[10].icoor[0] = 9.5280; //   C1*
  AA[p].atom[10].icoor[1] = 29.0570; //   C1*
  AA[p].atom[10].icoor[2] = 14.1090; //   C1*
  AA[p].atom[11].icoor[0] = 12.6670; //   N1
  AA[p].atom[11].icoor[1] = 25.0780; //   N1
  AA[p].atom[11].icoor[2] = 13.9220; //   N1
  AA[p].atom[12].icoor[0] = 11.3780; //   C2
  AA[p].atom[12].icoor[1] = 25.1620; //   C2
  AA[p].atom[12].icoor[2] = 14.0860; //   C2
  AA[p].atom[13].icoor[0] = 10.6170; //   N3
  AA[p].atom[13].icoor[1] = 26.2230; //   N3
  AA[p].atom[13].icoor[2] = 14.1600; //   N3
  AA[p].atom[14].icoor[0] = 11.3440; //   C4
  AA[p].atom[14].icoor[1] = 27.3070; //   C4
  AA[p].atom[14].icoor[2] = 14.0150; //   C4
  AA[p].atom[15].icoor[0] = 12.6860; //   C5
  AA[p].atom[15].icoor[1] = 27.3870; //   C5
  AA[p].atom[15].icoor[2] = 13.8430; //   C5
  AA[p].atom[16].icoor[0] = 13.3790; //   C6
  AA[p].atom[16].icoor[1] = 26.1970; //   C6
  AA[p].atom[16].icoor[2] = 13.7930; //   C6
  AA[p].atom[17].icoor[0] = 14.7040; //   N6
  AA[p].atom[17].icoor[1] = 26.1270; //   N6
  AA[p].atom[17].icoor[2] = 13.6150; //   N6
  AA[p].atom[18].icoor[0] = 13.1050; //   N7
  AA[p].atom[18].icoor[1] = 28.7080; //   N7
  AA[p].atom[18].icoor[2] = 13.7160; //   N7
  AA[p].atom[19].icoor[0] = 11.9830; //   C8
  AA[p].atom[19].icoor[1] = 29.3920; //   C8
  AA[p].atom[19].icoor[2] = 13.7990; //   C8
  AA[p].atom[20].icoor[0] = 10.8950; //   N9
  AA[p].atom[20].icoor[1] = 28.6050; //   N9
  AA[p].atom[20].icoor[2] = 14.0150; //   N9
  AA[p].atom[21].icoor[0] = 9.4080; //  1H5*
  AA[p].atom[21].icoor[1] = 31.8180; //  1H5*
  AA[p].atom[21].icoor[2] = 11.4540; //  1H5*
  AA[p].atom[22].icoor[0] = 8.2830; //  2H5*
  AA[p].atom[22].icoor[1] = 32.8660; //  2H5*
  AA[p].atom[22].icoor[2] = 12.4100; //  2H5*
  AA[p].atom[23].icoor[0] = 7.4110; //   H4*
  AA[p].atom[23].icoor[1] = 30.6070; //   H4*
  AA[p].atom[23].icoor[2] = 12.5340; //   H4*
  AA[p].atom[24].icoor[0] = 8.1410; //   H3*
  AA[p].atom[24].icoor[1] = 31.9720; //   H3*
  AA[p].atom[24].icoor[2] = 14.8830; //   H3*
  AA[p].atom[25].icoor[0] = 10.1520; //  1H2*
  AA[p].atom[25].icoor[1] = 30.7590; //  1H2*
  AA[p].atom[25].icoor[2] = 15.3240; //  1H2'
  AA[p].atom[26].icoor[0] = 8.9640; //  2H2*
  AA[p].atom[26].icoor[1] = 29.6610; //  2H2*
  AA[p].atom[26].icoor[2] = 16.1290; //  2H2*
  AA[p].atom[27].icoor[0] = 8.8640; //   H1*
  AA[p].atom[27].icoor[1] = 28.1780; //   H1*
  AA[p].atom[27].icoor[2] = 14.2460; //   H1*
  AA[p].atom[28].icoor[0] = 10.8530; //   H2
  AA[p].atom[28].icoor[1] = 24.2040; //   H2
  AA[p].atom[28].icoor[2] = 14.1880; //   H2
  AA[p].atom[29].icoor[0] = 15.1400; //  1H6
  AA[p].atom[29].icoor[1] = 25.2200; //  1H6
  AA[p].atom[29].icoor[2] = 13.5970; //  1H6
  AA[p].atom[30].icoor[0] = 15.2550; //  2H6
  AA[p].atom[30].icoor[1] = 26.9530; //  2H6
  AA[p].atom[30].icoor[2] = 13.5800; //  2H6
  AA[p].atom[31].icoor[0] = 11.8780; //   H8
  AA[p].atom[31].icoor[1] = 30.4800; //   H8
  AA[p].atom[31].icoor[2] = 13.7400; //   H8

  //FIN DE ADE


  //CYT
//  p = CYT;
  p = DCYT;

//  strcpy( AA[p].aa_name3, "  C" );
  strcpy( AA[p].aa_name3, " DC" );
//  AA[p].aa_name1 = 'X';
  AA[p].aa_name1 = 'C';

  AA[p].mass = 0;

  AA[p].is_protein = false;
  AA[p].is_RNA = false;
  AA[p].is_DNA = true;

  AA[p].aa_is_polar = true;
  AA[p].aa_is_nonpolar = false;
  AA[p].aa_is_aromatic = true;
  AA[p].aa_is_charged = true;

  AA[p].natoms = 30;
  AA[p].nheavyatoms = 19;
  AA[p].nchi = 0;

  AA[p].chi_atoms = NULL;
  AA[p].atom = ( t_aa_atom_p * ) malloc( AA[p].natoms * sizeof( t_aa_atom_p ) );


  strcpy( AA[p].atom[0].atom_name, " P  " );
//  strcpy( AA[p].atom[1].atom_name, " O1P" );
//  strcpy( AA[p].atom[2].atom_name, " O2P" );

  strcpy( AA[p].atom[1].atom_name, " OP1" );
  strcpy( AA[p].atom[2].atom_name, " OP2" );

//  strcpy( AA[p].atom[3].atom_name, " O5*" );
//  strcpy( AA[p].atom[4].atom_name, " C5*" );
//  strcpy( AA[p].atom[5].atom_name, " C4*" );
//  strcpy( AA[p].atom[6].atom_name, " O4*" );
//  strcpy( AA[p].atom[7].atom_name, " C3*" );
//  strcpy( AA[p].atom[8].atom_name, " O3*" );
//  strcpy( AA[p].atom[9].atom_name, " C2*" );
//  strcpy( AA[p].atom[10].atom_name, " C1*" );

  strcpy( AA[p].atom[3].atom_name, " O5'" );
  strcpy( AA[p].atom[4].atom_name, " C5'" );
  strcpy( AA[p].atom[5].atom_name, " C4'" );
  strcpy( AA[p].atom[6].atom_name, " O4'" );
  strcpy( AA[p].atom[7].atom_name, " C3'" );
  strcpy( AA[p].atom[8].atom_name, " O3'" );
  strcpy( AA[p].atom[9].atom_name, " C2'" );
  strcpy( AA[p].atom[10].atom_name, " C1'" );

  strcpy( AA[p].atom[11].atom_name, " N1 " );
  strcpy( AA[p].atom[12].atom_name, " C2 " );
  strcpy( AA[p].atom[13].atom_name, " O2 " );
  strcpy( AA[p].atom[14].atom_name, " N3 " );
  strcpy( AA[p].atom[15].atom_name, " C4 " );
  strcpy( AA[p].atom[16].atom_name, " N4 " );
  strcpy( AA[p].atom[17].atom_name, " C5 " );
  strcpy( AA[p].atom[18].atom_name, " C6 " );

//  strcpy( AA[p].atom[19].atom_name, "1H5*" );
//  strcpy( AA[p].atom[20].atom_name, "2H5*" );
//  strcpy( AA[p].atom[21].atom_name, " H4*" );
//  strcpy( AA[p].atom[22].atom_name, " H3*" );
//  strcpy( AA[p].atom[23].atom_name, "1H2*" );
//  strcpy( AA[p].atom[24].atom_name, "2H2*" );
//  strcpy( AA[p].atom[25].atom_name, " H1*" );

  strcpy( AA[p].atom[19].atom_name, "1H5'" );
  strcpy( AA[p].atom[20].atom_name, "2H5'" );
  strcpy( AA[p].atom[21].atom_name, " H4'" );
  strcpy( AA[p].atom[22].atom_name, " H3'" );
  strcpy( AA[p].atom[23].atom_name, "1H2'" );
  strcpy( AA[p].atom[24].atom_name, "2H2'" );
  strcpy( AA[p].atom[25].atom_name, " H1'" );

  strcpy( AA[p].atom[26].atom_name, "1H4 " );
  strcpy( AA[p].atom[27].atom_name, "2H4 " );
  strcpy( AA[p].atom[28].atom_name, " H5 " );
  strcpy( AA[p].atom[29].atom_name, " H6 " );

  AA[p].atom[0].charge = 1.500; // P
  AA[p].atom[1].charge = -.780; // O1P
  AA[p].atom[2].charge = -.780; // O2P
  AA[p].atom[3].charge = -.570; // O5*
  AA[p].atom[4].charge = -.080; // C5*
  AA[p].atom[5].charge =.160; // C4*
  AA[p].atom[6].charge = -.500; // O4*
  AA[p].atom[7].charge =.010; // C3*
  AA[p].atom[8].charge = -.570; // O3*
  AA[p].atom[9].charge = -.180; // C2*
  AA[p].atom[10].charge =.160; // C1*
  AA[p].atom[11].charge = -.130; // N1
  AA[p].atom[12].charge =.520; // C2
  AA[p].atom[13].charge = -.490; // O2
  AA[p].atom[14].charge = -.660; // N3
  AA[p].atom[15].charge =.650; // C4
  AA[p].atom[16].charge = -.750; // N4
  AA[p].atom[17].charge = -.130; // C5
  AA[p].atom[18].charge =.050; // C6
  AA[p].atom[19].charge =.090; // 1H5*
  AA[p].atom[20].charge =.090; // 2H5*
  AA[p].atom[21].charge =.090; // H4*
  AA[p].atom[22].charge =.090; // H3*
  AA[p].atom[23].charge =.090; // 1H2*
  AA[p].atom[24].charge =.090; // 2H2*
  AA[p].atom[25].charge =.090; // H1*
  AA[p].atom[26].charge =.370; // 1H4
  AA[p].atom[27].charge =.330; // 2H4
  AA[p].atom[28].charge =.070; // H5
  AA[p].atom[29].charge =.170; // H6


  switch(opt)
  {
  case Sybil:
	  AA[p].atom[0].fullatom_type = 12; // Phos   P
	  AA[p].atom[1].fullatom_type = 10; // OOC    O1P
	  AA[p].atom[2].fullatom_type = 10; // OOC    O2P
	  AA[p].atom[3].fullatom_type = 9; // ONH2   O5*
	  AA[p].atom[4].fullatom_type = 2; // CH2    C5*
	  AA[p].atom[5].fullatom_type = 2; // CH1    C4*
	  AA[p].atom[6].fullatom_type = 9; // OH     O4*
	  AA[p].atom[7].fullatom_type = 2; // CH1    C3*
	  AA[p].atom[8].fullatom_type = 9; // ONH2   O3*
	  AA[p].atom[9].fullatom_type = 2; // CH2    C2*
	  AA[p].atom[10].fullatom_type = 2; // CH1    C1*
	  AA[p].atom[11].fullatom_type = 5; // Nhis   N1
	  AA[p].atom[12].fullatom_type = 1; // CObb   C2
	  AA[p].atom[13].fullatom_type = 8; // OCbb   O2
	  AA[p].atom[14].fullatom_type = 5; // Nhis   N3
	  AA[p].atom[15].fullatom_type = 3; // aroC   C4
	  AA[p].atom[16].fullatom_type = 5; // NH2O   N4
	  AA[p].atom[17].fullatom_type = 3; // aroC   C5
	  AA[p].atom[18].fullatom_type = 3; // aroC   C6
	  AA[p].atom[19].fullatom_type = 14; // Hapo  1H5*
	  AA[p].atom[20].fullatom_type = 14; // Hapo  2H5*
	  AA[p].atom[21].fullatom_type = 14; // Hapo   H4*
	  AA[p].atom[22].fullatom_type = 14; // Hapo   H3*
	  AA[p].atom[23].fullatom_type = 14; // Hapo  1H2*
	  AA[p].atom[24].fullatom_type = 14; // Hapo  2H2*
	  AA[p].atom[25].fullatom_type = 14; // Hapo   H1*
	  AA[p].atom[26].fullatom_type = 13; // Hpol  1H4
	  AA[p].atom[27].fullatom_type = 13; // Hpol  2H4
	  AA[p].atom[28].fullatom_type = 15; // Haro   H5
	  AA[p].atom[29].fullatom_type = 15; // Haro   H6
	  break;
  default:
	  AA[p].atom[0].fullatom_type = 21; // Phos   P
	  AA[p].atom[1].fullatom_type = 15; // OOC    O1P
	  AA[p].atom[2].fullatom_type = 15; // OOC    O2P
	  AA[p].atom[3].fullatom_type = 14; // ONH2   O5*
	  AA[p].atom[4].fullatom_type = 4; // CH2    C5*
	  AA[p].atom[5].fullatom_type = 3; // CH1    C4*
	  AA[p].atom[6].fullatom_type = 13; // OH     O4*
	  AA[p].atom[7].fullatom_type = 3; // CH1    C3*
	  AA[p].atom[8].fullatom_type = 14; // ONH2   O3*
	  AA[p].atom[9].fullatom_type = 4; // CH2    C2*
	  AA[p].atom[10].fullatom_type = 3; // CH1    C1*
	  AA[p].atom[11].fullatom_type = 8; // Nhis   N1
	  AA[p].atom[12].fullatom_type = 19; // CObb   C2
	  AA[p].atom[13].fullatom_type = 20; // OCbb   O2
	  AA[p].atom[14].fullatom_type = 8; // Nhis   N3
	  AA[p].atom[15].fullatom_type = 6; // aroC   C4
	  AA[p].atom[16].fullatom_type = 9; // NH2O   N4
	  AA[p].atom[17].fullatom_type = 6; // aroC   C5
	  AA[p].atom[18].fullatom_type = 6; // aroC   C6
	  AA[p].atom[19].fullatom_type = 23; // Hapo  1H5*
	  AA[p].atom[20].fullatom_type = 23; // Hapo  2H5*
	  AA[p].atom[21].fullatom_type = 23; // Hapo   H4*
	  AA[p].atom[22].fullatom_type = 23; // Hapo   H3*
	  AA[p].atom[23].fullatom_type = 23; // Hapo  1H2*
	  AA[p].atom[24].fullatom_type = 23; // Hapo  2H2*
	  AA[p].atom[25].fullatom_type = 23; // Hapo   H1*
	  AA[p].atom[26].fullatom_type = 22; // Hpol  1H4
	  AA[p].atom[27].fullatom_type = 22; // Hpol  2H4
	  AA[p].atom[28].fullatom_type = 24; // Haro   H5
	  AA[p].atom[29].fullatom_type = 24; // Haro   H6
	  break;
  }
  //jjh intra residue bonding

  AA[p].atom[0].nbonded_neighbors = 3; // P
  AA[p].atom[0].bonded_neighbor[0] = 1; // P--O1P
  AA[p].atom[0].bonded_neighbor[1] = 2; // P--O2P
  AA[p].atom[0].bonded_neighbor[2] = 3; // P--O5*
  AA[p].atom[1].nbonded_neighbors = 1; // O1P
  AA[p].atom[1].bonded_neighbor[0] = 0; // O1P--P
  AA[p].atom[2].nbonded_neighbors = 1; // O2P
  AA[p].atom[2].bonded_neighbor[0] = 0; // O2P--P
  AA[p].atom[3].nbonded_neighbors = 2; // O5*
  AA[p].atom[3].bonded_neighbor[0] = 0; // O5*--P
  AA[p].atom[3].bonded_neighbor[1] = 4; // O5*--C5*
  AA[p].atom[4].nbonded_neighbors = 4; // C5*
  AA[p].atom[4].bonded_neighbor[0] = 3; // C5*--O5*
  AA[p].atom[4].bonded_neighbor[1] = 5; // C5*--C4*
  AA[p].atom[4].bonded_neighbor[2] = 19; // C5*--1H5*
  AA[p].atom[4].bonded_neighbor[3] = 20; // C5*--2H5*
  AA[p].atom[5].nbonded_neighbors = 4; // C4*
  AA[p].atom[5].bonded_neighbor[0] = 4; // C4*--C5*
  AA[p].atom[5].bonded_neighbor[1] = 6; // C4*--O4*
  AA[p].atom[5].bonded_neighbor[2] = 7; // C4*--C3*
  AA[p].atom[5].bonded_neighbor[3] = 21; // C4*--H4*
  AA[p].atom[6].nbonded_neighbors = 2; // O4*
  AA[p].atom[6].bonded_neighbor[0] = 5; // O4*--C4*
  AA[p].atom[6].bonded_neighbor[1] = 10; // O4*--C1*
  AA[p].atom[7].nbonded_neighbors = 4; // C3*
  AA[p].atom[7].bonded_neighbor[0] = 5; // C3*--C4*
  AA[p].atom[7].bonded_neighbor[1] = 8; // C3*--O3*
  AA[p].atom[7].bonded_neighbor[2] = 9; // C3*--C2*
  AA[p].atom[7].bonded_neighbor[3] = 22; // C3*--H3*
  AA[p].atom[8].nbonded_neighbors = 1; // O3*
  AA[p].atom[8].bonded_neighbor[0] = 7; // O3*--C3*
  AA[p].atom[9].nbonded_neighbors = 4; // C2*
  AA[p].atom[9].bonded_neighbor[0] = 7; // C2*--C3*
  AA[p].atom[9].bonded_neighbor[1] = 10; // C2*--C1*
  AA[p].atom[9].bonded_neighbor[2] = 23; // C2*--1H2*
  AA[p].atom[9].bonded_neighbor[3] = 24; // C2*--2H2*
  AA[p].atom[10].nbonded_neighbors = 4; // C1*
  AA[p].atom[10].bonded_neighbor[0] = 6; // C1*--O4*
  AA[p].atom[10].bonded_neighbor[1] = 9; // C1*--C2*
  AA[p].atom[10].bonded_neighbor[2] = 11; // C1*--N1
  AA[p].atom[10].bonded_neighbor[3] = 25; // C1*--H1*
  AA[p].atom[11].nbonded_neighbors = 3; // N1
  AA[p].atom[11].bonded_neighbor[0] = 10; // N1--C1*
  AA[p].atom[11].bonded_neighbor[1] = 12; // N1--C2
  AA[p].atom[11].bonded_neighbor[2] = 18; // N1--C6
  AA[p].atom[12].nbonded_neighbors = 3; // C2
  AA[p].atom[12].bonded_neighbor[0] = 11; // C2--N1
  AA[p].atom[12].bonded_neighbor[1] = 13; // C2--O2
  AA[p].atom[12].bonded_neighbor[2] = 14; // C2--N3
  AA[p].atom[13].nbonded_neighbors = 1; // O2
  AA[p].atom[13].bonded_neighbor[0] = 12; // O2--C2
  AA[p].atom[14].nbonded_neighbors = 2; // N3
  AA[p].atom[14].bonded_neighbor[0] = 12; // N3--C2
  AA[p].atom[14].bonded_neighbor[1] = 15; // N3--C4
  AA[p].atom[15].nbonded_neighbors = 3; // C4
  AA[p].atom[15].bonded_neighbor[0] = 14; // C4--N3
  AA[p].atom[15].bonded_neighbor[1] = 16; // C4--N4
  AA[p].atom[15].bonded_neighbor[2] = 17; // C4--C5
  AA[p].atom[16].nbonded_neighbors = 3; // N4
  AA[p].atom[16].bonded_neighbor[0] = 15; // N4--C4
  AA[p].atom[16].bonded_neighbor[1] = 26; // N4--1H4
  AA[p].atom[16].bonded_neighbor[2] = 27; // N4--2H4
  AA[p].atom[17].nbonded_neighbors = 3; // C5
  AA[p].atom[17].bonded_neighbor[0] = 15; // C5--C4
  AA[p].atom[17].bonded_neighbor[1] = 18; // C5--C6
  AA[p].atom[17].bonded_neighbor[2] = 28; // C5--H5
  AA[p].atom[18].nbonded_neighbors = 3; // C6
  AA[p].atom[18].bonded_neighbor[0] = 11; // C6--N1
  AA[p].atom[18].bonded_neighbor[1] = 17; // C6--C5
  AA[p].atom[18].bonded_neighbor[2] = 29; // C6--H6
  AA[p].atom[19].nbonded_neighbors = 1; //1H5*
  AA[p].atom[19].bonded_neighbor[0] = 4; //1H5*--C5*
  AA[p].atom[20].nbonded_neighbors = 1; //2H5*
  AA[p].atom[20].bonded_neighbor[0] = 4; //2H5*--C5*
  AA[p].atom[21].nbonded_neighbors = 1; // H4*
  AA[p].atom[21].bonded_neighbor[0] = 5; // H4*--C4*
  AA[p].atom[22].nbonded_neighbors = 1; // H3*
  AA[p].atom[22].bonded_neighbor[0] = 7; // H3*--C3*
  AA[p].atom[23].nbonded_neighbors = 1; //1H2*
  AA[p].atom[23].bonded_neighbor[0] = 9; // 1H2*--C2*
  AA[p].atom[24].nbonded_neighbors = 1; //2H2*
  AA[p].atom[24].bonded_neighbor[0] = 9; // 2H2*--C2*
  AA[p].atom[25].nbonded_neighbors = 1; // H1*
  AA[p].atom[25].bonded_neighbor[0] = 10; // H1*--C1*
  AA[p].atom[26].nbonded_neighbors = 1; //1H4
  AA[p].atom[26].bonded_neighbor[0] = 16; //1H4--N4
  AA[p].atom[27].nbonded_neighbors = 1; //2H4
  AA[p].atom[27].bonded_neighbor[0] = 16; //2H4--N4
  AA[p].atom[28].nbonded_neighbors = 1; // H5
  AA[p].atom[28].bonded_neighbor[0] = 17; // H5--C5
  AA[p].atom[29].nbonded_neighbors = 1; // H6
  AA[p].atom[29].bonded_neighbor[0] = 18; // H6--C6


  for ( k = 0; k < 3; k++ )
    AA[p].atom[k].hastemplate = true;

  for ( k = 3; k < 21; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 21; k < 30; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template atoms used for placing atoms CYT

  //bk   template for building P
  AA[p].atom[0].ta[0] = 3; //   O5*
  AA[p].atom[0].ta[1] = 4; //   C5*
  AA[p].atom[0].ta[2] = 5; //   C4*

  //bk   template for building O1P
  AA[p].atom[1].ta[0] = 0; //   P
  AA[p].atom[1].ta[1] = 4; //   C5*
  AA[p].atom[1].ta[2] = 3; //   O5*

  //bk   template for building  O2P
  AA[p].atom[2].ta[0] = 0; //   P
  AA[p].atom[2].ta[1] = 1; //   O1P
  AA[p].atom[2].ta[2] = 3; //   O5*

  //bk   template for building 1H5*
  AA[p].atom[21].ta[0] = 4; //   C5*
  AA[p].atom[21].ta[1] = 3; //   O5*
  AA[p].atom[21].ta[2] = 5; //   C4*

  //bk   template for building 2H5*
  AA[p].atom[22].ta[0] = 4; //   C5*
  AA[p].atom[22].ta[1] = 3; //   O5*
  AA[p].atom[22].ta[2] = 5; //   C4*

  //bk   template for building H4*
  AA[p].atom[23].ta[0] = 5; //   C4*
  AA[p].atom[23].ta[1] = 6; //   O4*
  AA[p].atom[23].ta[2] = 4; //   C5*

  //bk   template for building H3*
  AA[p].atom[24].ta[0] = 7; //   C3*
  AA[p].atom[24].ta[1] = 9; //   C2*
  AA[p].atom[24].ta[2] = 5; //   C4*

  //bk   template for building 1H2*
  AA[p].atom[25].ta[0] = 9; //   C2*
  AA[p].atom[25].ta[1] = 7; //   C3*
  AA[p].atom[25].ta[2] = 10; //   C1*

  //bk   template for building 2H2*
  AA[p].atom[26].ta[0] = 9; //   C2*
  AA[p].atom[26].ta[1] = 7; //   C3*
  AA[p].atom[26].ta[2] = 10; //   C1*

  //bk   template for building H1*
  AA[p].atom[27].ta[0] = 10; //   C1*
  AA[p].atom[27].ta[1] = 9; //   C2*
  AA[p].atom[27].ta[2] = 6; //   O4*

  //bk   template for building 1H4
  AA[p].atom[28].ta[0] = 16; //   N4
  AA[p].atom[28].ta[1] = 15; //   C4
  AA[p].atom[28].ta[2] = 14; //   N3

  //bk   template for building 2H4
  AA[p].atom[29].ta[0] = 16; //   N4
  AA[p].atom[29].ta[1] = 15; //   C4
  AA[p].atom[29].ta[2] = 14; //   N3

// Mon fixed (11/12/2009)
//  //bk   template for building H5
//  AA[p].atom[30].ta[0] = 17; //   C5
//  AA[p].atom[30].ta[1] = 18; //   C6
//  AA[p].atom[30].ta[2] = 11; //   N1
//
//  //bk   template for building H6
//  AA[p].atom[31].ta[0] = 18; //   C6
//  AA[p].atom[31].ta[1] = 17; //   C5
//  AA[p].atom[31].ta[2] = 15; //   C4


  //bk   template coordinates for the amino acid
  //jjh new template from Kevin Karplus

  AA[p].atom[0].icoor[0] = 21.1970; //   P
  AA[p].atom[0].icoor[1] = 31.2910; //   P
  AA[p].atom[0].icoor[2] = 5.2860; //   P
  AA[p].atom[1].icoor[0] = 22.0680; //   O1P
  AA[p].atom[1].icoor[1] = 32.4840; //   O1P
  AA[p].atom[1].icoor[2] = 5.3940; //   O1P
  AA[p].atom[2].icoor[0] = 21.5980; //   O2P
  AA[p].atom[2].icoor[1] = 30.0330; //   O2P
  AA[p].atom[2].icoor[2] = 5.9560; //   O2P
  AA[p].atom[3].icoor[0] = 19.7430; //   O5*
  AA[p].atom[3].icoor[1] = 31.6880; //   O5*
  AA[p].atom[3].icoor[2] = 5.7980; //   O5*
  AA[p].atom[4].icoor[0] = 18.7410; //   C5*
  AA[p].atom[4].icoor[1] = 32.0840; //   C5*
  AA[p].atom[4].icoor[2] = 5.0180; //   C5*
  AA[p].atom[5].icoor[0] = 17.5950; //   C4*
  AA[p].atom[5].icoor[1] = 32.7830; //   C4*
  AA[p].atom[5].icoor[2] = 5.7080; //   C4*
  AA[p].atom[6].icoor[0] = 16.7700; //   O4*
  AA[p].atom[6].icoor[1] = 31.6020; //   O4*
  AA[p].atom[6].icoor[2] = 5.5630; //   O4*
  AA[p].atom[7].icoor[0] = 17.6890; //   C3*
  AA[p].atom[7].icoor[1] = 33.0610; //   C3*
  AA[p].atom[7].icoor[2] = 7.2080; //   C3*
  AA[p].atom[8].icoor[0] = 16.7660; //   O3*
  AA[p].atom[8].icoor[1] = 34.0930; //   O3*
  AA[p].atom[8].icoor[2] = 7.5790; //   O3*
  AA[p].atom[9].icoor[0] = 17.3280; //   C2*
  AA[p].atom[9].icoor[1] = 31.7270; //   C2*
  AA[p].atom[9].icoor[2] = 7.8400; //   C2*
  AA[p].atom[10].icoor[0] = 16.3760; //   C1*
  AA[p].atom[10].icoor[1] = 31.1040; //   C1*
  AA[p].atom[10].icoor[2] = 6.8310; //   C1*
  AA[p].atom[11].icoor[0] = 16.4250; //   N1
  AA[p].atom[11].icoor[1] = 29.6330; //   N1
  AA[p].atom[11].icoor[2] = 6.7640; //   N1
  AA[p].atom[12].icoor[0] = 15.2760; //   C2
  AA[p].atom[12].icoor[1] = 28.8900; //   C2
  AA[p].atom[12].icoor[2] = 6.8540; //   C2
  AA[p].atom[13].icoor[0] = 14.1820; //   O2
  AA[p].atom[13].icoor[1] = 29.4550; //   O2
  AA[p].atom[13].icoor[2] = 6.9990; //   O2
  AA[p].atom[14].icoor[0] = 15.3600; //   N3
  AA[p].atom[14].icoor[1] = 27.5760; //   N3
  AA[p].atom[14].icoor[2] = 6.7900; //   N3
  AA[p].atom[15].icoor[0] = 16.5320; //   C4
  AA[p].atom[15].icoor[1] = 27.0030; //   C4
  AA[p].atom[15].icoor[2] = 6.5920; //   C4
  AA[p].atom[16].icoor[0] = 16.5640; //   N4
  AA[p].atom[16].icoor[1] = 25.6760; //   N4
  AA[p].atom[16].icoor[2] = 6.5150; //   N4
  AA[p].atom[17].icoor[0] = 17.7100; //   C5
  AA[p].atom[17].icoor[1] = 27.7410; //   C5
  AA[p].atom[17].icoor[2] = 6.4480; //   C5
  AA[p].atom[18].icoor[0] = 17.6100; //   C6
  AA[p].atom[18].icoor[1] = 29.0570; //   C6
  AA[p].atom[18].icoor[2] = 6.5230; //   C6
  AA[p].atom[19].icoor[0] = 18.7410; //  1H5*
  AA[p].atom[19].icoor[1] = 32.0840; //  1H5*
  AA[p].atom[19].icoor[2] = 4.0250; //  1H5*
  AA[p].atom[20].icoor[0] = 19.4390; //  2H5*
  AA[p].atom[20].icoor[1] = 33.5210; //  2H5*
  AA[p].atom[20].icoor[2] = 4.8770; //  2H5*
  AA[p].atom[21].icoor[0] = 17.0610; //   H4*
  AA[p].atom[21].icoor[1] = 33.6400; //   H4*
  AA[p].atom[21].icoor[2] = 5.2440; //   H4*
  AA[p].atom[22].icoor[0] = 18.7200; //   H3*
  AA[p].atom[22].icoor[1] = 33.3710; //   H3*
  AA[p].atom[22].icoor[2] = 7.4730; //   H3*
  AA[p].atom[23].icoor[0] = 18.2300; //  1H2*
  AA[p].atom[23].icoor[1] = 31.0940; //  1H2*
  AA[p].atom[23].icoor[2] = 7.9820; //  1H2*
  AA[p].atom[24].icoor[0] = 16.8320; //  2H2*
  AA[p].atom[24].icoor[1] = 31.8710; //  2H2*
  AA[p].atom[24].icoor[2] = 8.8230; //  2H2*
  AA[p].atom[25].icoor[0] = 15.3260; //   H1*
  AA[p].atom[25].icoor[1] = 31.4080; //   H1*
  AA[p].atom[25].icoor[2] = 7.0300; //   H1*
  AA[p].atom[26].icoor[0] = 15.6870; //  1H4
  AA[p].atom[26].icoor[1] = 25.1890; //  1H4
  AA[p].atom[26].icoor[2] = 6.5910; //  1H4
  AA[p].atom[27].icoor[0] = 17.3570; //  2H4
  AA[p].atom[27].icoor[1] = 25.1800; //  2H4
  AA[p].atom[27].icoor[2] = 6.1750; //  2H4
  AA[p].atom[28].icoor[0] = 18.6520; //   H5
  AA[p].atom[28].icoor[1] = 27.2310; //   H5
  AA[p].atom[28].icoor[2] = 6.2550; //   H5
  AA[p].atom[29].icoor[0] = 18.4730; //   H6
  AA[p].atom[29].icoor[1] = 29.7230; //   H6
  AA[p].atom[29].icoor[2] = 6.3990; //   H6

  //FIN DE CYT

  //THY
//  p = THY;
  p = DTHY;

//  strcpy( AA[p].aa_name3, "  T" );
  strcpy( AA[p].aa_name3, " DT" );
//  AA[p].aa_name1 = 'X';
  AA[p].aa_name1 = 'T';

  AA[p].mass = 0;

  AA[p].is_protein = false;
  AA[p].is_RNA = false;
  AA[p].is_DNA = true;

  AA[p].aa_is_polar = true;
  AA[p].aa_is_nonpolar = false;
  AA[p].aa_is_aromatic = true;
  AA[p].aa_is_charged = true;

  AA[p].natoms = 32;
  AA[p].nheavyatoms = 20;
  AA[p].nchi = 0;

  AA[p].chi_atoms = NULL;
  AA[p].atom = ( t_aa_atom_p * ) malloc( AA[p].natoms * sizeof( t_aa_atom_p ) );


  strcpy( AA[p].atom[0].atom_name, " P  " );
//  strcpy( AA[p].atom[1].atom_name, " O1P" );
//  strcpy( AA[p].atom[2].atom_name, " O2P" );

  strcpy( AA[p].atom[1].atom_name, " OP1" );
  strcpy( AA[p].atom[2].atom_name, " OP2" );

//  strcpy( AA[p].atom[3].atom_name, " O5*" );
//  strcpy( AA[p].atom[4].atom_name, " C5*" );
//  strcpy( AA[p].atom[5].atom_name, " C4*" );
//  strcpy( AA[p].atom[6].atom_name, " O4*" );
//  strcpy( AA[p].atom[7].atom_name, " C3*" );
//  strcpy( AA[p].atom[8].atom_name, " O3*" );
//  strcpy( AA[p].atom[9].atom_name, " C2*" );
//  strcpy( AA[p].atom[10].atom_name, " C1*" );

  strcpy( AA[p].atom[3].atom_name, " O5'" );
  strcpy( AA[p].atom[4].atom_name, " C5'" );
  strcpy( AA[p].atom[5].atom_name, " C4'" );
  strcpy( AA[p].atom[6].atom_name, " O4'" );
  strcpy( AA[p].atom[7].atom_name, " C3'" );
  strcpy( AA[p].atom[8].atom_name, " O3'" );
  strcpy( AA[p].atom[9].atom_name, " C2'" );
  strcpy( AA[p].atom[10].atom_name, " C1'" );

  strcpy( AA[p].atom[11].atom_name, " N1 " );
  strcpy( AA[p].atom[12].atom_name, " C2 " );
  strcpy( AA[p].atom[13].atom_name, " O2 " );
  strcpy( AA[p].atom[14].atom_name, " N3 " );
  strcpy( AA[p].atom[15].atom_name, " C4 " );
  strcpy( AA[p].atom[16].atom_name, " O4 " );
  strcpy( AA[p].atom[17].atom_name, " C5 " );
//strcpy( AA[p].atom[18].atom_name, " C5M" ); // alias atom DT C7 C5M
  strcpy( AA[p].atom[18].atom_name, " C7 " ); // alias atom DT C7 C5M
  strcpy( AA[p].atom[19].atom_name, " C6 " );

//  strcpy( AA[p].atom[20].atom_name, "1H5*" );
//  strcpy( AA[p].atom[21].atom_name, "2H5*" );
//  strcpy( AA[p].atom[22].atom_name, " H4*" );
//  strcpy( AA[p].atom[23].atom_name, " H3*" );
//  strcpy( AA[p].atom[24].atom_name, "1H2*" );
//  strcpy( AA[p].atom[25].atom_name, "2H2*" );
//  strcpy( AA[p].atom[26].atom_name, " H1*" );

  strcpy( AA[p].atom[20].atom_name, "1H5'" );
  strcpy( AA[p].atom[21].atom_name, "2H5'" );
  strcpy( AA[p].atom[22].atom_name, " H4'" );
  strcpy( AA[p].atom[23].atom_name, " H3'" );
  strcpy( AA[p].atom[24].atom_name, "1H2'" );
  strcpy( AA[p].atom[25].atom_name, "2H2'" );
  strcpy( AA[p].atom[26].atom_name, " H1'" );

  strcpy( AA[p].atom[27].atom_name, " H3 " );
//  strcpy( AA[p].atom[28].atom_name, "1H5M" ); // alias atom T 1H5M H51
//  strcpy( AA[p].atom[29].atom_name, "2H5M" );
//  strcpy( AA[p].atom[30].atom_name, "3H5M" );
  strcpy( AA[p].atom[28].atom_name, " H51" ); // alias atom T 1H5M H51
  strcpy( AA[p].atom[29].atom_name, " H52" );
  strcpy( AA[p].atom[30].atom_name, " H53" );
  strcpy( AA[p].atom[31].atom_name, " H6 " );

  AA[p].atom[0].charge = 1.500; // P
  AA[p].atom[1].charge = -.780; // O1P
  AA[p].atom[2].charge = -.780; // O2P
  AA[p].atom[3].charge = -.570; // O5*
  AA[p].atom[4].charge = -.080; // C5*
  AA[p].atom[5].charge =.160; // C4*
  AA[p].atom[6].charge = -.500; // O4*
  AA[p].atom[7].charge =.010; // C3*
  AA[p].atom[8].charge = -.570; // O3*
  AA[p].atom[9].charge = -.180; // C2*
  AA[p].atom[10].charge =.160; // C1*
  AA[p].atom[11].charge = -.340; // N1
  AA[p].atom[12].charge =.510; // C2
  AA[p].atom[13].charge = -.410; // O2
  AA[p].atom[14].charge = -.460; // N3
  AA[p].atom[15].charge =.500; // C4
  AA[p].atom[16].charge = -.450; // O4
  AA[p].atom[17].charge = -.150; // C5
  AA[p].atom[18].charge = -.110; // C5M
  AA[p].atom[19].charge =.170; // C6
  AA[p].atom[20].charge =.090; // 1H5*
  AA[p].atom[21].charge =.090; // 2H5*
  AA[p].atom[22].charge =.090; // H4*
  AA[p].atom[23].charge =.090; // H3*
  AA[p].atom[24].charge =.090; // 1H2*
  AA[p].atom[25].charge =.090; // 2H2*
  AA[p].atom[26].charge =.090; // H1*
  AA[p].atom[27].charge =.360; // H3
  AA[p].atom[28].charge =.070; // 1H5M
  AA[p].atom[29].charge =.070; // 2H5M
  AA[p].atom[30].charge =.070; // 3H5M
  AA[p].atom[31].charge =.170; // H6


  switch(opt)
  {
  default:
	  AA[p].atom[0].fullatom_type = 12; // Phos   P
	  AA[p].atom[1].fullatom_type = 10; // OOC    O1P
	  AA[p].atom[2].fullatom_type = 10; // OOC    O2P
	  AA[p].atom[3].fullatom_type = 9; // ONH2   O5*
	  AA[p].atom[4].fullatom_type = 2; // CH2    C5*
	  AA[p].atom[5].fullatom_type = 2; // CH1    C4*
	  AA[p].atom[6].fullatom_type = 9; // OH     O4*
	  AA[p].atom[7].fullatom_type = 2; // CH1    C3*
	  AA[p].atom[8].fullatom_type = 9; // ONH2   O3*
	  AA[p].atom[9].fullatom_type = 2; // CH2    C2*
	  AA[p].atom[10].fullatom_type = 2; // CH1    C1*
	  AA[p].atom[11].fullatom_type = 5; // Nhis   N1
	  AA[p].atom[12].fullatom_type = 1; // CObb   C2
	  AA[p].atom[13].fullatom_type = 8; // OCbb   O2
	  AA[p].atom[14].fullatom_type = 5; // Ntrp   N3
	  AA[p].atom[15].fullatom_type = 1; // CObb   C4
	  AA[p].atom[16].fullatom_type = 8; // OCbb   O4
	  AA[p].atom[17].fullatom_type = 3; // aroC   C5
	  AA[p].atom[18].fullatom_type = 2; // CH3    C5M
	  AA[p].atom[19].fullatom_type = 3; // aroC   C6
	  AA[p].atom[20].fullatom_type = 14; // Hapo  1H5*
	  AA[p].atom[21].fullatom_type = 14; // Hapo  2H5*
	  AA[p].atom[22].fullatom_type = 14; // Hapo   H4*
	  AA[p].atom[23].fullatom_type = 14; // Hapo   H3*
	  AA[p].atom[24].fullatom_type = 14; // Hapo  1H2*
	  AA[p].atom[25].fullatom_type = 14; // Hapo  2H2*
	  AA[p].atom[26].fullatom_type = 14; // Hapo   H1*
	  AA[p].atom[27].fullatom_type = 13; // Hpol   H3
	  AA[p].atom[28].fullatom_type = 14; // Hapo  1H5M
	  AA[p].atom[29].fullatom_type = 14; // Hapo  2H5M
	  AA[p].atom[30].fullatom_type = 14; // Hapo  3H5M
	  AA[p].atom[31].fullatom_type = 15; // Haro   H6
	  break;
  }
  //jjh intra residue bonding
  AA[p].atom[0].nbonded_neighbors = 3; // P
  AA[p].atom[0].bonded_neighbor[0] = 1; // P--O1P
  AA[p].atom[0].bonded_neighbor[1] = 2; // P--O2P
  AA[p].atom[0].bonded_neighbor[2] = 3; // P--O5*
  AA[p].atom[1].nbonded_neighbors = 1; // O1P
  AA[p].atom[1].bonded_neighbor[0] = 0; // O1P--P
  AA[p].atom[2].nbonded_neighbors = 1; // O2P
  AA[p].atom[2].bonded_neighbor[0] = 0; // O2P--P
  AA[p].atom[3].nbonded_neighbors = 2; // O5*
  AA[p].atom[3].bonded_neighbor[0] = 0; // O5*--P
  AA[p].atom[3].bonded_neighbor[1] = 4; // O5*--C5*
  AA[p].atom[4].nbonded_neighbors = 4; // C5*
  AA[p].atom[4].bonded_neighbor[0] = 3; // C5*--O5*
  AA[p].atom[4].bonded_neighbor[1] = 5; // C5*--C4*
  AA[p].atom[4].bonded_neighbor[2] = 20; // C5*--1H5
  AA[p].atom[4].bonded_neighbor[3] = 21; // C5*--2H5
  AA[p].atom[5].nbonded_neighbors = 4; // C4*
  AA[p].atom[5].bonded_neighbor[0] = 4; // C4*--C5*
  AA[p].atom[5].bonded_neighbor[1] = 6; // C4*--O4*
  AA[p].atom[5].bonded_neighbor[2] = 7; // C4*--C3*
  AA[p].atom[5].bonded_neighbor[3] = 22; // C4*--H4*
  AA[p].atom[6].nbonded_neighbors = 2; // O4*
  AA[p].atom[6].bonded_neighbor[0] = 5; // O4*--C4*
  AA[p].atom[6].bonded_neighbor[1] = 10; // O4*--C1*
  AA[p].atom[7].nbonded_neighbors = 4; // C3*
  AA[p].atom[7].bonded_neighbor[0] = 5; // C3*--C4*
  AA[p].atom[7].bonded_neighbor[1] = 8; // C3*--O3*
  AA[p].atom[7].bonded_neighbor[2] = 9; // C3*--C2*
  AA[p].atom[7].bonded_neighbor[3] = 23; // C3*--H3*
  AA[p].atom[8].nbonded_neighbors = 1; // O3*
  AA[p].atom[8].bonded_neighbor[0] = 7; // O3*--C3*
  AA[p].atom[9].nbonded_neighbors = 4; // C2*
  AA[p].atom[9].bonded_neighbor[0] = 7; // C2*--C3*
  AA[p].atom[9].bonded_neighbor[1] = 10; // C2*--C1*
  AA[p].atom[9].bonded_neighbor[2] = 24; // C2*--1H2*
  AA[p].atom[9].bonded_neighbor[3] = 25; // C2*--2H2*
  AA[p].atom[10].nbonded_neighbors = 4; // C1*
  AA[p].atom[10].bonded_neighbor[0] = 6; // C1*--O4*
  AA[p].atom[10].bonded_neighbor[1] = 9; // C1*--C2*
  AA[p].atom[10].bonded_neighbor[2] = 11; // C1*--N1
  AA[p].atom[10].bonded_neighbor[3] = 26; // C1*--H1*
  AA[p].atom[11].nbonded_neighbors = 3; // N1
  AA[p].atom[11].bonded_neighbor[0] = 10; // N1--C1*
  AA[p].atom[11].bonded_neighbor[1] = 12; // N1--C2
  AA[p].atom[11].bonded_neighbor[2] = 19; // N1--C6
  AA[p].atom[12].nbonded_neighbors = 3; // C2
  AA[p].atom[12].bonded_neighbor[0] = 11; // C2--N1
  AA[p].atom[12].bonded_neighbor[1] = 13; // C2--O2
  AA[p].atom[12].bonded_neighbor[2] = 14; // C2--N3
  AA[p].atom[13].nbonded_neighbors = 1; // O2
  AA[p].atom[13].bonded_neighbor[0] = 12; // O2--C2
  AA[p].atom[14].nbonded_neighbors = 3; // N3
  AA[p].atom[14].bonded_neighbor[0] = 12; // N3--C2
  AA[p].atom[14].bonded_neighbor[1] = 15; // N3--C4
  AA[p].atom[14].bonded_neighbor[2] = 27; // N3--H3
  AA[p].atom[15].nbonded_neighbors = 3; // C4
  AA[p].atom[15].bonded_neighbor[0] = 14; // C4--N3
  AA[p].atom[15].bonded_neighbor[1] = 16; // C4--O4
  AA[p].atom[15].bonded_neighbor[2] = 17; // C4--C5
  AA[p].atom[16].nbonded_neighbors = 1; // O4
  AA[p].atom[16].bonded_neighbor[0] = 15; // O4--C4
  AA[p].atom[17].nbonded_neighbors = 3; // C5
  AA[p].atom[17].bonded_neighbor[0] = 15; // C5--C4
  AA[p].atom[17].bonded_neighbor[1] = 18; // C5--C5M
  AA[p].atom[17].bonded_neighbor[2] = 19; // C5--C6
  AA[p].atom[18].nbonded_neighbors = 4; // C5M
  AA[p].atom[18].bonded_neighbor[0] = 17; // C5M--C5
  AA[p].atom[18].bonded_neighbor[1] = 28; // C5M--1H5M
  AA[p].atom[18].bonded_neighbor[2] = 29; // C5M--2H5M
  AA[p].atom[18].bonded_neighbor[3] = 30; // C5M--3H5M
  AA[p].atom[19].nbonded_neighbors = 3; // C6
  AA[p].atom[19].bonded_neighbor[0] = 11; // C6--N1
  AA[p].atom[19].bonded_neighbor[1] = 17; // C6--C5
  AA[p].atom[19].bonded_neighbor[2] = 31; // C6--H6
  AA[p].atom[20].nbonded_neighbors = 1; //1H5*
  AA[p].atom[20].bonded_neighbor[0] = 4; //1H5*--C5*
  AA[p].atom[21].nbonded_neighbors = 1; //2H5*
  AA[p].atom[21].bonded_neighbor[0] = 4; //2H5*--C5*
  AA[p].atom[22].nbonded_neighbors = 1; // H4*
  AA[p].atom[22].bonded_neighbor[0] = 5; // H4*--C4*
  AA[p].atom[23].nbonded_neighbors = 1; // H3*
  AA[p].atom[23].bonded_neighbor[0] = 7; // H3*--C3*
  AA[p].atom[24].nbonded_neighbors = 1; //1H2*
  AA[p].atom[24].bonded_neighbor[0] = 9; //1H2*--C2*
  AA[p].atom[25].nbonded_neighbors = 1; //2H2*
  AA[p].atom[25].bonded_neighbor[0] = 9; //2H2*--C2*
  AA[p].atom[26].nbonded_neighbors = 1; // H1*
  AA[p].atom[26].bonded_neighbor[0] = 10; // H1*--C1*
  AA[p].atom[27].nbonded_neighbors = 1; // H3
  AA[p].atom[27].bonded_neighbor[0] = 14; // H3--N3
  AA[p].atom[28].nbonded_neighbors = 1; //1H5M
  AA[p].atom[28].bonded_neighbor[0] = 18; //1H5M--C5M
  AA[p].atom[29].nbonded_neighbors = 1; //2H5M
  AA[p].atom[29].bonded_neighbor[0] = 18; //2H5M--C5M
  AA[p].atom[30].nbonded_neighbors = 1; //3H5M
  AA[p].atom[30].bonded_neighbor[0] = 18; //3H5M--C5M
  AA[p].atom[31].nbonded_neighbors = 1; // H6
  AA[p].atom[31].bonded_neighbor[0] = 19; // H6--C6

  for ( k = 0; k < 3; k++ )
    AA[p].atom[k].hastemplate = true;

  for ( k = 3; k < 20; k++ )
    AA[p].atom[k].hastemplate = false;

  for ( k = 20; k < 32; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template atoms used for placing atoms THY

  //bk   template for building P
  AA[p].atom[0].ta[0] = 3; //   O5*
  AA[p].atom[0].ta[1] = 4; //   C5*
  AA[p].atom[0].ta[2] = 5; //   C4*

  //bk   template for building O1P
  AA[p].atom[1].ta[0] = 0; //   P
  AA[p].atom[1].ta[1] = 4; //   C5*
  AA[p].atom[1].ta[2] = 3; //   O5*

  //bk   template for building  O2P
  AA[p].atom[2].ta[0] = 0; //   P
  AA[p].atom[2].ta[1] = 1; //   O1P
  AA[p].atom[2].ta[2] = 3; //   O5*

  //bk   template for building 1H5*
  AA[p].atom[20].ta[0] = 4; //   C5*
  AA[p].atom[20].ta[1] = 3; //   O5*
  AA[p].atom[20].ta[2] = 5; //   C4*

  //bk   template for building 2H5*
  AA[p].atom[21].ta[0] = 4; //   C5*
  AA[p].atom[21].ta[1] = 3; //   O5*
  AA[p].atom[21].ta[2] = 5; //   C4*

  //bk   template for building H4*
  AA[p].atom[22].ta[0] = 5; //   C4*
  AA[p].atom[22].ta[1] = 6; //   O4*
  AA[p].atom[22].ta[2] = 4; //   C5*

  //bk   template for building H3*
  AA[p].atom[23].ta[0] = 7; //   C3*
  AA[p].atom[23].ta[1] = 9; //   C2*
  AA[p].atom[23].ta[2] = 5; //   C4*

  //bk   template for building 1H2*
  AA[p].atom[24].ta[0] = 9; //   C2*
  AA[p].atom[24].ta[1] = 7; //   C3*
  AA[p].atom[24].ta[2] = 10; //   C1*

  //bk   template for building 2H2*
  AA[p].atom[25].ta[0] = 9; //   C2*
  AA[p].atom[25].ta[1] = 7; //   C3*
  AA[p].atom[25].ta[2] = 10; //   C1*

  //bk   template for building H1*
  AA[p].atom[26].ta[0] = 10; //   C1*
  AA[p].atom[26].ta[1] = 9; //   C2*
  AA[p].atom[26].ta[2] = 6; //   O4*

  //bk   template for building H3
  AA[p].atom[27].ta[0] = 14; //   N3
  AA[p].atom[27].ta[1] = 12; //   C2
  AA[p].atom[27].ta[2] = 15; //   C4

  //bk   template for building 1H5M
  AA[p].atom[28].ta[0] = 18; //   C5M
  AA[p].atom[28].ta[1] = 17; //   C5
  AA[p].atom[28].ta[2] = 15; //   C4

  //bk   template for building 2H5M
  AA[p].atom[29].ta[0] = 18; //   C5M
  AA[p].atom[29].ta[1] = 17; //   C5
  AA[p].atom[29].ta[2] = 15; //   C4

  //bk   template for building 3H5M
  AA[p].atom[30].ta[0] = 18; //   C5M
  AA[p].atom[30].ta[1] = 17; //   C5
  AA[p].atom[30].ta[2] = 15; //   C4

  //bk   template for building H6
  AA[p].atom[31].ta[0] = 19; //   C6
  AA[p].atom[31].ta[1] = 17; //   C5
  AA[p].atom[31].ta[2] = 15; //   C4


  //bk   template coordinates for the amino acid
  //jjh new template from Kevin Karplus

  AA[p].atom[0].icoor[0] = 16.9890; //   P
  AA[p].atom[0].icoor[1] = 16.3550; //   P
  AA[p].atom[0].icoor[2] = 15.5690; //   P
  AA[p].atom[1].icoor[0] = 17.1860; //   O1P
  AA[p].atom[1].icoor[1] = 14.8960; //   O1P
  AA[p].atom[1].icoor[2] = 15.7680; //   O1P
  AA[p].atom[2].icoor[0] = 17.9460; //   O2P
  AA[p].atom[2].icoor[1] = 17.1130; //   O2P
  AA[p].atom[2].icoor[2] = 14.7250; //   O2P
  AA[p].atom[3].icoor[0] = 15.5280; //   O5*
  AA[p].atom[3].icoor[1] = 16.6070; //   O5*
  AA[p].atom[3].icoor[2] = 14.9980; //   O5*
  AA[p].atom[4].icoor[0] = 14.3820; //   C5*
  AA[p].atom[4].icoor[1] = 16.2560; //   C5*
  AA[p].atom[4].icoor[2] = 15.7660; //   C5*
  AA[p].atom[5].icoor[0] = 13.1310; //   C4*
  AA[p].atom[5].icoor[1] = 16.7580; //   C4*
  AA[p].atom[5].icoor[2] = 15.0890; //   C4*
  AA[p].atom[6].icoor[0] = 13.1040; //   O4*
  AA[p].atom[6].icoor[1] = 18.2040; //   O4*
  AA[p].atom[6].icoor[2] = 15.1450; //   O4*
  AA[p].atom[7].icoor[0] = 13.0500; //   C3*
  AA[p].atom[7].icoor[1] = 16.3880; //   C3*
  AA[p].atom[7].icoor[2] = 13.6110; //   C3*
  AA[p].atom[8].icoor[0] = 11.6910; //   O3*
  AA[p].atom[8].icoor[1] = 16.1090; //   O3*
  AA[p].atom[8].icoor[2] = 13.2590; //   O3*
  AA[p].atom[9].icoor[0] = 13.5600; //   C2*
  AA[p].atom[9].icoor[1] = 17.6380; //   C2*
  AA[p].atom[9].icoor[2] = 12.9160; //   C2*
  AA[p].atom[10].icoor[0] = 13.0650; //   C1*
  AA[p].atom[10].icoor[1] = 18.7410; //   C1*
  AA[p].atom[10].icoor[2] = 13.8350; //   C1*
  AA[p].atom[11].icoor[0] = 13.8560; //   N1
  AA[p].atom[11].icoor[1] = 19.9940; //   N1
  AA[p].atom[11].icoor[2] = 13.8270; //   N1
  AA[p].atom[12].icoor[0] = 13.2050; //   C2
  AA[p].atom[12].icoor[1] = 21.1310; //   C2
  AA[p].atom[12].icoor[2] = 14.1060; //   C2
  AA[p].atom[13].icoor[0] = 12.0200; //   O2
  AA[p].atom[13].icoor[1] = 21.1840; //   O2
  AA[p].atom[13].icoor[2] = 14.4140; //   O2
  AA[p].atom[14].icoor[0] = 13.9830; //   N3
  AA[p].atom[14].icoor[1] = 22.2190; //   N3
  AA[p].atom[14].icoor[2] = 14.0330; //   N3
  AA[p].atom[15].icoor[0] = 15.2960; //   C4
  AA[p].atom[15].icoor[1] = 22.3150; //   C4
  AA[p].atom[15].icoor[2] = 13.7480; //   C4
  AA[p].atom[16].icoor[0] = 15.8300; //   O4
  AA[p].atom[16].icoor[1] = 23.4230; //   O4
  AA[p].atom[16].icoor[2] = 13.6760; //   O4
  AA[p].atom[17].icoor[0] = 15.9400; //   C5
  AA[p].atom[17].icoor[1] = 21.0810; //   C5
  AA[p].atom[17].icoor[2] = 13.5940; //   C5
  AA[p].atom[18].icoor[0] = 17.3820; //   C5M
  AA[p].atom[18].icoor[1] = 21.0400; //   C5M
  AA[p].atom[18].icoor[2] = 13.2180; //   C5M
  AA[p].atom[19].icoor[0] = 15.1870; //   C6
  AA[p].atom[19].icoor[1] = 19.9950; //   C6
  AA[p].atom[19].icoor[2] = 13.6210; //   C6
  AA[p].atom[20].icoor[0] = 14.4670; //  1H5*
  AA[p].atom[20].icoor[1] = 16.7100; //  1H5*
  AA[p].atom[20].icoor[2] = 16.7770; //  1H5*
  AA[p].atom[21].icoor[0] = 14.3340; //  2H5*
  AA[p].atom[21].icoor[1] = 15.1510; //  2H5*
  AA[p].atom[21].icoor[2] = 15.8660; //  2H5*
  AA[p].atom[22].icoor[0] = 12.2280; //   H4*
  AA[p].atom[22].icoor[1] = 16.3730; //   H4*
  AA[p].atom[22].icoor[2] = 15.6100; //   H4*
  AA[p].atom[23].icoor[0] = 13.6960; //   H3*
  AA[p].atom[23].icoor[1] = 15.5100; //   H3*
  AA[p].atom[23].icoor[2] = 13.4000; //   H3*
  AA[p].atom[24].icoor[0] = 14.6680; //  1H2*
  AA[p].atom[24].icoor[1] = 17.6380; //  1H2*
  AA[p].atom[24].icoor[2] = 12.8430; //  1H2*
  AA[p].atom[25].icoor[0] = 13.1320; //  2H2*
  AA[p].atom[25].icoor[1] = 17.7380; //  2H2*
  AA[p].atom[25].icoor[2] = 11.8960; //  2H2*
  AA[p].atom[26].icoor[0] = 12.0200; //   H1*
  AA[p].atom[26].icoor[1] = 19.0270; //   H1*
  AA[p].atom[26].icoor[2] = 13.5910; //   H1*
  AA[p].atom[27].icoor[0] = 13.4880; //   H3
  AA[p].atom[27].icoor[1] = 23.0920; //   H3
  AA[p].atom[27].icoor[2] = 14.1350; //   H3
  AA[p].atom[28].icoor[0] = 17.8460; //  1H5M
  AA[p].atom[28].icoor[1] = 20.0790; //  1H5M
  AA[p].atom[28].icoor[2] = 13.5380; //  1H5M
  AA[p].atom[29].icoor[0] = 17.9520; //  2H5M
  AA[p].atom[29].icoor[1] = 21.8710; //  2H5M
  AA[p].atom[29].icoor[2] = 13.6860; //  2H5M
  AA[p].atom[30].icoor[0] = 17.4950; //  3H5M
  AA[p].atom[30].icoor[1] = 21.1200; //  3H5M
  AA[p].atom[30].icoor[2] = 12.1130; //  3H5M
  AA[p].atom[31].icoor[0] = 15.6760; //   H6
  AA[p].atom[31].icoor[1] = 19.0160; //   H6
  AA[p].atom[31].icoor[2] = 13.6610; //   H6


  //FIN DE DTHY


  //rGUA
//  p = RGU;
  p = GUA;

//  strcpy( AA[p].aa_name3, "  rG" );
  strcpy( AA[p].aa_name3, "  G" );
//  AA[p].aa_name1 = 'X';
  AA[p].aa_name1 = 'g';

  AA[p].mass = 0;

  AA[p].is_protein = false;
  AA[p].is_RNA = true;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = true;
  AA[p].aa_is_nonpolar = false;
  AA[p].aa_is_aromatic = true;
  AA[p].aa_is_charged = true;

  AA[p].natoms = 34;
  AA[p].nheavyatoms = 23;
  AA[p].nchi = 1;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );
  AA[p].atom = ( t_aa_atom_p * ) malloc( AA[p].natoms * sizeof( t_aa_atom_p ) );


  strcpy( AA[p].atom[0].atom_name, " P  " );
//  strcpy( AA[p].atom[1].atom_name, " O1P" );
//  strcpy( AA[p].atom[2].atom_name, " O2P" );

  strcpy( AA[p].atom[1].atom_name, " OP1" );
  strcpy( AA[p].atom[2].atom_name, " OP2" );

//  strcpy( AA[p].atom[3].atom_name, " O5*" );
//  strcpy( AA[p].atom[4].atom_name, " C5*" );
//  strcpy( AA[p].atom[5].atom_name, " C4*" );
//  strcpy( AA[p].atom[6].atom_name, " O4*" );
//  strcpy( AA[p].atom[7].atom_name, " C3*" );
//  strcpy( AA[p].atom[8].atom_name, " O3*" );
//  strcpy( AA[p].atom[9].atom_name, " C2*" );
//  strcpy( AA[p].atom[10].atom_name, " O2*" );
//  strcpy( AA[p].atom[11].atom_name, " C1*" );

  strcpy( AA[p].atom[3].atom_name, " O5'" );
  strcpy( AA[p].atom[4].atom_name, " C5'" );
  strcpy( AA[p].atom[5].atom_name, " C4'" );
  strcpy( AA[p].atom[6].atom_name, " O4'" );
  strcpy( AA[p].atom[7].atom_name, " C3'" );
  strcpy( AA[p].atom[8].atom_name, " O3'" );
  strcpy( AA[p].atom[9].atom_name, " C2'" );
  strcpy( AA[p].atom[10].atom_name, " O2'" );
  strcpy( AA[p].atom[11].atom_name, " C1'" );

  strcpy( AA[p].atom[12].atom_name, " N1 " );
  strcpy( AA[p].atom[13].atom_name, " C2 " );
  strcpy( AA[p].atom[14].atom_name, " N2 " );
  strcpy( AA[p].atom[15].atom_name, " N3 " );
  strcpy( AA[p].atom[16].atom_name, " C4 " );
  strcpy( AA[p].atom[17].atom_name, " C5 " );
  strcpy( AA[p].atom[18].atom_name, " C6 " );
  strcpy( AA[p].atom[19].atom_name, " O6 " );
  strcpy( AA[p].atom[20].atom_name, " N7 " );
  strcpy( AA[p].atom[21].atom_name, " C8 " );
  strcpy( AA[p].atom[22].atom_name, " N9 " );

//  strcpy( AA[p].atom[23].atom_name, "1H5*" );
//  strcpy( AA[p].atom[24].atom_name, "2H5*" );
//  strcpy( AA[p].atom[25].atom_name, " H4*" );
//  strcpy( AA[p].atom[26].atom_name, " H3*" );
//  strcpy( AA[p].atom[27].atom_name, "1H2*" );
//  strcpy( AA[p].atom[28].atom_name, " H1*" );
//  strcpy( AA[p].atom[29].atom_name, "2HO*" );

  strcpy( AA[p].atom[23].atom_name, "1H5'" );
  strcpy( AA[p].atom[24].atom_name, "2H5'" );
  strcpy( AA[p].atom[25].atom_name, " H4'" );
  strcpy( AA[p].atom[26].atom_name, " H3'" );
  strcpy( AA[p].atom[27].atom_name, "1H2'" );
  strcpy( AA[p].atom[28].atom_name, " H1'" );
  strcpy( AA[p].atom[29].atom_name, "2HO'" );

  strcpy( AA[p].atom[30].atom_name, " H1 " );
  strcpy( AA[p].atom[31].atom_name, "1H2 " );
  strcpy( AA[p].atom[32].atom_name, "2H2 " );
  strcpy( AA[p].atom[33].atom_name, " H8 " );


  AA[p].atom[0].charge = 1.500; // P
  AA[p].atom[1].charge = -.780; // O1P
  AA[p].atom[2].charge = -.780; // O2P
  AA[p].atom[3].charge = -.570; // O5*
  AA[p].atom[4].charge = -.080; // C5*
  AA[p].atom[5].charge =.160; // C4*
  AA[p].atom[6].charge = -.500; // O4*
  AA[p].atom[7].charge =.010; // C3*
  AA[p].atom[8].charge = -.570; // O3*
  AA[p].atom[9].charge =.140; // C2*
  AA[p].atom[10].charge = -.660; // O2*
  AA[p].atom[11].charge =.160; // C1*
  AA[p].atom[12].charge = -.340; // N1
  AA[p].atom[13].charge =.750; // C2
  AA[p].atom[14].charge = -.680; // N2
  AA[p].atom[15].charge = -.740; // N3
  AA[p].atom[16].charge =.260; // C4
  AA[p].atom[17].charge =.000; // C5
  AA[p].atom[18].charge =.540; // C6
  AA[p].atom[19].charge = -.510; // O6
  AA[p].atom[20].charge = -.600; // N7
  AA[p].atom[21].charge =.250; // C8
  AA[p].atom[22].charge = -.020; // N9
  AA[p].atom[23].charge =.090; // 1H5*
  AA[p].atom[24].charge =.090; // 2H5*
  AA[p].atom[25].charge =.090; // H4*
  AA[p].atom[26].charge =.090; // H3*
  AA[p].atom[27].charge =.090; // 1H2*
  AA[p].atom[28].charge =.090; // H1*
  AA[p].atom[29].charge =.430; // 2HO*
  AA[p].atom[30].charge =.260; // H1
  AA[p].atom[31].charge =.320; // 1H2
  AA[p].atom[32].charge =.350; // 2H2
  AA[p].atom[33].charge =.160; // H8

  switch(opt)
  {
  case Sybil:
	  AA[p].atom[0].fullatom_type = 12; // Phos   P    P
	  AA[p].atom[1].fullatom_type = 10; // 69; // OOC    O1P  ON3
	  AA[p].atom[2].fullatom_type = 10; // 69; // OOC    O2P  ON3
	  AA[p].atom[3].fullatom_type = 9; // 68; // ONH2   O5*  ON2
	  AA[p].atom[4].fullatom_type = 2; // 61; // CH2    C5*  CN8 (CN8B)
	  AA[p].atom[5].fullatom_type = 2; // 60; // CH1    C4*  CN7
	  AA[p].atom[6].fullatom_type = 9; // 72; // OH     O4*  ON6 (ON6B)
	  AA[p].atom[7].fullatom_type = 2; // 60; // CH1    C3*  CN7
	  AA[p].atom[8].fullatom_type = 9; // 68; // ONH2   O3*  ON2
	  AA[p].atom[9].fullatom_type = 2; // 60; // CH2    C2*  CN7 (CN7B)
	  AA[p].atom[10].fullatom_type = 9; // 71; // OH     O2*  ON5
	  AA[p].atom[11].fullatom_type = 2; // 60; // CH1    C1*  CN7 (CN7B)
	  AA[p].atom[12].fullatom_type = 5; // 64; // Ntrp   N1   NN2 (NN2G)
	  AA[p].atom[13].fullatom_type = 3; // 56; // aroC   C2   CN2
	  AA[p].atom[14].fullatom_type = 5; // 63; // NH2O   N2   NN1
	  AA[p].atom[15].fullatom_type = 5; // 65; // Nhis   N3   NN3 (NN3G)
	  AA[p].atom[16].fullatom_type = 3; // 59; // aroC   C4   CN5
	  AA[p].atom[17].fullatom_type = 3; // 59; // aroC   C5   CN5 (CN5G)
	  AA[p].atom[18].fullatom_type = 1; // 55; // CObb   C6   CN1
	  AA[p].atom[19].fullatom_type = 8; // 67; // OCbb   O6   ON1
	  AA[p].atom[20].fullatom_type = 5; // 66; // Nhis   N7   NN4
	  AA[p].atom[21].fullatom_type = 3; // 58; // aroC   C8   CN4
	  AA[p].atom[22].fullatom_type = 5; // 64; // NHis   N9   NN2 (NN2B)
	  AA[p].atom[23].fullatom_type = 14; // 74; // Hapo  1H5*  HNP
	  AA[p].atom[24].fullatom_type = 14; // 74; // Hapo  2H5*  HNP
	  AA[p].atom[25].fullatom_type = 14; // 74; // Hapo   H4*  HNP
	  AA[p].atom[26].fullatom_type = 14; // 74; // Hapo   H3*  HNP
	  AA[p].atom[27].fullatom_type = 14; // 74; // Hapo  1H2*  HNP
	  AA[p].atom[28].fullatom_type = 14; // 74; // Hapo   H1*  HNP
	  AA[p].atom[29].fullatom_type = 13; // 73; // Hpol  2HO*  HP
	  AA[p].atom[30].fullatom_type = 13; // 73; // Hpol   H1   HP
	  AA[p].atom[31].fullatom_type = 13; // 73; // Hpol  1H2   HP
	  AA[p].atom[32].fullatom_type = 13; // 73; // Hpol  2H2   HP
	  AA[p].atom[33].fullatom_type = 15; // 74; // Haro   H8   HNP
	  break;

  default:
	  AA[p].atom[0].fullatom_type = 21; // Phos   P    P
	  AA[p].atom[1].fullatom_type = 15; // 69; // OOC    O1P  ON3
	  AA[p].atom[2].fullatom_type = 15; // 69; // OOC    O2P  ON3
	  AA[p].atom[3].fullatom_type = 14; // 68; // ONH2   O5*  ON2
	  AA[p].atom[4].fullatom_type = 4; // 61; // CH2    C5*  CN8 (CN8B)
	  AA[p].atom[5].fullatom_type = 3; // 60; // CH1    C4*  CN7
	  AA[p].atom[6].fullatom_type = 13; // 72; // OH     O4*  ON6 (ON6B)
	  AA[p].atom[7].fullatom_type = 3; // 60; // CH1    C3*  CN7
	  AA[p].atom[8].fullatom_type = 14; // 68; // ONH2   O3*  ON2
	  AA[p].atom[9].fullatom_type = 4; // 60; // CH2    C2*  CN7 (CN7B)
	  AA[p].atom[10].fullatom_type = 13; // 71; // OH     O2*  ON5
	  AA[p].atom[11].fullatom_type = 3; // 60; // CH1    C1*  CN7 (CN7B)
	  AA[p].atom[12].fullatom_type = 7; // 64; // Ntrp   N1   NN2 (NN2G)
	  AA[p].atom[13].fullatom_type = 6; // 56; // aroC   C2   CN2
	  AA[p].atom[14].fullatom_type = 9; // 63; // NH2O   N2   NN1
	  AA[p].atom[15].fullatom_type = 8; // 65; // Nhis   N3   NN3 (NN3G)
	  AA[p].atom[16].fullatom_type = 6; // 59; // aroC   C4   CN5
	  AA[p].atom[17].fullatom_type = 6; // 59; // aroC   C5   CN5 (CN5G)
	  AA[p].atom[18].fullatom_type = 19; // 55; // CObb   C6   CN1
	  AA[p].atom[19].fullatom_type = 20; // 67; // OCbb   O6   ON1
	  AA[p].atom[20].fullatom_type = 8; // 66; // Nhis   N7   NN4
	  AA[p].atom[21].fullatom_type = 6; // 58; // aroC   C8   CN4
	  AA[p].atom[22].fullatom_type = 8; // 64; // NHis   N9   NN2 (NN2B)
	  AA[p].atom[23].fullatom_type = 23; // 74; // Hapo  1H5*  HNP
	  AA[p].atom[24].fullatom_type = 23; // 74; // Hapo  2H5*  HNP
	  AA[p].atom[25].fullatom_type = 23; // 74; // Hapo   H4*  HNP
	  AA[p].atom[26].fullatom_type = 23; // 74; // Hapo   H3*  HNP
	  AA[p].atom[27].fullatom_type = 23; // 74; // Hapo  1H2*  HNP
	  AA[p].atom[28].fullatom_type = 23; // 74; // Hapo   H1*  HNP
	  AA[p].atom[29].fullatom_type = 22; // 73; // Hpol  2HO*  HP
	  AA[p].atom[30].fullatom_type = 22; // 73; // Hpol   H1   HP
	  AA[p].atom[31].fullatom_type = 22; // 73; // Hpol  1H2   HP
	  AA[p].atom[32].fullatom_type = 22; // 73; // Hpol  2H2   HP
	  AA[p].atom[33].fullatom_type = 24; // 74; // Haro   H8   HNP
	  break;
  }
  //jjh intra residue bonding

  AA[p].atom[0].nbonded_neighbors = 3; // P
  AA[p].atom[0].bonded_neighbor[0] = 1; // P--O1P
  AA[p].atom[0].bonded_neighbor[1] = 2; // P--O2P
  AA[p].atom[0].bonded_neighbor[2] = 3; // P--O5*
  AA[p].atom[1].nbonded_neighbors = 1; // O1P
  AA[p].atom[1].bonded_neighbor[0] = 0; // O1P--P
  AA[p].atom[2].nbonded_neighbors = 1; // O2P
  AA[p].atom[2].bonded_neighbor[0] = 0; // O2P--P
  AA[p].atom[3].nbonded_neighbors = 2; // O5*
  AA[p].atom[3].bonded_neighbor[0] = 0; // O5*--P
  AA[p].atom[3].bonded_neighbor[1] = 4; // O5*--C5*
  AA[p].atom[4].nbonded_neighbors = 4; // C5*
  AA[p].atom[4].bonded_neighbor[0] = 3; // C5*--O5*
  AA[p].atom[4].bonded_neighbor[1] = 5; // C5*--C4*
  AA[p].atom[4].bonded_neighbor[2] = 23; // C5*--1H5*
  AA[p].atom[4].bonded_neighbor[3] = 24; // C5*--2H5*
  AA[p].atom[5].nbonded_neighbors = 4; // C4*
  AA[p].atom[5].bonded_neighbor[0] = 4; // C4*--C5*
  AA[p].atom[5].bonded_neighbor[1] = 6; // C4*--O4*
  AA[p].atom[5].bonded_neighbor[2] = 7; // C4*--C3*
  AA[p].atom[5].bonded_neighbor[3] = 25; // C4*--H4*
  AA[p].atom[6].nbonded_neighbors = 2; // O4*
  AA[p].atom[6].bonded_neighbor[0] = 5; // O4*--C4*
  AA[p].atom[6].bonded_neighbor[1] = 11; // O4*--C1*
  AA[p].atom[7].nbonded_neighbors = 4; // C3*
  AA[p].atom[7].bonded_neighbor[0] = 5; // C3*--C4*
  AA[p].atom[7].bonded_neighbor[1] = 8; // C3*--O3*
  AA[p].atom[7].bonded_neighbor[2] = 9; // C3*--C2*
  AA[p].atom[7].bonded_neighbor[3] = 26; // C3*--H3*
  AA[p].atom[8].nbonded_neighbors = 1; // O3*
  AA[p].atom[8].bonded_neighbor[0] = 7; // O3*--C3*
  AA[p].atom[9].nbonded_neighbors = 4; // C2*
  AA[p].atom[9].bonded_neighbor[0] = 7; // C2*--C3*
  AA[p].atom[9].bonded_neighbor[1] = 10; // C2*--O2*
  AA[p].atom[9].bonded_neighbor[2] = 11; // C2*--C1*
  AA[p].atom[9].bonded_neighbor[3] = 27; // C2*--1H2*
  AA[p].atom[10].nbonded_neighbors = 2; // O2*
  AA[p].atom[10].bonded_neighbor[0] = 9; // O2*--C2*
  AA[p].atom[10].bonded_neighbor[1] = 29; // 02*--2HO*
  AA[p].atom[11].nbonded_neighbors = 4; // C1*
  AA[p].atom[11].bonded_neighbor[0] = 6; // C1*--O4*
  AA[p].atom[11].bonded_neighbor[1] = 9; // C1*--C2*
  AA[p].atom[11].bonded_neighbor[2] = 22; // C1*--N9
  AA[p].atom[11].bonded_neighbor[3] = 28; // C1*--H1*
  AA[p].atom[12].nbonded_neighbors = 3; // N1
  AA[p].atom[12].bonded_neighbor[0] = 13; // N1--C2
  AA[p].atom[12].bonded_neighbor[1] = 18; // N1--C6
  AA[p].atom[12].bonded_neighbor[2] = 30; // N1--H1
  AA[p].atom[13].nbonded_neighbors = 3; // C2
  AA[p].atom[13].bonded_neighbor[0] = 12; // C2--N1
  AA[p].atom[13].bonded_neighbor[1] = 14; // C2--N2
  AA[p].atom[13].bonded_neighbor[2] = 15; // C2--N3
  AA[p].atom[14].nbonded_neighbors = 3; // N2
  AA[p].atom[14].bonded_neighbor[0] = 13; // N2--C2
  AA[p].atom[14].bonded_neighbor[1] = 31; // N2--1H2
  AA[p].atom[14].bonded_neighbor[2] = 32; // N2--2H2
  AA[p].atom[15].nbonded_neighbors = 2; // N3
  AA[p].atom[15].bonded_neighbor[0] = 13; // N3--C2
  AA[p].atom[15].bonded_neighbor[1] = 16; // N3--C4
  AA[p].atom[16].nbonded_neighbors = 3; // C4
  AA[p].atom[16].bonded_neighbor[0] = 15; // C4--N3
  AA[p].atom[16].bonded_neighbor[1] = 17; // C4--C5
  AA[p].atom[16].bonded_neighbor[2] = 22; // C4--N9
  AA[p].atom[17].nbonded_neighbors = 3; // C5
  AA[p].atom[17].bonded_neighbor[0] = 16; // C5--C4
  AA[p].atom[17].bonded_neighbor[1] = 18; // C5--C6
  AA[p].atom[17].bonded_neighbor[2] = 20; // C5--N7
  AA[p].atom[18].nbonded_neighbors = 3; // C6
  AA[p].atom[18].bonded_neighbor[0] = 12; // C6--N1
  AA[p].atom[18].bonded_neighbor[1] = 17; // C6--C5
  AA[p].atom[18].bonded_neighbor[2] = 19; // C6--O6
  AA[p].atom[19].nbonded_neighbors = 1; // O6
  AA[p].atom[19].bonded_neighbor[0] = 18; // O6--C6
  AA[p].atom[20].nbonded_neighbors = 2; // N7
  AA[p].atom[20].bonded_neighbor[0] = 17; // N7--C5
  AA[p].atom[20].bonded_neighbor[1] = 21; // N7--C8
  AA[p].atom[21].nbonded_neighbors = 3; //C8
  AA[p].atom[21].bonded_neighbor[0] = 20; // C8--N7
  AA[p].atom[21].bonded_neighbor[1] = 22; // C8--N9
  AA[p].atom[21].bonded_neighbor[2] = 33; // C8--H8
  AA[p].atom[22].nbonded_neighbors = 3; //N9
  AA[p].atom[22].bonded_neighbor[0] = 11; // N9--C1*
  AA[p].atom[22].bonded_neighbor[1] = 16; // N9--C4
  AA[p].atom[22].bonded_neighbor[2] = 21; // N9--C8
  AA[p].atom[23].nbonded_neighbors = 1; // 1H5*
  AA[p].atom[23].bonded_neighbor[0] = 4; //1H5*--C5*
  AA[p].atom[24].nbonded_neighbors = 1; //2H5*
  AA[p].atom[24].bonded_neighbor[0] = 4; //2H5*--C5*
  AA[p].atom[25].nbonded_neighbors = 1; // H4*
  AA[p].atom[25].bonded_neighbor[0] = 5; // H4*--C4*
  AA[p].atom[26].nbonded_neighbors = 1; // H3*
  AA[p].atom[26].bonded_neighbor[0] = 7; // H3*--C3*
  AA[p].atom[27].nbonded_neighbors = 1; //1H2*
  AA[p].atom[27].bonded_neighbor[0] = 9; //1H2*--C2*
  AA[p].atom[28].nbonded_neighbors = 1; // H1*
  AA[p].atom[28].bonded_neighbor[0] = 11; // H1*--C1*
  AA[p].atom[29].nbonded_neighbors = 1; //2HO*
  AA[p].atom[29].bonded_neighbor[0] = 10; //2HO*--O2*
  AA[p].atom[30].nbonded_neighbors = 1; // H1
  AA[p].atom[30].bonded_neighbor[0] = 12; // H1--N1
  AA[p].atom[31].nbonded_neighbors = 1; //1H2
  AA[p].atom[31].bonded_neighbor[0] = 14; //!H2--N2
  AA[p].atom[32].nbonded_neighbors = 1; //2H2
  AA[p].atom[32].bonded_neighbor[0] = 14; //2H2--N2
  AA[p].atom[33].nbonded_neighbors = 1; // H8
  AA[p].atom[33].bonded_neighbor[0] = 21; // H8--C8

  for ( k = 0; k < 3; k++ )
    AA[p].atom[k].hastemplate = true;

  for ( k = 3; k < 23; k++ )
    AA[p].atom[k].hastemplate = false;

  AA[p].atom[10].hastemplate = true;

  for ( k = 23; k < 34; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template atoms used for placing atoms rGUA

  //bk   template for building P
  AA[p].atom[0].ta[0] = 3; //   O5*
  AA[p].atom[0].ta[1] = 4; //   C5*
  AA[p].atom[0].ta[2] = 5; //   C4*

  //bk   template for building O1P
  AA[p].atom[1].ta[0] = 0; //   P
  AA[p].atom[1].ta[1] = 4; //   C5*
  AA[p].atom[1].ta[2] = 3; //   O5*

  //bk   template for building  O2P
  AA[p].atom[2].ta[0] = 0; //   P
  AA[p].atom[2].ta[1] = 1; //   O1P
  AA[p].atom[2].ta[2] = 3; //   O5*

  //bk   template for building  O2*
  AA[p].atom[10].ta[0] = 9; //   C2*
  AA[p].atom[10].ta[1] = 7; //   C3*
  AA[p].atom[10].ta[2] = 11; //   C1*


  //bk   template for building 1H5*
  AA[p].atom[23].ta[0] = 4; //   C5*
  AA[p].atom[23].ta[1] = 3; //   O5*
  AA[p].atom[23].ta[2] = 5; //   C4*

  //bk   template for building 2H5*
  AA[p].atom[24].ta[0] = 4; //   C5*
  AA[p].atom[24].ta[1] = 3; //   O5*
  AA[p].atom[24].ta[2] = 5; //   C4*

  //bk   template for building H4*
  AA[p].atom[25].ta[0] = 5; //   C4*
  AA[p].atom[25].ta[1] = 6; //   O4*
  AA[p].atom[25].ta[2] = 4; //   C5*

  //bk   template for building H3*
  AA[p].atom[26].ta[0] = 7; //   C3*
  AA[p].atom[26].ta[1] = 9; //   C2*
  AA[p].atom[26].ta[2] = 5; //   C4*

  //bk   template for building 1H2*
  AA[p].atom[27].ta[0] = 9; //   C2*
  AA[p].atom[27].ta[1] = 7; //   C3*
  AA[p].atom[27].ta[2] = 11; //   C1*

  //bk   template for building H1*
  AA[p].atom[28].ta[0] = 11; //   C1*
  AA[p].atom[28].ta[1] = 9; //   C2*
  AA[p].atom[28].ta[2] = 6; //   O4*

  //bk   template for building 2HO*
  AA[p].atom[29].ta[0] = 10; //   O2*
  AA[p].atom[29].ta[1] = 9; //   C2*
  AA[p].atom[29].ta[2] = 11; //   C1*


  //bk   template for building H1
  AA[p].atom[30].ta[0] = 12; //   N1
  AA[p].atom[30].ta[1] = 13; //   C2
  AA[p].atom[30].ta[2] = 18; //   C6

  //bk   template for building 1H2
  AA[p].atom[31].ta[0] = 14; //   N2
  AA[p].atom[31].ta[1] = 13; //   C2
  AA[p].atom[31].ta[2] = 15; //   N3

  //bk   template for building 2H2
  AA[p].atom[32].ta[0] = 14; //   N2
  AA[p].atom[32].ta[1] = 13; //   C2
  AA[p].atom[32].ta[2] = 15; //   N3

  //bk   template for building H8
  AA[p].atom[33].ta[0] = 21; //   C8
  AA[p].atom[33].ta[1] = 20; //   N7
  AA[p].atom[33].ta[2] = 22; //   N9


  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 11; //   C1*
  AA[p].chi_atoms[0] [1] = 9; //   C2*
  AA[p].chi_atoms[0] [2] = 10; //   O2*
  AA[p].chi_atoms[0] [3] = 29; //   2HO*



  //bk   template coordinates for the amino acid
  //jjh new template from Kevin Karplus

  AA[p].atom[0].icoor[0] = 58.6140; //   P
  AA[p].atom[0].icoor[1] = 25.3250; //   P
  AA[p].atom[0].icoor[2] = 8.0490; //   P
  AA[p].atom[1].icoor[0] = 58.8580; //   O1P
  AA[p].atom[1].icoor[1] = 25.8950; //   O1P
  AA[p].atom[1].icoor[2] = 6.7060; //   O1P
  AA[p].atom[2].icoor[0] = 58.1010; //   O2P
  AA[p].atom[2].icoor[1] = 23.9420; //   O2P
  AA[p].atom[2].icoor[2] = 8.1820; //   O2P
  AA[p].atom[3].icoor[0] = 59.8330; //   O5*
  AA[p].atom[3].icoor[1] = 25.5970; //   O5*
  AA[p].atom[3].icoor[2] = 9.0270; //   O5*
  AA[p].atom[4].icoor[0] = 59.3530; //   C5*
  AA[p].atom[4].icoor[1] = 25.6870; //   C5*
  AA[p].atom[4].icoor[2] = 10.3700; //   C5*
  AA[p].atom[5].icoor[0] = 60.5470; //   C4*
  AA[p].atom[5].icoor[1] = 25.7860; //   C4*
  AA[p].atom[5].icoor[2] = 11.3010; //   C4*
  AA[p].atom[6].icoor[0] = 60.0670; //   O4*
  AA[p].atom[6].icoor[1] = 25.9590; //   O4*
  AA[p].atom[6].icoor[2] = 12.6500; //   O4*
  AA[p].atom[7].icoor[0] = 61.3310; //   C3*
  AA[p].atom[7].icoor[1] = 24.4800; //   C3*
  AA[p].atom[7].icoor[2] = 11.3660; //   C3*
  AA[p].atom[8].icoor[0] = 62.6840; //   O3*
  AA[p].atom[8].icoor[1] = 24.7240; //   O3*
  AA[p].atom[8].icoor[2] = 11.7700; //   O3*
  AA[p].atom[9].icoor[0] = 60.5970; //   C2*
  AA[p].atom[9].icoor[1] = 23.6920; //   C2*
  AA[p].atom[9].icoor[2] = 12.4050; //   C2*
  AA[p].atom[10].icoor[0] = 61.4780; //   O2*
  AA[p].atom[10].icoor[1] = 22.7840; //   O2*
  AA[p].atom[10].icoor[2] = 13.0020; //   O2*
  AA[p].atom[11].icoor[0] = 60.2380; //   C1*
  AA[p].atom[11].icoor[1] = 24.7610; //   C1*
  AA[p].atom[11].icoor[2] = 13.3810; //   C1*
  AA[p].atom[12].icoor[0] = 58.1080; //   N1
  AA[p].atom[12].icoor[1] = 23.9630; //   N1
  AA[p].atom[12].icoor[2] = 17.8260; //   N1
  AA[p].atom[13].icoor[0] = 59.4380; //   C2
  AA[p].atom[13].icoor[1] = 24.2550; //   C2
  AA[p].atom[13].icoor[2] = 17.5830; //   C2
  AA[p].atom[14].icoor[0] = 60.2420; //   N2
  AA[p].atom[14].icoor[1] = 24.3490; //   N2
  AA[p].atom[14].icoor[2] = 18.6840; //   N2
  AA[p].atom[15].icoor[0] = 59.8800; //   N3
  AA[p].atom[15].icoor[1] = 24.4280; //   N3
  AA[p].atom[15].icoor[2] = 16.3520; //   N3
  AA[p].atom[16].icoor[0] = 58.8950; //   C4
  AA[p].atom[16].icoor[1] = 24.2740; //   C4
  AA[p].atom[16].icoor[2] = 15.4830; //   C4
  AA[p].atom[17].icoor[0] = 57.5750; //   C5
  AA[p].atom[17].icoor[1] = 23.9660; //   C5
  AA[p].atom[17].icoor[2] = 15.6670; //   C5
  AA[p].atom[18].icoor[0] = 57.0700; //   C6
  AA[p].atom[18].icoor[1] = 23.7910; //   C6
  AA[p].atom[18].icoor[2] = 16.9210; //   C6
  AA[p].atom[19].icoor[0] = 55.8920; //   O6
  AA[p].atom[19].icoor[1] = 23.5230; //   O6
  AA[p].atom[19].icoor[2] = 17.1730; //   O6
  AA[p].atom[20].icoor[0] = 56.8980; //   N7
  AA[p].atom[20].icoor[1] = 23.8960; //   N7
  AA[p].atom[20].icoor[2] = 14.4790; //   N7
  AA[p].atom[21].icoor[0] = 57.7940; //   C8
  AA[p].atom[21].icoor[1] = 24.1680; //   C8
  AA[p].atom[21].icoor[2] = 13.5770; //   C8
  AA[p].atom[22].icoor[0] = 59.0270; //   N9
  AA[p].atom[22].icoor[1] = 24.4150; //   N9
  AA[p].atom[22].icoor[2] = 14.1220; //   N9
  AA[p].atom[23].icoor[0] = 58.7270; //  1H5*
  AA[p].atom[23].icoor[1] = 24.8020; //  1H5*
  AA[p].atom[23].icoor[2] = 10.6070; //  1H5*
  AA[p].atom[24].icoor[0] = 58.7220; //  2H5*
  AA[p].atom[24].icoor[1] = 26.5980; //  2H5*
  AA[p].atom[24].icoor[2] = 10.4660; //  2H5*
  AA[p].atom[25].icoor[0] = 61.2140; //   H4*
  AA[p].atom[25].icoor[1] = 26.6310; //   H4*
  AA[p].atom[25].icoor[2] = 11.0270; //   H4*
  AA[p].atom[26].icoor[0] = 61.3260; //   H3*
  AA[p].atom[26].icoor[1] = 23.9690; //   H3*
  AA[p].atom[26].icoor[2] = 10.3790; //   H3*
  AA[p].atom[27].icoor[0] = 59.7060; //  1H2*
  AA[p].atom[27].icoor[1] = 23.1770; //  1H2*
  AA[p].atom[27].icoor[2] = 11.9880; //  1H2*
  AA[p].atom[28].icoor[0] = 61.0920; //   H1*
  AA[p].atom[28].icoor[1] = 24.9440; //   H1*
  AA[p].atom[28].icoor[2] = 14.0720; //   H1*
  AA[p].atom[29].icoor[0] = 61.0020; //  2HO*
  AA[p].atom[29].icoor[1] = 22.2790; //  2HO*
  AA[p].atom[29].icoor[2] = 13.6660; //  2HO*
  AA[p].atom[30].icoor[0] = 57.8380; //   H1
  AA[p].atom[30].icoor[1] = 23.8540; //   H1
  AA[p].atom[30].icoor[2] = 18.7930; //   H1
  AA[p].atom[31].icoor[0] = 59.8770; //  1H2
  AA[p].atom[31].icoor[1] = 24.2300; //  1H2
  AA[p].atom[31].icoor[2] = 19.6070; //  1H2
  AA[p].atom[32].icoor[0] = 61.2220; //  2H2
  AA[p].atom[32].icoor[1] = 24.5280; //  2H2
  AA[p].atom[32].icoor[2] = 18.6030; //  2H2
  AA[p].atom[33].icoor[0] = 57.6080; //   H8
  AA[p].atom[33].icoor[1] = 24.2100; //   H8
  AA[p].atom[33].icoor[2] = 12.5040; //   H8

  //FIN DE rGUA


  //rADE
//  p = RAD;
  p = ADE;

//  strcpy( AA[p].aa_name3, "  rA" );
  strcpy( AA[p].aa_name3, "  A" );
//  AA[p].aa_name1 = 'X';
  AA[p].aa_name1 = 'a';

  AA[p].mass = 0;

  AA[p].is_protein = false;
  AA[p].is_RNA = true;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = true;
  AA[p].aa_is_nonpolar = false;
  AA[p].aa_is_aromatic = true;
  AA[p].aa_is_charged = true;

  AA[p].natoms = 33;
  AA[p].nheavyatoms = 22;
  AA[p].nchi = 1;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );
  AA[p].atom = ( t_aa_atom_p * ) malloc( AA[p].natoms * sizeof( t_aa_atom_p ) );


  strcpy( AA[p].atom[0].atom_name, " P  " );
//  strcpy( AA[p].atom[1].atom_name, " O1P" );
//  strcpy( AA[p].atom[2].atom_name, " O2P" );

  strcpy( AA[p].atom[1].atom_name, " OP1" );
  strcpy( AA[p].atom[2].atom_name, " OP2" );

//  strcpy( AA[p].atom[3].atom_name, " O5*" );
//  strcpy( AA[p].atom[4].atom_name, " C5*" );
//  strcpy( AA[p].atom[5].atom_name, " C4*" );
//  strcpy( AA[p].atom[6].atom_name, " O4*" );
//  strcpy( AA[p].atom[7].atom_name, " C3*" );
//  strcpy( AA[p].atom[8].atom_name, " O3*" );
//  strcpy( AA[p].atom[9].atom_name, " C2*" );
//  strcpy( AA[p].atom[10].atom_name, " O2*" );
//  strcpy( AA[p].atom[11].atom_name, " C1*" );

  strcpy( AA[p].atom[3].atom_name, " O5'" );
  strcpy( AA[p].atom[4].atom_name, " C5'" );
  strcpy( AA[p].atom[5].atom_name, " C4'" );
  strcpy( AA[p].atom[6].atom_name, " O4'" );
  strcpy( AA[p].atom[7].atom_name, " C3'" );
  strcpy( AA[p].atom[8].atom_name, " O3'" );
  strcpy( AA[p].atom[9].atom_name, " C2'" );
  strcpy( AA[p].atom[10].atom_name, " O2'" );
  strcpy( AA[p].atom[11].atom_name, " C1'" );

  strcpy( AA[p].atom[12].atom_name, " N1 " );
  strcpy( AA[p].atom[13].atom_name, " C2 " );
  strcpy( AA[p].atom[14].atom_name, " N3 " );
  strcpy( AA[p].atom[15].atom_name, " C4 " );
  strcpy( AA[p].atom[16].atom_name, " C5 " );
  strcpy( AA[p].atom[17].atom_name, " C6 " );
  strcpy( AA[p].atom[18].atom_name, " N6 " );
  strcpy( AA[p].atom[19].atom_name, " N7 " );
  strcpy( AA[p].atom[20].atom_name, " C8 " );
  strcpy( AA[p].atom[21].atom_name, " N9 " );

//  strcpy( AA[p].atom[22].atom_name, "1H5*" );
//  strcpy( AA[p].atom[23].atom_name, "2H5*" );
//  strcpy( AA[p].atom[24].atom_name, " H4*" );
//  strcpy( AA[p].atom[25].atom_name, " H3*" );
//  strcpy( AA[p].atom[26].atom_name, "1H2*" );
//  strcpy( AA[p].atom[27].atom_name, "2HO*" );
//  strcpy( AA[p].atom[28].atom_name, " H1*" );

  strcpy( AA[p].atom[22].atom_name, "1H5'" );
  strcpy( AA[p].atom[23].atom_name, "2H5'" );
  strcpy( AA[p].atom[24].atom_name, " H4'" );
  strcpy( AA[p].atom[25].atom_name, " H3'" );
  strcpy( AA[p].atom[26].atom_name, "1H2'" );
  strcpy( AA[p].atom[27].atom_name, "2HO'" );
  strcpy( AA[p].atom[28].atom_name, " H1'" );

  strcpy( AA[p].atom[29].atom_name, " H2 " );
  strcpy( AA[p].atom[30].atom_name, "1H6 " );
  strcpy( AA[p].atom[31].atom_name, "2H6 " );
  strcpy( AA[p].atom[32].atom_name, " H8 " );


  AA[p].atom[0].charge = 1.500; // P
  AA[p].atom[1].charge = -.780; // O1P
  AA[p].atom[2].charge = -.780; // O2P
  AA[p].atom[3].charge = -.570; // O5*
  AA[p].atom[4].charge = -.080; // C5*
  AA[p].atom[5].charge =.160; // C4*
  AA[p].atom[6].charge = -.500; // O4*
  AA[p].atom[7].charge =.010; // C3*
  AA[p].atom[8].charge = -.570; // O3*
  AA[p].atom[9].charge =.140; // C2*
  AA[p].atom[10].charge = -.660; // O2*
  AA[p].atom[11].charge =.160; // C1*
  AA[p].atom[12].charge = -.740; // N1
  AA[p].atom[13].charge =.500; // C2
  AA[p].atom[14].charge = -.750; // N3
  AA[p].atom[15].charge =.430; // C4
  AA[p].atom[16].charge =.280; // C5
  AA[p].atom[17].charge =.460; // C6
  AA[p].atom[18].charge = -.770; // N6
  AA[p].atom[19].charge = -.710; // N7
  AA[p].atom[20].charge =.340; // C8
  AA[p].atom[21].charge = -.050; // N9
  AA[p].atom[22].charge =.090; // 1H5*
  AA[p].atom[23].charge =.090; // 2H5*
  AA[p].atom[24].charge =.090; // H4*
  AA[p].atom[25].charge =.090; // H3*
  AA[p].atom[26].charge =.090; // 1H2*
  AA[p].atom[27].charge =.430; // 2HO*
  AA[p].atom[28].charge =.090; // H1*
  AA[p].atom[29].charge =.130; // H2
  AA[p].atom[30].charge =.380; // 1H6
  AA[p].atom[31].charge =.380; // 2H6
  AA[p].atom[32].charge =.120; // H8


  switch(opt)
  {
  case Sybil:
	  AA[p].atom[0].fullatom_type = 12; // Phos   P    P
	  AA[p].atom[1].fullatom_type = 10; // 69; // OOC    O1P  ON3
	  AA[p].atom[2].fullatom_type = 10; // 69; // OOC    O2P  ON3
	  AA[p].atom[3].fullatom_type = 9; // 68; // ONH2   O5*  ON2
	  AA[p].atom[4].fullatom_type = 2; // 61; // CH2    C5*  CN8 (CN8B)
	  AA[p].atom[5].fullatom_type = 2; // 60; // CH1    C4*  CN7
	  AA[p].atom[6].fullatom_type = 9; // 72; // OH     O4*  ON6 (ON6B)
	  AA[p].atom[7].fullatom_type = 2; // 60; // CH1    C3*  CN7
	  AA[p].atom[8].fullatom_type = 9; // 68; // ONH2   O3*  ON2
	  AA[p].atom[9].fullatom_type = 2; // 60; // CH2    C2*  CN7 (CN7B)
	  AA[p].atom[10].fullatom_type = 9; // 71; // OH     O2*  ON5
	  AA[p].atom[11].fullatom_type = 2; // 60; // CH1    C1*  CN7 (CN7B)
	  AA[p].atom[12].fullatom_type = 5; // 65; // Nhis   N1   NN3 (NN3A)
	  AA[p].atom[13].fullatom_type = 3; // 58; // aroC   C2   CN4
	  AA[p].atom[14].fullatom_type = 5; // 65; // Nhis   N3   NN3 (NN3A)
	  AA[p].atom[15].fullatom_type = 3; // 59; // aroC   C4   CN5
	  AA[p].atom[16].fullatom_type = 3; // 59; // aroC   C5   CN5
	  AA[p].atom[17].fullatom_type = 3; // 56; // aroC   C6   CN2
	  AA[p].atom[18].fullatom_type = 5; // 63; // NH2O   N6   NN1
	  AA[p].atom[19].fullatom_type = 5; // 66; // Nhis   N7   NN4
	  AA[p].atom[20].fullatom_type = 3; // 58; // aroC   C8   CN4
	  AA[p].atom[21].fullatom_type = 5; // 64; // Nhis   N9   NN2
	  AA[p].atom[22].fullatom_type = 14; // 74; // Hapo  1H5*  HNP
	  AA[p].atom[23].fullatom_type = 14; // 74; // Hapo  2H5*  HNP
	  AA[p].atom[24].fullatom_type = 14; // 74; // Hapo   H4*  HNP
	  AA[p].atom[25].fullatom_type = 14; // 74; // Hapo   H3*  HNP
	  AA[p].atom[26].fullatom_type = 14; // 74; // Hapo  1H2*  HNP
	  AA[p].atom[27].fullatom_type = 13; // 73; // Hpol  2HO*  HP
	  AA[p].atom[28].fullatom_type = 14; // 74; // Hapo   H1*  HNP
	  AA[p].atom[29].fullatom_type = 15; // 74; // Haro   H2   HNP
	  AA[p].atom[30].fullatom_type = 13; // 73; // Hpol  1H6   HP
	  AA[p].atom[31].fullatom_type = 13; // 73; // Hpol  2H6   HP
	  AA[p].atom[32].fullatom_type = 15; // 74; // Haro   H8   HNP
	  break;
  default:
	  AA[p].atom[0].fullatom_type = 21; // Phos   P    P
	  AA[p].atom[1].fullatom_type = 15; // 69; // OOC    O1P  ON3
	  AA[p].atom[2].fullatom_type = 15; // 69; // OOC    O2P  ON3
	  AA[p].atom[3].fullatom_type = 14; // 68; // ONH2   O5*  ON2
	  AA[p].atom[4].fullatom_type = 4; // 61; // CH2    C5*  CN8 (CN8B)
	  AA[p].atom[5].fullatom_type = 3; // 60; // CH1    C4*  CN7
	  AA[p].atom[6].fullatom_type = 13; // 72; // OH     O4*  ON6 (ON6B)
	  AA[p].atom[7].fullatom_type = 3; // 60; // CH1    C3*  CN7
	  AA[p].atom[8].fullatom_type = 14; // 68; // ONH2   O3*  ON2
	  AA[p].atom[9].fullatom_type = 4; // 60; // CH2    C2*  CN7 (CN7B)
	  AA[p].atom[10].fullatom_type = 13; // 71; // OH     O2*  ON5
	  AA[p].atom[11].fullatom_type = 3; // 60; // CH1    C1*  CN7 (CN7B)
	  AA[p].atom[12].fullatom_type = 8; // 65; // Nhis   N1   NN3 (NN3A)
	  AA[p].atom[13].fullatom_type = 6; // 58; // aroC   C2   CN4
	  AA[p].atom[14].fullatom_type = 8; // 65; // Nhis   N3   NN3 (NN3A)
	  AA[p].atom[15].fullatom_type = 6; // 59; // aroC   C4   CN5
	  AA[p].atom[16].fullatom_type = 6; // 59; // aroC   C5   CN5
	  AA[p].atom[17].fullatom_type = 6; // 56; // aroC   C6   CN2
	  AA[p].atom[18].fullatom_type = 9; // 63; // NH2O   N6   NN1
	  AA[p].atom[19].fullatom_type = 8; // 66; // Nhis   N7   NN4
	  AA[p].atom[20].fullatom_type = 6; // 58; // aroC   C8   CN4
	  AA[p].atom[21].fullatom_type = 8; // 64; // Nhis   N9   NN2
	  AA[p].atom[22].fullatom_type = 23; // 74; // Hapo  1H5*  HNP
	  AA[p].atom[23].fullatom_type = 23; // 74; // Hapo  2H5*  HNP
	  AA[p].atom[24].fullatom_type = 23; // 74; // Hapo   H4*  HNP
	  AA[p].atom[25].fullatom_type = 23; // 74; // Hapo   H3*  HNP
	  AA[p].atom[26].fullatom_type = 23; // 74; // Hapo  1H2*  HNP
	  AA[p].atom[27].fullatom_type = 23; // 73; // Hpol  2HO*  HP
	  AA[p].atom[28].fullatom_type = 23; // 74; // Hapo   H1*  HNP
	  AA[p].atom[29].fullatom_type = 24; // 74; // Haro   H2   HNP
	  AA[p].atom[30].fullatom_type = 22; // 73; // Hpol  1H6   HP
	  AA[p].atom[31].fullatom_type = 22; // 73; // Hpol  2H6   HP
	  AA[p].atom[32].fullatom_type = 24; // 74; // Haro   H8   HNP
	  break;

  }

  //jjh intra residue bonding

  AA[p].atom[0].nbonded_neighbors = 3; // P
  AA[p].atom[0].bonded_neighbor[0] = 1; // P--O1P
  AA[p].atom[0].bonded_neighbor[1] = 2; // P--O2P
  AA[p].atom[0].bonded_neighbor[2] = 3; // P--O5*
  AA[p].atom[1].nbonded_neighbors = 1; // O1P
  AA[p].atom[1].bonded_neighbor[0] = 0; // O1P--P
  AA[p].atom[2].nbonded_neighbors = 1; // O2P
  AA[p].atom[2].bonded_neighbor[0] = 0; // O2P--P
  AA[p].atom[3].nbonded_neighbors = 2; // O5*
  AA[p].atom[3].bonded_neighbor[0] = 0; // O5*--P
  AA[p].atom[3].bonded_neighbor[1] = 4; // O5*--C5*
  AA[p].atom[4].nbonded_neighbors = 4; // C5*
  AA[p].atom[4].bonded_neighbor[0] = 3; // C5*--O5*
  AA[p].atom[4].bonded_neighbor[1] = 5; // C5*--C4*
  AA[p].atom[4].bonded_neighbor[2] = 22; // C5*--1H5*
  AA[p].atom[4].bonded_neighbor[3] = 23; // C5*--2H5*
  AA[p].atom[5].nbonded_neighbors = 4; // C4*
  AA[p].atom[5].bonded_neighbor[0] = 4; // C4*--C5*
  AA[p].atom[5].bonded_neighbor[1] = 6; // C4*--O4*
  AA[p].atom[5].bonded_neighbor[2] = 7; // C4*--C3*
  AA[p].atom[5].bonded_neighbor[3] = 24; // C4*--H4*
  AA[p].atom[6].nbonded_neighbors = 2; // O4*
  AA[p].atom[6].bonded_neighbor[0] = 5; // O4*--C4*
  AA[p].atom[6].bonded_neighbor[1] = 11; // O4*--C1*
  AA[p].atom[7].nbonded_neighbors = 4; // C3*
  AA[p].atom[7].bonded_neighbor[0] = 5; // C3*--C4*
  AA[p].atom[7].bonded_neighbor[1] = 8; // C3*--O3*
  AA[p].atom[7].bonded_neighbor[2] = 9; // C3*--C2*
  AA[p].atom[7].bonded_neighbor[3] = 25; // C3*--H3*
  AA[p].atom[8].nbonded_neighbors = 1; // O3*
  AA[p].atom[8].bonded_neighbor[0] = 7; // O3*--C3*
  AA[p].atom[9].nbonded_neighbors = 4; // C2*
  AA[p].atom[9].bonded_neighbor[0] = 7; // C2*--C3*
  AA[p].atom[9].bonded_neighbor[1] = 10; // C2*--C1*
  AA[p].atom[9].bonded_neighbor[2] = 11; // C2*--1H2*
  AA[p].atom[9].bonded_neighbor[3] = 26; // C2*--2H2*
  AA[p].atom[10].nbonded_neighbors = 2; // 02*
  AA[p].atom[10].bonded_neighbor[0] = 9; // 02*--C2*
  AA[p].atom[10].bonded_neighbor[1] = 27; // O2*--2HO*
  AA[p].atom[11].nbonded_neighbors = 4; // C1*
  AA[p].atom[11].bonded_neighbor[0] = 6; // C1*--O4*
  AA[p].atom[11].bonded_neighbor[1] = 9; // C1*--C2*
  AA[p].atom[11].bonded_neighbor[2] = 21; // C1*--N9
  AA[p].atom[11].bonded_neighbor[3] = 28; // C1*--H1*
  AA[p].atom[12].nbonded_neighbors = 2; // N1
  AA[p].atom[12].bonded_neighbor[0] = 13; // N1--C2
  AA[p].atom[12].bonded_neighbor[1] = 17; // N1--C6
  AA[p].atom[13].nbonded_neighbors = 3; // C2
  AA[p].atom[13].bonded_neighbor[0] = 12; // C2--N1
  AA[p].atom[13].bonded_neighbor[1] = 14; // C2--N3
  AA[p].atom[13].bonded_neighbor[2] = 29; // C2--H2
  AA[p].atom[14].nbonded_neighbors = 2; // N3
  AA[p].atom[14].bonded_neighbor[0] = 13; // N3--C2
  AA[p].atom[14].bonded_neighbor[1] = 15; // N3--C4
  AA[p].atom[15].nbonded_neighbors = 3; // C4
  AA[p].atom[15].bonded_neighbor[0] = 14; // C4--N3
  AA[p].atom[15].bonded_neighbor[1] = 16; // C4--C5
  AA[p].atom[15].bonded_neighbor[2] = 21; // C4--N9
  AA[p].atom[16].nbonded_neighbors = 3; // C5
  AA[p].atom[16].bonded_neighbor[0] = 15; // C5--C4
  AA[p].atom[16].bonded_neighbor[1] = 17; // C5--C6
  AA[p].atom[16].bonded_neighbor[2] = 21; // C5--N7
  AA[p].atom[17].nbonded_neighbors = 3; // C6
  AA[p].atom[17].bonded_neighbor[0] = 12; // C6--N1
  AA[p].atom[17].bonded_neighbor[1] = 16; // C6--C5
  AA[p].atom[17].bonded_neighbor[2] = 18; // C6--N6
  AA[p].atom[18].nbonded_neighbors = 3; // N6
  AA[p].atom[18].bonded_neighbor[0] = 17; // N6--C6
  AA[p].atom[18].bonded_neighbor[1] = 30; // N6--1H6
  AA[p].atom[18].bonded_neighbor[2] = 31; // N6--2H6
  AA[p].atom[19].nbonded_neighbors = 2; // N7
  AA[p].atom[19].bonded_neighbor[0] = 16; // N7--C5
  AA[p].atom[19].bonded_neighbor[1] = 20; // N7--C8
  AA[p].atom[20].nbonded_neighbors = 3; // C8
  AA[p].atom[20].bonded_neighbor[0] = 19; // C8--N7
  AA[p].atom[20].bonded_neighbor[1] = 21; // C8--N9
  AA[p].atom[20].bonded_neighbor[2] = 32; // C8--H8
  AA[p].atom[21].nbonded_neighbors = 3; // N9
  AA[p].atom[21].bonded_neighbor[0] = 11; // N9--C1*
  AA[p].atom[21].bonded_neighbor[1] = 15; // N9--C4
  AA[p].atom[21].bonded_neighbor[2] = 20; // N9--C8
  AA[p].atom[22].nbonded_neighbors = 1; //1H5*
  AA[p].atom[22].bonded_neighbor[0] = 4; //1H5*--C5*
  AA[p].atom[23].nbonded_neighbors = 1; //2H5*
  AA[p].atom[23].bonded_neighbor[0] = 4; //2H5*--C5*
  AA[p].atom[24].nbonded_neighbors = 1; // H4*
  AA[p].atom[24].bonded_neighbor[0] = 5; // H4*--C4*
  AA[p].atom[25].nbonded_neighbors = 1; // H3*
  AA[p].atom[25].bonded_neighbor[0] = 7; // H3*--C3*
  AA[p].atom[26].nbonded_neighbors = 1; //1H2*
  AA[p].atom[26].bonded_neighbor[0] = 9; //1H2*--C2*
  AA[p].atom[27].nbonded_neighbors = 1; //2HO*
  AA[p].atom[27].bonded_neighbor[0] = 10; //2H2*--C2*
  AA[p].atom[28].nbonded_neighbors = 1; // H1*
  AA[p].atom[28].bonded_neighbor[0] = 11; // H1*--C1*
  AA[p].atom[29].nbonded_neighbors = 1; // H2
  AA[p].atom[29].bonded_neighbor[0] = 13; // H2--C2
  AA[p].atom[30].nbonded_neighbors = 1; //1H6
  AA[p].atom[30].bonded_neighbor[0] = 18; //1H6--N6
  AA[p].atom[31].nbonded_neighbors = 1; //2H6
  AA[p].atom[31].bonded_neighbor[0] = 18; //2H6--N6
  AA[p].atom[32].nbonded_neighbors = 1; // H8
  AA[p].atom[32].bonded_neighbor[0] = 20; // H8--C8


  for ( k = 0; k < 3; k++ )
    AA[p].atom[k].hastemplate = true;

  for ( k = 3; k < 22; k++ )
    AA[p].atom[k].hastemplate = false;

  AA[p].atom[10].hastemplate = true;

  for ( k = 22; k < 33; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template atoms used for placing atoms rADE


  //bk   template for building P
  AA[p].atom[0].ta[0] = 3; //   O5*
  AA[p].atom[0].ta[1] = 4; //   C5*
  AA[p].atom[0].ta[2] = 5; //   C4*

  //bk   template for building O1P
  AA[p].atom[1].ta[0] = 0; //   P
  AA[p].atom[1].ta[1] = 4; //   C5*
  AA[p].atom[1].ta[2] = 3; //   O5*

  //bk   template for building  O2P
  AA[p].atom[2].ta[0] = 0; //   P
  AA[p].atom[2].ta[1] = 1; //   O1P
  AA[p].atom[2].ta[2] = 3; //   O5*


  //bk   template for building  O2*
  AA[p].atom[10].ta[0] = 9; //   C2*
  AA[p].atom[10].ta[1] = 7; //   C3*
  AA[p].atom[10].ta[2] = 11; //   C1*


  //bk   template for building 1H5*
  AA[p].atom[22].ta[0] = 4; //   C5*
  AA[p].atom[22].ta[1] = 3; //   O5*
  AA[p].atom[22].ta[2] = 5; //   C4*

  //bk   template for building 2H5*
  AA[p].atom[23].ta[0] = 4; //   C5*
  AA[p].atom[23].ta[1] = 3; //   O5*
  AA[p].atom[23].ta[2] = 5; //   C4*

  //bk   template for building H4*
  AA[p].atom[24].ta[0] = 5; //   C4*
  AA[p].atom[24].ta[1] = 6; //   O4*
  AA[p].atom[24].ta[2] = 4; //   C5*

  //bk   template for building H3*
  AA[p].atom[25].ta[0] = 7; //   C3*
  AA[p].atom[25].ta[1] = 9; //   C2*
  AA[p].atom[25].ta[2] = 5; //   C4*

  //bk   template for building 1H2*
  AA[p].atom[26].ta[0] = 9; //   C2*
  AA[p].atom[26].ta[1] = 7; //   C3*
  AA[p].atom[26].ta[2] = 11; //   C1*

  //bk   template for building 2HO*
  AA[p].atom[27].ta[0] = 10; //   C2*
  AA[p].atom[27].ta[1] = 9; //   C3*
  AA[p].atom[27].ta[2] = 11; //   C1*

  //bk   template for building H1*
  AA[p].atom[28].ta[0] = 11; //   C1*
  AA[p].atom[28].ta[1] = 9; //   C2*
  AA[p].atom[28].ta[2] = 6; //   O4*

  //bk   template for building H2
  AA[p].atom[29].ta[0] = 13; //   C2
  AA[p].atom[29].ta[1] = 12; //   N1
  AA[p].atom[29].ta[2] = 14; //   N3

  //bk   template for building 1H6
  AA[p].atom[30].ta[0] = 18; //   N6
  AA[p].atom[30].ta[1] = 17; //   C6
  AA[p].atom[30].ta[2] = 12; //   N1

  //bk   template for building 2H6
  AA[p].atom[31].ta[0] = 18; //   N6
  AA[p].atom[31].ta[1] = 17; //   C6
  AA[p].atom[31].ta[2] = 12; //   N1

  //bk   template for building H8
  AA[p].atom[32].ta[0] = 20; //   C8
  AA[p].atom[32].ta[1] = 19; //   N7
  AA[p].atom[32].ta[2] = 21; //   N9


  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 11; //   C1*
  AA[p].chi_atoms[0] [1] = 9; //   C2*
  AA[p].chi_atoms[0] [2] = 10; //   O2*
  AA[p].chi_atoms[0] [3] = 27; //   2HO*



  //bk   template coordinates for the amino acid
  //jjh new template from Kevin Karplus

  AA[p].atom[0].icoor[0] = 63.4550; //   P
  AA[p].atom[0].icoor[1] = 24.2390; //   P
  AA[p].atom[0].icoor[2] = 11.6530; //   P
  AA[p].atom[1].icoor[0] = 64.6390; //   O1P
  AA[p].atom[1].icoor[1] = 24.3570; //   O1P
  AA[p].atom[1].icoor[2] = 10.7740; //   O1P
  AA[p].atom[2].icoor[0] = 62.1120; //   O2P
  AA[p].atom[2].icoor[1] = 24.0280; //   O2P
  AA[p].atom[2].icoor[2] = 11.0630; //   O2P
  AA[p].atom[3].icoor[0] = 63.7200; //   O5*
  AA[p].atom[3].icoor[1] = 23.2450; //   O5*
  AA[p].atom[3].icoor[2] = 12.8570; //   O5*
  AA[p].atom[4].icoor[0] = 64.8710; //   C5*
  AA[p].atom[4].icoor[1] = 23.6550; //   C5*
  AA[p].atom[4].icoor[2] = 13.5960; //   C5*
  AA[p].atom[5].icoor[0] = 64.8100; //   C4*
  AA[p].atom[5].icoor[1] = 23.0350; //   C4*
  AA[p].atom[5].icoor[2] = 14.9800; //   C4*
  AA[p].atom[6].icoor[0] = 63.4270; //   O4*
  AA[p].atom[6].icoor[1] = 22.9230; //   O4*
  AA[p].atom[6].icoor[2] = 15.3940; //   O4*
  AA[p].atom[7].icoor[0] = 65.3320; //   C3*
  AA[p].atom[7].icoor[1] = 21.5960; //   C3*
  AA[p].atom[7].icoor[2] = 15.0160; //   C3*
  AA[p].atom[8].icoor[0] = 65.9880; //   O3*
  AA[p].atom[8].icoor[1] = 21.3280; //   O3*
  AA[p].atom[8].icoor[2] = 16.2600; //   O3*
  AA[p].atom[9].icoor[0] = 64.0920; //   C2*
  AA[p].atom[9].icoor[1] = 20.7690; //   C2*
  AA[p].atom[9].icoor[2] = 14.8940; //   C2*
  AA[p].atom[10].icoor[0] = 64.2970; //   O2*
  AA[p].atom[10].icoor[1] = 19.5250; //   O2*
  AA[p].atom[10].icoor[2] = 15.5020; //   O2*
  AA[p].atom[11].icoor[0] = 63.1460; //   C1*
  AA[p].atom[11].icoor[1] = 21.5780; //   C1*
  AA[p].atom[11].icoor[2] = 15.7040; //   C1*
  AA[p].atom[12].icoor[0] = 58.5220; //   N1
  AA[p].atom[12].icoor[1] = 20.7450; //   N1
  AA[p].atom[12].icoor[2] = 17.6490; //   N1
  AA[p].atom[13].icoor[0] = 59.6770; //   C2
  AA[p].atom[13].icoor[1] = 21.1050; //   C2
  AA[p].atom[13].icoor[2] = 18.1780; //   C2
  AA[p].atom[14].icoor[0] = 60.8560; //   N3
  AA[p].atom[14].icoor[1] = 21.3330; //   N3
  AA[p].atom[14].icoor[2] = 17.6420; //   N3
  AA[p].atom[15].icoor[0] = 60.7580; //   C4
  AA[p].atom[15].icoor[1] = 21.1330; //   C4
  AA[p].atom[15].icoor[2] = 16.3550; //   C4
  AA[p].atom[16].icoor[0] = 59.6310; //   C5
  AA[p].atom[16].icoor[1] = 20.7760; //   C5
  AA[p].atom[16].icoor[2] = 15.6710; //   C5
  AA[p].atom[17].icoor[0] = 58.4460; //   C6
  AA[p].atom[17].icoor[1] = 20.5610; //   C6
  AA[p].atom[17].icoor[2] = 16.3160; //   C6
  AA[p].atom[18].icoor[0] = 57.3710; //   N6
  AA[p].atom[18].icoor[1] = 20.2000; //   N6
  AA[p].atom[18].icoor[2] = 15.5580; //   N6
  AA[p].atom[19].icoor[0] = 59.8960; //   N7
  AA[p].atom[19].icoor[1] = 20.6600; //   N7
  AA[p].atom[19].icoor[2] = 14.3430; //   N7
  AA[p].atom[20].icoor[0] = 61.1660; //   C8
  AA[p].atom[20].icoor[1] = 20.9340; //   C8
  AA[p].atom[20].icoor[2] = 14.2250; //   C8
  AA[p].atom[21].icoor[0] = 61.7590; //   N9
  AA[p].atom[21].icoor[1] = 21.2270; //   N9
  AA[p].atom[21].icoor[2] = 15.4240; //   N9
  AA[p].atom[22].icoor[0] = 64.8640; //  1H5*
  AA[p].atom[22].icoor[1] = 24.7620; //  1H5*
  AA[p].atom[22].icoor[2] = 13.6820; //  1H5*
  AA[p].atom[23].icoor[0] = 65.7960; //  2H5*
  AA[p].atom[23].icoor[1] = 23.3450; //  2H5*
  AA[p].atom[23].icoor[2] = 13.0650; //  2H5*
  AA[p].atom[24].icoor[0] = 65.3760; //   H4*
  AA[p].atom[24].icoor[1] = 23.6420; //   H4*
  AA[p].atom[24].icoor[2] = 15.7170; //   H4*
  AA[p].atom[25].icoor[0] = 66.0350; //   H3*
  AA[p].atom[25].icoor[1] = 21.4230; //   H3*
  AA[p].atom[25].icoor[2] = 14.1740; //   H3*
  AA[p].atom[26].icoor[0] = 63.7730; //  1H2*
  AA[p].atom[26].icoor[1] = 20.6530; //  1H2*
  AA[p].atom[26].icoor[2] = 13.8360; //  1H2*
  AA[p].atom[27].icoor[0] = 63.4860; //  2HO*
  AA[p].atom[27].icoor[1] = 19.0140; //  2HO*
  AA[p].atom[27].icoor[2] = 15.4420; //  2HO*
  AA[p].atom[28].icoor[0] = 63.3700; //   H1*
  AA[p].atom[28].icoor[1] = 21.4680; //   H1*
  AA[p].atom[28].icoor[2] = 16.7890; //   H1*
  AA[p].atom[29].icoor[0] = 59.6540; //   H2
  AA[p].atom[29].icoor[1] = 21.2400; //   H2
  AA[p].atom[29].icoor[2] = 19.2590; //   H2
  AA[p].atom[30].icoor[0] = 57.4310; //  1H6
  AA[p].atom[30].icoor[1] = 20.1300; //  1H6
  AA[p].atom[30].icoor[2] = 14.5620; //  1H6
  AA[p].atom[31].icoor[0] = 56.4810; //  2H6
  AA[p].atom[31].icoor[1] = 20.0110; //  2H6
  AA[p].atom[31].icoor[2] = 15.9740; //  2H6
  AA[p].atom[32].icoor[0] = 61.7180; //   H8
  AA[p].atom[32].icoor[1] = 20.9440; //   H8
  AA[p].atom[32].icoor[2] = 13.2860; //   H8

  //FIN DE rADE

  //rCYT
//  p = RCY;
  p = CYT;

//  strcpy( AA[p].aa_name3, "  rC" );
  strcpy( AA[p].aa_name3, "  C" );
//  AA[p].aa_name1 = 'X';
  AA[p].aa_name1 = 'c';

  AA[p].mass = 0;

  AA[p].is_protein = false;
  AA[p].is_RNA = true;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = true;
  AA[p].aa_is_nonpolar = false;
  AA[p].aa_is_aromatic = true;
  AA[p].aa_is_charged = true;

  AA[p].natoms = 31;
  AA[p].nheavyatoms = 20;
  AA[p].nchi = 1;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );
  AA[p].atom = ( t_aa_atom_p * ) malloc( AA[p].natoms * sizeof( t_aa_atom_p ) );


  strcpy( AA[p].atom[0].atom_name, " P  " );
//  strcpy( AA[p].atom[1].atom_name, " O1P" );
//  strcpy( AA[p].atom[2].atom_name, " O2P" );

  strcpy( AA[p].atom[1].atom_name, " OP1" );
  strcpy( AA[p].atom[2].atom_name, " OP2" );

//  strcpy( AA[p].atom[3].atom_name, " O5*" );
//  strcpy( AA[p].atom[4].atom_name, " C5*" );
//  strcpy( AA[p].atom[5].atom_name, " C4*" );
//  strcpy( AA[p].atom[6].atom_name, " O4*" );
//  strcpy( AA[p].atom[7].atom_name, " C3*" );
//  strcpy( AA[p].atom[8].atom_name, " O3*" );
//  strcpy( AA[p].atom[9].atom_name, " C2*" );
//  strcpy( AA[p].atom[10].atom_name, " O2*" );
//  strcpy( AA[p].atom[11].atom_name, " C1*" );

  strcpy( AA[p].atom[3].atom_name, " O5'" );
  strcpy( AA[p].atom[4].atom_name, " C5'" );
  strcpy( AA[p].atom[5].atom_name, " C4'" );
  strcpy( AA[p].atom[6].atom_name, " O4'" );
  strcpy( AA[p].atom[7].atom_name, " C3'" );
  strcpy( AA[p].atom[8].atom_name, " O3'" );
  strcpy( AA[p].atom[9].atom_name, " C2'" );
  strcpy( AA[p].atom[10].atom_name, " O2'" );
  strcpy( AA[p].atom[11].atom_name, " C1'" );

  strcpy( AA[p].atom[12].atom_name, " N1 " );
  strcpy( AA[p].atom[13].atom_name, " C2 " );
  strcpy( AA[p].atom[14].atom_name, " O2 " );
  strcpy( AA[p].atom[15].atom_name, " N3 " );
  strcpy( AA[p].atom[16].atom_name, " C4 " );
  strcpy( AA[p].atom[17].atom_name, " N4 " );
  strcpy( AA[p].atom[18].atom_name, " C5 " );
  strcpy( AA[p].atom[19].atom_name, " C6 " );

//  strcpy( AA[p].atom[20].atom_name, "1H5*" );
//  strcpy( AA[p].atom[21].atom_name, "2H5*" );
//  strcpy( AA[p].atom[22].atom_name, " H4*" );
//  strcpy( AA[p].atom[23].atom_name, " H3*" );
//  strcpy( AA[p].atom[24].atom_name, "1H2*" );
//  strcpy( AA[p].atom[25].atom_name, "2HO*" );
//  strcpy( AA[p].atom[26].atom_name, " H1*" );

  strcpy( AA[p].atom[20].atom_name, "1H5'" );
  strcpy( AA[p].atom[21].atom_name, "2H5'" );
  strcpy( AA[p].atom[22].atom_name, " H4'" );
  strcpy( AA[p].atom[23].atom_name, " H3'" );
  strcpy( AA[p].atom[24].atom_name, "1H2'" );
  strcpy( AA[p].atom[25].atom_name, "2HO'" );
  strcpy( AA[p].atom[26].atom_name, " H1'" );

  strcpy( AA[p].atom[27].atom_name, "1H4 " );
  strcpy( AA[p].atom[28].atom_name, "2H4 " );
  strcpy( AA[p].atom[29].atom_name, " H5 " );
  strcpy( AA[p].atom[30].atom_name, " H6 " );

  AA[p].atom[0].charge = 1.500; // P
  AA[p].atom[1].charge = -.780; // O1P
  AA[p].atom[2].charge = -.780; // O2P
  AA[p].atom[3].charge = -.570; // O5*
  AA[p].atom[4].charge = -.080; // C5*
  AA[p].atom[5].charge =.160; // C4*
  AA[p].atom[6].charge = -.500; // O4*
  AA[p].atom[7].charge =.010; // C3*
  AA[p].atom[8].charge = -.570; // O3*
  AA[p].atom[9].charge =.140; // C2*
  AA[p].atom[10].charge = -.660; // O2*
  AA[p].atom[11].charge =.160; // C1*
  AA[p].atom[12].charge = -.130; // N1
  AA[p].atom[13].charge =.520; // C2
  AA[p].atom[14].charge = -.490; // O2
  AA[p].atom[15].charge = -.660; // N3
  AA[p].atom[16].charge =.650; // C4
  AA[p].atom[17].charge = -.750; // N4
  AA[p].atom[18].charge = -.130; // C5
  AA[p].atom[19].charge =.050; // C6
  AA[p].atom[20].charge =.090; // 1H5*
  AA[p].atom[21].charge =.090; // 2H5*
  AA[p].atom[22].charge =.090; // H4*
  AA[p].atom[23].charge =.090; // H3*
  AA[p].atom[24].charge =.090; // 1H2*
  AA[p].atom[25].charge =.430; // 2HO*
  AA[p].atom[26].charge =.090; // H1*
  AA[p].atom[27].charge =.370; // 1H4
  AA[p].atom[28].charge =.330; // 2H4
  AA[p].atom[29].charge =.070; // H5
  AA[p].atom[30].charge =.170; // H6


  switch(opt)
  {
  case Sybil:
	  AA[p].atom[0].fullatom_type = 12; // Phos   P     P
	  AA[p].atom[1].fullatom_type = 10; // 69; // OOC    O1P   ON3
	  AA[p].atom[2].fullatom_type = 10; // 69; // OOC    O2P   ON3
	  AA[p].atom[3].fullatom_type = 9; // 68; // ONH2   O5*   ON2
	  AA[p].atom[4].fullatom_type = 2; // 61; // CH2    C5*   CN8 (CN8B)
	  AA[p].atom[5].fullatom_type = 2; // 60; // CH1    C4*   CN7
	  AA[p].atom[6].fullatom_type = 9; // 72; // OH     O4*   ON6 (ON6B)
	  AA[p].atom[7].fullatom_type = 2; // 60; // CH1    C3*   CN7
	  AA[p].atom[8].fullatom_type = 9; // 68; // ONH2   O3*   ON2
	  AA[p].atom[9].fullatom_type = 2; // 60; // CH2    C2*   CN7 (CN7B)
	  AA[p].atom[10].fullatom_type = 9; // 71; // OH     O2*   ON5
	  AA[p].atom[11].fullatom_type = 2; // 60; // CH1    C1*   CN7 (CN7B)
	  AA[p].atom[12].fullatom_type = 5; // 64; // Nhis   N1    NN2
	  AA[p].atom[13].fullatom_type = 1; // 55; // CObb   C2    CN1
	  AA[p].atom[14].fullatom_type = 8; // 67; // OCbb   O2    ON1 (ON1C)
	  AA[p].atom[15].fullatom_type = 5; // 65; // Nhis   N3    NN3
	  AA[p].atom[16].fullatom_type = 3; // 56; // aroC   C4    CN2
	  AA[p].atom[17].fullatom_type = 5; // 63; // NH2O   N4    NN1
	  AA[p].atom[18].fullatom_type = 3; // 57; // aroC   C5    CN3
	  AA[p].atom[19].fullatom_type = 3; // 57; // aroC   C6    CN3
	  AA[p].atom[20].fullatom_type = 14; // 74; // Hapo  1H5*   HNP
	  AA[p].atom[21].fullatom_type = 14; // 74; // Hapo  2H5*   HNP
	  AA[p].atom[22].fullatom_type = 14; // 74; // Hapo   H4*   HNP
	  AA[p].atom[23].fullatom_type = 14; // 74; // Hapo   H3*   HNP
	  AA[p].atom[24].fullatom_type = 14; // 74; // Hapo  1H2*   HNP
	  AA[p].atom[25].fullatom_type = 13; // 73; // Hpol  2HO*   HP
	  AA[p].atom[26].fullatom_type = 14; // 74; // Hapo   H1*   HNP
	  AA[p].atom[27].fullatom_type = 13; // 73; // Hpol  1H4    HP
	  AA[p].atom[28].fullatom_type = 13; // 73; // Hpol  2H4    HP
	  AA[p].atom[29].fullatom_type = 15; // 74; // Haro   H5    HNP
	  AA[p].atom[30].fullatom_type = 15; // 74; // Haro   H6    HNP
	  break;

  default:
	  AA[p].atom[0].fullatom_type = 21; // Phos   P     P
	  AA[p].atom[1].fullatom_type = 15; // 69; // OOC    O1P   ON3
	  AA[p].atom[2].fullatom_type = 15; // 69; // OOC    O2P   ON3
	  AA[p].atom[3].fullatom_type = 14; // 68; // ONH2   O5*   ON2
	  AA[p].atom[4].fullatom_type = 4; // 61; // CH2    C5*   CN8 (CN8B)
	  AA[p].atom[5].fullatom_type = 3; // 60; // CH1    C4*   CN7
	  AA[p].atom[6].fullatom_type = 13; // 72; // OH     O4*   ON6 (ON6B)
	  AA[p].atom[7].fullatom_type = 3; // 60; // CH1    C3*   CN7
	  AA[p].atom[8].fullatom_type = 14; // 68; // ONH2   O3*   ON2
	  AA[p].atom[9].fullatom_type = 4; // 60; // CH2    C2*   CN7 (CN7B)
	  AA[p].atom[10].fullatom_type = 13; // 71; // OH     O2*   ON5
	  AA[p].atom[11].fullatom_type = 3; // 60; // CH1    C1*   CN7 (CN7B)
	  AA[p].atom[12].fullatom_type = 8; // 64; // Nhis   N1    NN2
	  AA[p].atom[13].fullatom_type = 19; // 55; // CObb   C2    CN1
	  AA[p].atom[14].fullatom_type = 20; // 67; // OCbb   O2    ON1 (ON1C)
	  AA[p].atom[15].fullatom_type = 8; // 65; // Nhis   N3    NN3
	  AA[p].atom[16].fullatom_type = 6; // 56; // aroC   C4    CN2
	  AA[p].atom[17].fullatom_type = 9; // 63; // NH2O   N4    NN1
	  AA[p].atom[18].fullatom_type = 6; // 57; // aroC   C5    CN3
	  AA[p].atom[19].fullatom_type = 6; // 57; // aroC   C6    CN3
	  AA[p].atom[20].fullatom_type = 23; // 74; // Hapo  1H5*   HNP
	  AA[p].atom[21].fullatom_type = 23; // 74; // Hapo  2H5*   HNP
	  AA[p].atom[22].fullatom_type = 23; // 74; // Hapo   H4*   HNP
	  AA[p].atom[23].fullatom_type = 23; // 74; // Hapo   H3*   HNP
	  AA[p].atom[24].fullatom_type = 23; // 74; // Hapo  1H2*   HNP
	  AA[p].atom[25].fullatom_type = 22; // 73; // Hpol  2HO*   HP
	  AA[p].atom[26].fullatom_type = 23; // 74; // Hapo   H1*   HNP
	  AA[p].atom[27].fullatom_type = 22; // 73; // Hpol  1H4    HP
	  AA[p].atom[28].fullatom_type = 22; // 73; // Hpol  2H4    HP
	  AA[p].atom[29].fullatom_type = 24; // 74; // Haro   H5    HNP
	  AA[p].atom[30].fullatom_type = 24; // 74; // Haro   H6    HNP
	  break;
  }
  //jjh intra residue bonding

  AA[p].atom[0].nbonded_neighbors = 3; // P
  AA[p].atom[0].bonded_neighbor[0] = 1; // P--O1P
  AA[p].atom[0].bonded_neighbor[1] = 2; // P--O2P
  AA[p].atom[0].bonded_neighbor[2] = 3; // P--O5*
  AA[p].atom[1].nbonded_neighbors = 1; // O1P
  AA[p].atom[1].bonded_neighbor[0] = 0; // O1P--P
  AA[p].atom[2].nbonded_neighbors = 1; // O2P
  AA[p].atom[2].bonded_neighbor[0] = 0; // O2P--P
  AA[p].atom[3].nbonded_neighbors = 2; // O5*
  AA[p].atom[3].bonded_neighbor[0] = 0; // O5*--P
  AA[p].atom[3].bonded_neighbor[1] = 4; // O5*--C5*
  AA[p].atom[4].nbonded_neighbors = 4; // C5*
  AA[p].atom[4].bonded_neighbor[0] = 3; // C5*--O5*
  AA[p].atom[4].bonded_neighbor[1] = 5; // C5*--C4*
  AA[p].atom[4].bonded_neighbor[2] = 20; // C5*--1H5*
  AA[p].atom[4].bonded_neighbor[3] = 21; // C5*--2H5*
  AA[p].atom[5].nbonded_neighbors = 4; // C4*
  AA[p].atom[5].bonded_neighbor[0] = 4; // C4*--C5*
  AA[p].atom[5].bonded_neighbor[1] = 6; // C4*--O4*
  AA[p].atom[5].bonded_neighbor[2] = 7; // C4*--C3*
  AA[p].atom[5].bonded_neighbor[3] = 22; // C4*--H4*
  AA[p].atom[6].nbonded_neighbors = 2; // O4*
  AA[p].atom[6].bonded_neighbor[0] = 5; // O4*--C4*
  AA[p].atom[6].bonded_neighbor[1] = 11; // O4*--C1*
  AA[p].atom[7].nbonded_neighbors = 4; // C3*
  AA[p].atom[7].bonded_neighbor[0] = 5; // C3*--C4*
  AA[p].atom[7].bonded_neighbor[1] = 8; // C3*--O3*
  AA[p].atom[7].bonded_neighbor[2] = 9; // C3*--C2*
  AA[p].atom[7].bonded_neighbor[3] = 23; // C3*--H3*
  AA[p].atom[8].nbonded_neighbors = 1; // O3*
  AA[p].atom[8].bonded_neighbor[0] = 7; // O3*--C3*
  AA[p].atom[9].nbonded_neighbors = 4; // C2*
  AA[p].atom[9].bonded_neighbor[0] = 7; // C2*--C3*
  AA[p].atom[9].bonded_neighbor[1] = 10; // C2*--02*
  AA[p].atom[9].bonded_neighbor[2] = 11; // C2*--C1*
  AA[p].atom[9].bonded_neighbor[3] = 24; // C2*--1H2*
  AA[p].atom[10].nbonded_neighbors = 4; // 02*
  AA[p].atom[10].bonded_neighbor[0] = 9; // O2*--C2*
  AA[p].atom[10].bonded_neighbor[1] = 25; // 02*--2HO*
  AA[p].atom[11].nbonded_neighbors = 4; // C1*
  AA[p].atom[11].bonded_neighbor[0] = 6; // C1*--O4*
  AA[p].atom[11].bonded_neighbor[1] = 9; // C1*--C2*
  AA[p].atom[11].bonded_neighbor[2] = 12; // C1*--N1
  AA[p].atom[11].bonded_neighbor[3] = 26; // C1*--H1*
  AA[p].atom[12].nbonded_neighbors = 3; // N1
  AA[p].atom[12].bonded_neighbor[0] = 11; // N1--C1*
  AA[p].atom[12].bonded_neighbor[1] = 13; // N1--C2
  AA[p].atom[12].bonded_neighbor[2] = 19; // N1--C6
  AA[p].atom[13].nbonded_neighbors = 3; // C2
  AA[p].atom[13].bonded_neighbor[0] = 12; // C2--N1
  AA[p].atom[13].bonded_neighbor[1] = 14; // C2--O2
  AA[p].atom[13].bonded_neighbor[2] = 15; // C2--N3
  AA[p].atom[14].nbonded_neighbors = 1; // O2
  AA[p].atom[14].bonded_neighbor[0] = 13; // O2--C2
  AA[p].atom[15].nbonded_neighbors = 2; // N3
  AA[p].atom[15].bonded_neighbor[0] = 13; // N3--C2
  AA[p].atom[15].bonded_neighbor[1] = 16; // N3--C4
  AA[p].atom[16].nbonded_neighbors = 3; // C4
  AA[p].atom[16].bonded_neighbor[0] = 15; // C4--N3
  AA[p].atom[16].bonded_neighbor[1] = 17; // C4--N4
  AA[p].atom[16].bonded_neighbor[2] = 18; // C4--C5
  AA[p].atom[17].nbonded_neighbors = 3; // N4
  AA[p].atom[17].bonded_neighbor[0] = 16; // N4--C4
  AA[p].atom[17].bonded_neighbor[1] = 27; // N4--1H4
  AA[p].atom[17].bonded_neighbor[2] = 28; // N4--2H4
  AA[p].atom[18].nbonded_neighbors = 3; // C5
  AA[p].atom[18].bonded_neighbor[0] = 16; // C5--C4
  AA[p].atom[18].bonded_neighbor[1] = 19; // C5--C6
  AA[p].atom[18].bonded_neighbor[2] = 29; // C5--H5
  AA[p].atom[19].nbonded_neighbors = 3; // C6
  AA[p].atom[19].bonded_neighbor[0] = 12; // C6--N1
  AA[p].atom[19].bonded_neighbor[1] = 18; // C6--C5
  AA[p].atom[19].bonded_neighbor[2] = 30; // C6--H6
  AA[p].atom[20].nbonded_neighbors = 1; //1H5*
  AA[p].atom[20].bonded_neighbor[0] = 4; //1H5*--C5*
  AA[p].atom[21].nbonded_neighbors = 1; //2H5*
  AA[p].atom[21].bonded_neighbor[0] = 4; //2H5*--C5*
  AA[p].atom[22].nbonded_neighbors = 1; // H4*
  AA[p].atom[22].bonded_neighbor[0] = 5; // H4*--C4*
  AA[p].atom[23].nbonded_neighbors = 1; // H3*
  AA[p].atom[23].bonded_neighbor[0] = 7; // H3*--C3*
  AA[p].atom[24].nbonded_neighbors = 1; //1H2*
  AA[p].atom[24].bonded_neighbor[0] = 9; // 1H2*--C2*
  AA[p].atom[25].nbonded_neighbors = 1; //2HO*
  AA[p].atom[25].bonded_neighbor[0] = 10; // 2HO*--C2*
  AA[p].atom[26].nbonded_neighbors = 1; // H1*
  AA[p].atom[26].bonded_neighbor[0] = 11; // H1*--C1*
  AA[p].atom[27].nbonded_neighbors = 1; //1H4
  AA[p].atom[27].bonded_neighbor[0] = 17; //1H4--N4
  AA[p].atom[28].nbonded_neighbors = 1; //2H4
  AA[p].atom[28].bonded_neighbor[0] = 17; //2H4--N4
  AA[p].atom[29].nbonded_neighbors = 1; // H5
  AA[p].atom[29].bonded_neighbor[0] = 18; // H5--C5
  AA[p].atom[30].nbonded_neighbors = 1; // H6
  AA[p].atom[30].bonded_neighbor[0] = 19; // H6--C6

  for ( k = 0; k < 3; k++ )
    AA[p].atom[k].hastemplate = true;

  for ( k = 3; k < 20; k++ )
    AA[p].atom[k].hastemplate = false;

  AA[p].atom[10].hastemplate = false;

  for ( k = 30; k < 31; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template atoms used for placing atoms rCYT

  //bk   template for building P
  AA[p].atom[0].ta[0] = 3; //   O5*
  AA[p].atom[0].ta[1] = 4; //   C5*
  AA[p].atom[0].ta[2] = 5; //   C4*

  //bk   template for building O1P
  AA[p].atom[1].ta[0] = 0; //   P
  AA[p].atom[1].ta[1] = 4; //   C5*
  AA[p].atom[1].ta[2] = 3; //   O5*

  //bk   template for building  O2P
  AA[p].atom[2].ta[0] = 0; //   P
  AA[p].atom[2].ta[1] = 1; //   O1P
  AA[p].atom[2].ta[2] = 3; //   O5*

  //bk   template for building  O2*
  AA[p].atom[2].ta[0] = 9; //   C2*
  AA[p].atom[2].ta[1] = 7; //   C3*
  AA[p].atom[2].ta[2] = 11; //   C1*


  //bk   template for building 1H5*
  AA[p].atom[10].ta[0] = 4; //   C5*
  AA[p].atom[10].ta[1] = 3; //   O5*
  AA[p].atom[10].ta[2] = 5; //   C4*

  //bk   template for building 2H5*
  AA[p].atom[20].ta[0] = 4; //   C5*
  AA[p].atom[20].ta[1] = 3; //   O5*
  AA[p].atom[20].ta[2] = 5; //   C4*

  //bk   template for building H4*
  AA[p].atom[21].ta[0] = 5; //   C4*
  AA[p].atom[21].ta[1] = 6; //   O4*
  AA[p].atom[21].ta[2] = 4; //   C5*

  //bk   template for building H3*
  AA[p].atom[22].ta[0] = 7; //   C3*
  AA[p].atom[22].ta[1] = 9; //   C2*
  AA[p].atom[22].ta[2] = 5; //   C4*

  //bk   template for building 1H2*
  AA[p].atom[23].ta[0] = 9; //   C2*
  AA[p].atom[23].ta[1] = 7; //   C3*
  AA[p].atom[23].ta[2] = 11; //   C1*

  //bk   template for building 2HO*
  AA[p].atom[24].ta[0] = 10; //   C2*
  AA[p].atom[24].ta[1] = 9; //   C3*
  AA[p].atom[24].ta[2] = 11; //   C1*

  //bk   template for building H1*
  AA[p].atom[25].ta[0] = 11; //   C1*
  AA[p].atom[25].ta[1] = 9; //   C2*
  AA[p].atom[25].ta[2] = 6; //   O4*

  //bk   template for building 1H4
  AA[p].atom[26].ta[0] = 17; //   N4
  AA[p].atom[26].ta[1] = 16; //   C4
  AA[p].atom[26].ta[2] = 15; //   N3

  //bk   template for building 2H4
  AA[p].atom[27].ta[0] = 17; //   N4
  AA[p].atom[27].ta[1] = 16; //   C4
  AA[p].atom[27].ta[2] = 15; //   N3

  //bk   template for building H5
  AA[p].atom[28].ta[0] = 18; //   C5
  AA[p].atom[28].ta[1] = 19; //   C6
  AA[p].atom[28].ta[2] = 12; //   N1

  //bk   template for building H6
  AA[p].atom[29].ta[0] = 19; //   C6
  AA[p].atom[29].ta[1] = 18; //   C5
  AA[p].atom[29].ta[2] = 16; //   C4


  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 11; //   C1*
  AA[p].chi_atoms[0] [1] = 9; //   C2*
  AA[p].chi_atoms[0] [2] = 10; //   O2*
  AA[p].chi_atoms[0] [3] = 25; //   2HO*


  //bk   template coordinates for the amino acid
  //jjh new template from Kevin Karplus

  AA[p].atom[0].icoor[0] = 55.3470; //   P
  AA[p].atom[0].icoor[1] = 21.8570; //   P
  AA[p].atom[0].icoor[2] = 26.7560; //   P
  AA[p].atom[1].icoor[0] = 55.0590; //   O1P
  AA[p].atom[1].icoor[1] = 21.4890; //   O1P
  AA[p].atom[1].icoor[2] = 28.1600; //   O1P
  AA[p].atom[2].icoor[0] = 54.2580; //   O2P
  AA[p].atom[2].icoor[1] = 21.8630; //   O2P
  AA[p].atom[2].icoor[2] = 25.7540; //   O2P
  AA[p].atom[3].icoor[0] = 56.2350; //   O5*
  AA[p].atom[3].icoor[1] = 23.1670; //   O5*
  AA[p].atom[3].icoor[2] = 26.6600; //   O5*
  AA[p].atom[4].icoor[0] = 57.5350; //   C5*
  AA[p].atom[4].icoor[1] = 22.9130; //   C5*
  AA[p].atom[4].icoor[2] = 27.1930; //   C5*
  AA[p].atom[5].icoor[0] = 58.5920; //   C4*
  AA[p].atom[5].icoor[1] = 23.2500; //   C4*
  AA[p].atom[5].icoor[2] = 26.1530; //   C4*
  AA[p].atom[6].icoor[0] = 58.1460; //   O4*
  AA[p].atom[6].icoor[1] = 22.8270; //   O4*
  AA[p].atom[6].icoor[2] = 24.8440; //   O4*
  AA[p].atom[7].icoor[0] = 58.8920; //   C3*
  AA[p].atom[7].icoor[1] = 24.7420; //   C3*
  AA[p].atom[7].icoor[2] = 25.9920; //   C3*
  AA[p].atom[8].icoor[0] = 60.3010; //   O3*
  AA[p].atom[8].icoor[1] = 24.9480; //   O3*
  AA[p].atom[8].icoor[2] = 25.8380; //   O3*
  AA[p].atom[9].icoor[0] = 58.1830; //   C2*
  AA[p].atom[9].icoor[1] = 25.1250; //   C2*
  AA[p].atom[9].icoor[2] = 24.7290; //   C2*
  AA[p].atom[10].icoor[0] = 58.8740; //   O2*
  AA[p].atom[10].icoor[1] = 26.1750; //   O2*
  AA[p].atom[10].icoor[2] = 24.1130; //   O2*
  AA[p].atom[11].icoor[0] = 58.3430; //   C1*
  AA[p].atom[11].icoor[1] = 23.8790; //   C1*
  AA[p].atom[11].icoor[2] = 23.9330; //   C1*
  AA[p].atom[12].icoor[0] = 57.3960; //   N1
  AA[p].atom[12].icoor[1] = 23.8290; //   N1
  AA[p].atom[12].icoor[2] = 22.8140; //   N1
  AA[p].atom[13].icoor[0] = 57.8040; //   C2
  AA[p].atom[13].icoor[1] = 24.0050; //   C2
  AA[p].atom[13].icoor[2] = 21.5160; //   C2
  AA[p].atom[14].icoor[0] = 59.0080; //   O2
  AA[p].atom[14].icoor[1] = 24.2080; //   O2
  AA[p].atom[14].icoor[2] = 21.2750; //   O2
  AA[p].atom[15].icoor[0] = 56.8670; //   N3
  AA[p].atom[15].icoor[1] = 23.9460; //   N3
  AA[p].atom[15].icoor[2] = 20.5610; //   N3
  AA[p].atom[16].icoor[0] = 55.5690; //   C4
  AA[p].atom[16].icoor[1] = 23.7540; //   C4
  AA[p].atom[16].icoor[2] = 20.8620; //   C4
  AA[p].atom[17].icoor[0] = 54.6040; //   N4
  AA[p].atom[17].icoor[1] = 23.6870; //   N4
  AA[p].atom[17].icoor[2] = 19.8980; //   N4
  AA[p].atom[18].icoor[0] = 55.1640; //   C5
  AA[p].atom[18].icoor[1] = 23.6160; //   C5
  AA[p].atom[18].icoor[2] = 22.2010; //   C5
  AA[p].atom[19].icoor[0] = 56.1040; //   C6
  AA[p].atom[19].icoor[1] = 23.6650; //   C6
  AA[p].atom[19].icoor[2] = 23.1290; //   C6
  AA[p].atom[20].icoor[0] = 57.6090; //  1H5*
  AA[p].atom[20].icoor[1] = 21.8340; //  1H5*
  AA[p].atom[20].icoor[2] = 27.4530; //  1H5*
  AA[p].atom[21].icoor[0] = 57.6880; //  2H5*
  AA[p].atom[21].icoor[1] = 23.5070; //  2H5*
  AA[p].atom[21].icoor[2] = 28.1200; //  2H5*
  AA[p].atom[22].icoor[0] = 59.5480; //   H4*
  AA[p].atom[22].icoor[1] = 22.7350; //   H4*
  AA[p].atom[22].icoor[2] = 26.3910; //   H4*
  AA[p].atom[23].icoor[0] = 58.5210; //   H3*
  AA[p].atom[23].icoor[1] = 25.3170; //   H3*
  AA[p].atom[23].icoor[2] = 26.8680; //   H3*
  AA[p].atom[24].icoor[0] = 57.1170; //  1H2*
  AA[p].atom[24].icoor[1] = 25.3880; //  1H2*
  AA[p].atom[24].icoor[2] = 24.9050; //  1H2*
  AA[p].atom[25].icoor[0] = 58.4180; //  2HO*
  AA[p].atom[25].icoor[1] = 26.4160; //  2HO*
  AA[p].atom[25].icoor[2] = 23.3020; //  2HO*
  AA[p].atom[26].icoor[0] = 59.3960; //   H1*
  AA[p].atom[26].icoor[1] = 23.7810; //   H1*
  AA[p].atom[26].icoor[2] = 23.5900; //   H1*
  AA[p].atom[27].icoor[0] = 54.8240; //  1H4
  AA[p].atom[27].icoor[1] = 23.7780; //  1H4
  AA[p].atom[27].icoor[2] = 18.9260; //  1H4
  AA[p].atom[28].icoor[0] = 53.6410; //  2H4
  AA[p].atom[28].icoor[1] = 23.5620; //  2H4
  AA[p].atom[28].icoor[2] = 20.1350; //  2H4
  AA[p].atom[29].icoor[0] = 54.1130; //   H5
  AA[p].atom[29].icoor[1] = 23.4720; //   H5
  AA[p].atom[29].icoor[2] = 22.4470; //   H5
  AA[p].atom[29].icoor[0] = 55.8750; //   H6
  AA[p].atom[29].icoor[1] = 23.5650; //   H6
  AA[p].atom[29].icoor[2] = 24.1890; //   H6

  //FIN DE rCYT


  //rURA
  p = URA;

//  strcpy( AA[p].aa_name3, "  rU" );
  strcpy( AA[p].aa_name3, "  U" );
//  AA[p].aa_name1 = 'X';
  AA[p].aa_name1 = 'u';
  AA[p].mass = 0;

  AA[p].is_protein = false;
  AA[p].is_RNA = true;
  AA[p].is_DNA = false;

  AA[p].aa_is_polar = true;
  AA[p].aa_is_nonpolar = false;
  AA[p].aa_is_aromatic = true;
  AA[p].aa_is_charged = true;

  AA[p].natoms = 30;
  AA[p].nheavyatoms = 20;
  AA[p].nchi = 1;

  AA[p].chi_atoms = ( entero4 * ) malloc( AA[p].nchi * sizeof( entero4 ) );
  AA[p].atom = ( t_aa_atom_p * ) malloc( AA[p].natoms * sizeof( t_aa_atom_p ) );


  strcpy( AA[p].atom[0].atom_name, " P  " );
//  strcpy( AA[p].atom[1].atom_name, " O1P" );
//  strcpy( AA[p].atom[2].atom_name, " O2P" );

  strcpy( AA[p].atom[1].atom_name, " OP1" );
  strcpy( AA[p].atom[2].atom_name, " OP2" );

//  strcpy( AA[p].atom[3].atom_name, " O5*" );
//  strcpy( AA[p].atom[4].atom_name, " C5*" );
//  strcpy( AA[p].atom[5].atom_name, " C4*" );
//  strcpy( AA[p].atom[6].atom_name, " O4*" );
//  strcpy( AA[p].atom[7].atom_name, " C3*" );
//  strcpy( AA[p].atom[8].atom_name, " O3*" );
//  strcpy( AA[p].atom[9].atom_name, " C2*" );
//  strcpy( AA[p].atom[10].atom_name, " O2*" );
//  strcpy( AA[p].atom[11].atom_name, " C1*" );

  strcpy( AA[p].atom[3].atom_name, " O5'" );
  strcpy( AA[p].atom[4].atom_name, " C5'" );
  strcpy( AA[p].atom[5].atom_name, " C4'" );
  strcpy( AA[p].atom[6].atom_name, " O4'" );
  strcpy( AA[p].atom[7].atom_name, " C3'" );
  strcpy( AA[p].atom[8].atom_name, " O3'" );
  strcpy( AA[p].atom[9].atom_name, " C2'" );
  strcpy( AA[p].atom[10].atom_name, " O2'" );
  strcpy( AA[p].atom[11].atom_name, " C1'" );

  strcpy( AA[p].atom[12].atom_name, " N1 " );
  strcpy( AA[p].atom[13].atom_name, " C2 " );
  strcpy( AA[p].atom[14].atom_name, " O2 " );
  strcpy( AA[p].atom[15].atom_name, " N3 " );
  strcpy( AA[p].atom[16].atom_name, " C4 " );
  strcpy( AA[p].atom[17].atom_name, " O4 " );
  strcpy( AA[p].atom[18].atom_name, " C5 " );
  strcpy( AA[p].atom[19].atom_name, " C6 " );

//  strcpy( AA[p].atom[20].atom_name, "1H5*" );
//  strcpy( AA[p].atom[21].atom_name, "2H5*" );
//  strcpy( AA[p].atom[22].atom_name, " H4*" );
//  strcpy( AA[p].atom[23].atom_name, " H3*" );
//  strcpy( AA[p].atom[24].atom_name, "1H2*" );
//  strcpy( AA[p].atom[25].atom_name, "2HO*" );
//  strcpy( AA[p].atom[26].atom_name, " H1*" );

  strcpy( AA[p].atom[20].atom_name, "1H5'" );
  strcpy( AA[p].atom[21].atom_name, "2H5'" );
  strcpy( AA[p].atom[22].atom_name, " H4'" );
  strcpy( AA[p].atom[23].atom_name, " H3'" );
  strcpy( AA[p].atom[24].atom_name, "1H2'" );
  strcpy( AA[p].atom[25].atom_name, "2HO'" );
  strcpy( AA[p].atom[26].atom_name, " H1'" );

  strcpy( AA[p].atom[27].atom_name, " H3 " );
  strcpy( AA[p].atom[28].atom_name, " H5 " );
  strcpy( AA[p].atom[29].atom_name, " H6 " );

  AA[p].atom[0].charge = 1.500; // P
  AA[p].atom[1].charge = -.780; // O1P
  AA[p].atom[2].charge = -.780; // O2P
  AA[p].atom[3].charge = -.570; // O5*
  AA[p].atom[4].charge = -.080; // C5*
  AA[p].atom[5].charge =.160; // C4*
  AA[p].atom[6].charge = -.500; // O4*
  AA[p].atom[7].charge =.010; // C3*
  AA[p].atom[8].charge = -.570; // O3*
  AA[p].atom[9].charge =.140; // C2*
  AA[p].atom[10].charge = -.660; // O2*
  AA[p].atom[11].charge =.160; // C1*
  AA[p].atom[12].charge = -.340; // N1
  AA[p].atom[13].charge =.550; // C2
  AA[p].atom[14].charge = -.450; // O2
  AA[p].atom[15].charge = -.460; // N3
  AA[p].atom[16].charge =.530; // C4
  AA[p].atom[17].charge = -.480; // O4
  AA[p].atom[18].charge = -.150; // C5
  AA[p].atom[19].charge =.200; // C6
  AA[p].atom[20].charge =.090; // 1H5
  AA[p].atom[21].charge =.090; // 2H5
  AA[p].atom[22].charge =.090; // H4*
  AA[p].atom[23].charge =.090; // H3*
  AA[p].atom[24].charge =.090; // 1H2
  AA[p].atom[25].charge =.430; // 2HO
  AA[p].atom[26].charge =.090; // H1*
  AA[p].atom[27].charge =.360; // H3
  AA[p].atom[28].charge =.100; // H5
  AA[p].atom[29].charge =.140; // H6

  switch(opt)
  {
  case Sybil:
	  AA[p].atom[0].fullatom_type = 12; // Phos   P     P
	  AA[p].atom[1].fullatom_type = 10; // 69; // OOC    O1P   ON3
	  AA[p].atom[2].fullatom_type = 10; // 69; // OOC    O2P   ON3
	  AA[p].atom[3].fullatom_type = 9; // 68; // ONH2   O5*   ON2
	  AA[p].atom[4].fullatom_type = 2; // 61; // CH2    C5*   CN8 (CN8B)
	  AA[p].atom[5].fullatom_type = 2; // 60; // CH1    C4*   CN7
	  AA[p].atom[6].fullatom_type = 9; // 72; // OH     O4*   ON6 (ON6B)
	  AA[p].atom[7].fullatom_type = 2; // 60; // CH1    C3*   CN7
	  AA[p].atom[8].fullatom_type = 9; // 68; // ONH2   O3*   ON2
	  AA[p].atom[9].fullatom_type = 2; // 60; // CH2    C2*   CN7 (CN7B)
	  AA[p].atom[10].fullatom_type = 9; // 67; // OH     O2*   ON1
	  AA[p].atom[11].fullatom_type = 2; // 60; // CH1    C1*   CN7 (CN7B)
	  AA[p].atom[12].fullatom_type = 5; // 64; // Nhis   N1    NN2 (NN2B)
	  AA[p].atom[13].fullatom_type = 1; // 55; // CObb   C2    CN1 (CN1T)
	  AA[p].atom[14].fullatom_type = 8; // 67; // OCbb   O2    ON1
	  AA[p].atom[15].fullatom_type = 5; // 64; // Ntrp   N3    NN2 (NN2U)
	  AA[p].atom[16].fullatom_type = 1; // 55; // CObb   C4    CN1
	  AA[p].atom[17].fullatom_type = 8; // 67; // OCbb   O4    ON1
	  AA[p].atom[18].fullatom_type = 3; // 57; // aroC   C5    CN3
	  AA[p].atom[19].fullatom_type = 3; // 57; // aroC   C6    CN3
	  AA[p].atom[20].fullatom_type = 14; // 74; // Hapo  1H5*   HNP
	  AA[p].atom[21].fullatom_type = 14; // 74; // Hapo  2H5*   HNP
	  AA[p].atom[22].fullatom_type = 14; // 74; // Hapo   H4*   HNP
	  AA[p].atom[23].fullatom_type = 14; // 74; // Hapo   H3*   HNP
	  AA[p].atom[24].fullatom_type = 14; // 74; // Hapo  1H2*   HNP
	  AA[p].atom[25].fullatom_type = 13; // 73; // Hpol  2HO*   HP
	  AA[p].atom[26].fullatom_type = 14; // 74; // Hapo   H1*   HNP
	  AA[p].atom[27].fullatom_type = 13; // 73; // Hpol   H3    HP
	  AA[p].atom[28].fullatom_type = 15; // 74; // Haro   H5    HNP
	  AA[p].atom[29].fullatom_type = 15; // 74; // Haro   H6    HNP
	  break;

  default:
	  AA[p].atom[0].fullatom_type = 21; // Phos   P     P
	  AA[p].atom[1].fullatom_type = 15; // 69; // OOC    O1P   ON3
	  AA[p].atom[2].fullatom_type = 15; // 69; // OOC    O2P   ON3
	  AA[p].atom[3].fullatom_type = 14; // 68; // ONH2   O5*   ON2
	  AA[p].atom[4].fullatom_type = 4; // 61; // CH2    C5*   CN8 (CN8B)
	  AA[p].atom[5].fullatom_type = 3; // 60; // CH1    C4*   CN7
	  AA[p].atom[6].fullatom_type = 13; // 72; // OH     O4*   ON6 (ON6B)
	  AA[p].atom[7].fullatom_type = 3; // 60; // CH1    C3*   CN7
	  AA[p].atom[8].fullatom_type = 14; // 68; // ONH2   O3*   ON2
	  AA[p].atom[9].fullatom_type = 4; // 60; // CH2    C2*   CN7 (CN7B)
	  AA[p].atom[10].fullatom_type = 13; // 67; // OH     O2*   ON1
	  AA[p].atom[11].fullatom_type = 3; // 60; // CH1    C1*   CN7 (CN7B)
	  AA[p].atom[12].fullatom_type = 8; // 64; // Nhis   N1    NN2 (NN2B)
	  AA[p].atom[13].fullatom_type = 19; // 55; // CObb   C2    CN1 (CN1T)
	  AA[p].atom[14].fullatom_type = 20; // 67; // OCbb   O2    ON1
	  AA[p].atom[15].fullatom_type = 7; // 64; // Ntrp   N3    NN2 (NN2U)
	  AA[p].atom[16].fullatom_type = 19; // 55; // CObb   C4    CN1
	  AA[p].atom[17].fullatom_type = 20; // 67; // OCbb   O4    ON1
	  AA[p].atom[18].fullatom_type = 6; // 57; // aroC   C5    CN3
	  AA[p].atom[19].fullatom_type = 6; // 57; // aroC   C6    CN3
	  AA[p].atom[20].fullatom_type = 23; // 74; // Hapo  1H5*   HNP
	  AA[p].atom[21].fullatom_type = 23; // 74; // Hapo  2H5*   HNP
	  AA[p].atom[22].fullatom_type = 23; // 74; // Hapo   H4*   HNP
	  AA[p].atom[23].fullatom_type = 23; // 74; // Hapo   H3*   HNP
	  AA[p].atom[24].fullatom_type = 23; // 74; // Hapo  1H2*   HNP
	  AA[p].atom[25].fullatom_type = 22; // 73; // Hpol  2HO*   HP
	  AA[p].atom[26].fullatom_type = 23; // 74; // Hapo   H1*   HNP
	  AA[p].atom[27].fullatom_type = 22; // 73; // Hpol   H3    HP
	  AA[p].atom[28].fullatom_type = 24; // 74; // Haro   H5    HNP
	  AA[p].atom[29].fullatom_type = 24; // 74; // Haro   H6    HNP
	  break;
  }
  //jjh intra residue bonding
  AA[p].atom[0].nbonded_neighbors = 3; // P
  AA[p].atom[0].bonded_neighbor[0] = 1; // P--O1P
  AA[p].atom[0].bonded_neighbor[1] = 2; // P--O2P
  AA[p].atom[0].bonded_neighbor[2] = 3; // P--O5*
  AA[p].atom[1].nbonded_neighbors = 1; // O1P
  AA[p].atom[1].bonded_neighbor[0] = 0; // O1P--P
  AA[p].atom[2].nbonded_neighbors = 1; // O2P
  AA[p].atom[2].bonded_neighbor[0] = 0; // O2P--P
  AA[p].atom[3].nbonded_neighbors = 2; // O5*
  AA[p].atom[3].bonded_neighbor[0] = 0; // O5*--P
  AA[p].atom[3].bonded_neighbor[1] = 4; // O5*--C5*
  AA[p].atom[4].nbonded_neighbors = 4; // C5*
  AA[p].atom[4].bonded_neighbor[0] = 3; // C5*--O5*
  AA[p].atom[4].bonded_neighbor[1] = 5; // C5*--C4*
  AA[p].atom[4].bonded_neighbor[2] = 20; // C5*--1H5*
  AA[p].atom[4].bonded_neighbor[3] = 21; // C5*--2H5*
  AA[p].atom[5].nbonded_neighbors = 4; // C4*
  AA[p].atom[5].bonded_neighbor[0] = 4; // C4*--C5*
  AA[p].atom[5].bonded_neighbor[1] = 6; // C4*--O4*
  AA[p].atom[5].bonded_neighbor[2] = 7; // C4*--C3*
  AA[p].atom[5].bonded_neighbor[3] = 22; // C4*--H4*
  AA[p].atom[6].nbonded_neighbors = 2; // O4*
  AA[p].atom[6].bonded_neighbor[0] = 5; // O4*--C4*
  AA[p].atom[6].bonded_neighbor[1] = 11; // O4*--C1*
  AA[p].atom[7].nbonded_neighbors = 4; // C3*
  AA[p].atom[7].bonded_neighbor[0] = 5; // C3*--C4*
  AA[p].atom[7].bonded_neighbor[1] = 8; // C3*--O3*
  AA[p].atom[7].bonded_neighbor[2] = 9; // C3*--C2*
  AA[p].atom[7].bonded_neighbor[3] = 23; // C3*--H3*
  AA[p].atom[8].nbonded_neighbors = 1; // O3*
  AA[p].atom[8].bonded_neighbor[0] = 7; // O3*--C3*
  AA[p].atom[9].nbonded_neighbors = 4; // C2*
  AA[p].atom[9].bonded_neighbor[0] = 7; // C2*--C3*
  AA[p].atom[9].bonded_neighbor[1] = 10; // C2*--O2*
  AA[p].atom[9].bonded_neighbor[2] = 11; // C2*--C1*
  AA[p].atom[9].bonded_neighbor[3] = 24; // C2*--1H2*
  AA[p].atom[10].nbonded_neighbors = 2; // O2*
  AA[p].atom[10].bonded_neighbor[0] = 9; // O2*--C2*
  AA[p].atom[10].bonded_neighbor[1] = 25; // O2*--2HO*
  AA[p].atom[11].nbonded_neighbors = 4; // C1*
  AA[p].atom[11].bonded_neighbor[0] = 6; // C1*--O4*
  AA[p].atom[11].bonded_neighbor[1] = 9; // C1*--C2*
  AA[p].atom[11].bonded_neighbor[2] = 12; // C1*--N1
  AA[p].atom[11].bonded_neighbor[3] = 26; // C1*--H1*
  AA[p].atom[12].nbonded_neighbors = 3; // N1
  AA[p].atom[12].bonded_neighbor[0] = 11; // N1--C1*
  AA[p].atom[12].bonded_neighbor[1] = 13; // N1--C2
  AA[p].atom[12].bonded_neighbor[2] = 19; // N1--C6
  AA[p].atom[13].nbonded_neighbors = 3; // C2
  AA[p].atom[13].bonded_neighbor[0] = 12; // C2--N1
  AA[p].atom[13].bonded_neighbor[1] = 14; // C2--O2
  AA[p].atom[13].bonded_neighbor[2] = 15; // C2--N3
  AA[p].atom[14].nbonded_neighbors = 1; // O2
  AA[p].atom[14].bonded_neighbor[0] = 13; // O2--C2
  AA[p].atom[15].nbonded_neighbors = 3; // N3
  AA[p].atom[15].bonded_neighbor[0] = 13; // N3--C2
  AA[p].atom[15].bonded_neighbor[1] = 16; // N3--C4
  AA[p].atom[15].bonded_neighbor[2] = 27; // N3--H3
  AA[p].atom[16].nbonded_neighbors = 3; // C4
  AA[p].atom[16].bonded_neighbor[0] = 15; // C4--N3
  AA[p].atom[16].bonded_neighbor[1] = 17; // C4--O4
  AA[p].atom[16].bonded_neighbor[2] = 18; // C4--C5
  AA[p].atom[17].nbonded_neighbors = 1; // O4
  AA[p].atom[17].bonded_neighbor[0] = 16; // O4--C4
  AA[p].atom[18].nbonded_neighbors = 3; // C5
  AA[p].atom[18].bonded_neighbor[0] = 16; // C5--C4
  AA[p].atom[18].bonded_neighbor[1] = 19; // C5--C6
  AA[p].atom[18].bonded_neighbor[2] = 28; // C5--H5
  AA[p].atom[19].nbonded_neighbors = 3; // C6
  AA[p].atom[19].bonded_neighbor[0] = 12; // C6--N1
  AA[p].atom[19].bonded_neighbor[1] = 18; // C6--C5
  AA[p].atom[19].bonded_neighbor[2] = 29; // C6--H6
  AA[p].atom[20].nbonded_neighbors = 1; //1H5*
  AA[p].atom[20].bonded_neighbor[0] = 4; //1H5*--C5*
  AA[p].atom[21].nbonded_neighbors = 1; //2H5*
  AA[p].atom[21].bonded_neighbor[0] = 4; //2H5*--C5*
  AA[p].atom[22].nbonded_neighbors = 1; // H4*
  AA[p].atom[22].bonded_neighbor[0] = 5; // H4*--C4*
  AA[p].atom[23].nbonded_neighbors = 1; // H3*
  AA[p].atom[23].bonded_neighbor[0] = 7; // H3*--C3*
  AA[p].atom[24].nbonded_neighbors = 1; //1H2*
  AA[p].atom[24].bonded_neighbor[0] = 9; //1H2*--C2*
  AA[p].atom[25].nbonded_neighbors = 1; //2HO*
  AA[p].atom[25].bonded_neighbor[0] = 10; //2HO*--O2*
  AA[p].atom[26].nbonded_neighbors = 1; // H1*
  AA[p].atom[26].bonded_neighbor[0] = 11; // H1*--C1*
  AA[p].atom[27].nbonded_neighbors = 1; // H3
  AA[p].atom[27].bonded_neighbor[0] = 15; // H3--N3
  AA[p].atom[28].nbonded_neighbors = 1; // H5
  AA[p].atom[28].bonded_neighbor[0] = 18; // H5--C5
  AA[p].atom[29].nbonded_neighbors = 1; // H6
  AA[p].atom[29].bonded_neighbor[0] = 19; // H6--C6

  for ( k = 0; k < 3; k++ )
    AA[p].atom[k].hastemplate = true;

  for ( k = 3; k < 20; k++ )
    AA[p].atom[k].hastemplate = false;

  AA[p].atom[10].hastemplate = true;

  for ( k = 20; k < 30; k++ )
    AA[p].atom[k].hastemplate = true;


  //bk   template atoms used for placing atoms rURA

  //bk   template for building P
  AA[p].atom[0].ta[0] = 3; //   O5*
  AA[p].atom[0].ta[1] = 4; //   C5*
  AA[p].atom[0].ta[2] = 5; //   C4*

  //bk   template for building O1P
  AA[p].atom[1].ta[0] = 0; //   P
  AA[p].atom[1].ta[1] = 4; //   C5*
  AA[p].atom[1].ta[2] = 3; //   O5*

  //bk   template for building  O2P
  AA[p].atom[2].ta[0] = 0; //   P
  AA[p].atom[2].ta[1] = 1; //   O1P
  AA[p].atom[2].ta[2] = 3; //   O5*

  //bk   template for building  O2*
  AA[p].atom[10].ta[0] = 9; //   C2*
  AA[p].atom[10].ta[1] = 7; //   C3*
  AA[p].atom[10].ta[2] = 11; //   C1*

  //bk   template for building 1H5*
  AA[p].atom[20].ta[0] = 4; //   C5*
  AA[p].atom[20].ta[1] = 3; //   O5*
  AA[p].atom[20].ta[2] = 5; //   C4*

  //bk   template for building 2H5*
  AA[p].atom[21].ta[0] = 4; //   C5*
  AA[p].atom[21].ta[1] = 3; //   O5*
  AA[p].atom[21].ta[2] = 5; //   C4*

  //bk   template for building H4*
  AA[p].atom[22].ta[0] = 5; //   C4*
  AA[p].atom[22].ta[1] = 6; //   O4*
  AA[p].atom[22].ta[2] = 4; //   C5*

  //bk   template for building H3*
  AA[p].atom[23].ta[0] = 7; //   C3*
  AA[p].atom[23].ta[1] = 9; //   C2*
  AA[p].atom[23].ta[2] = 5; //   C4*

  //bk   template for building 1H2*
  AA[p].atom[24].ta[0] = 9; //   C2*
  AA[p].atom[24].ta[1] = 7; //   C3*
  AA[p].atom[24].ta[2] = 11; //   C1*

  //bk   template for building 2HO*
  AA[p].atom[25].ta[0] = 10; //   C2*
  AA[p].atom[25].ta[1] = 9; //   C3*
  AA[p].atom[25].ta[2] = 11; //   C1*

  //bk   template for building H1*
  AA[p].atom[26].ta[0] = 11; //   C1*
  AA[p].atom[26].ta[1] = 9; //   C2*
  AA[p].atom[26].ta[2] = 6; //   O4*

  //bk   template for building H3
  AA[p].atom[27].ta[0] = 15; //   N3
  AA[p].atom[27].ta[1] = 13; //   C2
  AA[p].atom[27].ta[2] = 16; //   C4

  //bk   template for building H5
  AA[p].atom[28].ta[0] = 18; //   C5
  AA[p].atom[28].ta[1] = 19; //   C6
  AA[p].atom[28].ta[2] = 16; //   C4

  //bk   template for building H6
  AA[p].atom[29].ta[0] = 19; //   C6
  AA[p].atom[29].ta[1] = 18; //   C5
  AA[p].atom[29].ta[2] = 16; //   C4


  //bk   four atoms that define chi angle  1
  AA[p].chi_atoms[0] [0] = 11; //   C1*
  AA[p].chi_atoms[0] [1] = 9; //   C2*
  AA[p].chi_atoms[0] [2] = 10; //   O2*
  AA[p].chi_atoms[0] [3] = 25; //   2HO*


  //bk   template coordinates for the amino acid
  //jjh new template from Kevin Karplus

  AA[p].atom[0].icoor[0] = 52.1620; //   P
  AA[p].atom[0].icoor[1] = 16.6620; //   P
  AA[p].atom[0].icoor[2] = 24.8610; //   P
  AA[p].atom[1].icoor[0] = 51.4030; //   O1P
  AA[p].atom[1].icoor[1] = 16.7770; //   O1P
  AA[p].atom[1].icoor[2] = 26.1290; //   O1P
  AA[p].atom[2].icoor[0] = 51.4500; //   O2P
  AA[p].atom[2].icoor[1] = 16.2990; //   O2P
  AA[p].atom[2].icoor[2] = 23.6140; //   O2P
  AA[p].atom[3].icoor[0] = 53.1210; //   O5*
  AA[p].atom[3].icoor[1] = 17.9070; //   O5*
  AA[p].atom[3].icoor[2] = 24.6510; //   O5*
  AA[p].atom[4].icoor[0] = 54.1200; //   C5*
  AA[p].atom[4].icoor[1] = 17.9070; //   C5*
  AA[p].atom[4].icoor[2] = 25.6740; //   C5*
  AA[p].atom[5].icoor[0] = 55.2470; //   C4*
  AA[p].atom[5].icoor[1] = 18.8460; //   C4*
  AA[p].atom[5].icoor[2] = 25.2810; //   C4*
  AA[p].atom[6].icoor[0] = 55.5460; //   O4*
  AA[p].atom[6].icoor[1] = 18.6710; //   O4*
  AA[p].atom[6].icoor[2] = 23.8770; //   O4*
  AA[p].atom[7].icoor[0] = 54.9300; //   C3*
  AA[p].atom[7].icoor[1] = 20.3420; //   C3*
  AA[p].atom[7].icoor[2] = 25.4220; //   C3*
  AA[p].atom[8].icoor[0] = 56.0750; //   O3*
  AA[p].atom[8].icoor[1] = 21.0620; //   O3*
  AA[p].atom[8].icoor[2] = 25.8960; //   O3*
  AA[p].atom[9].icoor[0] = 54.6180; //   C2*
  AA[p].atom[9].icoor[1] = 20.7580; //   C2*
  AA[p].atom[9].icoor[2] = 24.0200; //   C2*
  AA[p].atom[10].icoor[0] = 54.9030; //   O2*
  AA[p].atom[10].icoor[1] = 22.1140; //   O2*
  AA[p].atom[10].icoor[2] = 23.8470; //   O2*
  AA[p].atom[11].icoor[0] = 55.6040; //   C1*
  AA[p].atom[11].icoor[1] = 19.9410; //   C1*
  AA[p].atom[11].icoor[2] = 23.2710; //   C1*
  AA[p].atom[12].icoor[0] = 55.2750; //   N1
  AA[p].atom[12].icoor[1] = 19.9060; //   N1
  AA[p].atom[12].icoor[2] = 21.8410; //   N1
  AA[p].atom[13].icoor[0] = 56.1240; //   C2
  AA[p].atom[13].icoor[1] = 20.4660; //   C2
  AA[p].atom[13].icoor[2] = 20.9500; //   C2
  AA[p].atom[14].icoor[0] = 57.1940; //   O2
  AA[p].atom[14].icoor[1] = 20.9800; //   O2
  AA[p].atom[14].icoor[2] = 21.2820; //   O2
  AA[p].atom[15].icoor[0] = 55.6650; //   N3
  AA[p].atom[15].icoor[1] = 20.4090; //   N3
  AA[p].atom[15].icoor[2] = 19.6550; //   N3
  AA[p].atom[16].icoor[0] = 54.4500; //   C4
  AA[p].atom[16].icoor[1] = 19.9360; //   C4
  AA[p].atom[16].icoor[2] = 19.1870; //   C4
  AA[p].atom[17].icoor[0] = 54.1380; //   O4
  AA[p].atom[17].icoor[1] = 20.0260; //   O4
  AA[p].atom[17].icoor[2] = 17.9970; //   O4
  AA[p].atom[18].icoor[0] = 53.6420; //   C5
  AA[p].atom[18].icoor[1] = 19.3780; //   C5
  AA[p].atom[18].icoor[2] = 20.1960; //   C5
  AA[p].atom[19].icoor[0] = 54.0980; //   C6
  AA[p].atom[19].icoor[1] = 19.3780; //   C6
  AA[p].atom[19].icoor[2] = 21.4460; //   C6
  AA[p].atom[20].icoor[0] = 54.5200; //  1H5*
  AA[p].atom[20].icoor[1] = 16.8730; //  1H5*
  AA[p].atom[20].icoor[2] = 25.7690; //  1H5*
  AA[p].atom[21].icoor[0] = 53.6730; //  2H5*
  AA[p].atom[21].icoor[1] = 18.2050; //  2H5*
  AA[p].atom[21].icoor[2] = 26.6470; //  2H5*
  AA[p].atom[22].icoor[0] = 56.1630; //   H4*
  AA[p].atom[22].icoor[1] = 18.6370; //   H4*
  AA[p].atom[22].icoor[2] = 25.8730; //   H4*
  AA[p].atom[23].icoor[0] = 54.0720; //   H3*
  AA[p].atom[23].icoor[1] = 20.5050; //   H3*
  AA[p].atom[23].icoor[2] = 26.1120; //   H3*
  AA[p].atom[24].icoor[0] = 53.5690; //  1H2*
  AA[p].atom[24].icoor[1] = 20.5270; //  1H2*
  AA[p].atom[24].icoor[2] = 23.7370; //  1H2*
  AA[p].atom[25].icoor[0] = 54.8320; //  2HO*
  AA[p].atom[25].icoor[1] = 22.3100; //  2HO*
  AA[p].atom[25].icoor[2] = 22.9080; //  2HO*
  AA[p].atom[26].icoor[0] = 56.6330; //   H1*
  AA[p].atom[26].icoor[1] = 20.3160; //   H1*
  AA[p].atom[26].icoor[2] = 23.4580; //   H1*
  AA[p].atom[27].icoor[0] = 56.2830; //   H3
  AA[p].atom[27].icoor[1] = 20.8060; //   H3
  AA[p].atom[27].icoor[2] = 18.9610; //   H3
  AA[p].atom[28].icoor[0] = 52.6590; //   H5
  AA[p].atom[28].icoor[1] = 18.9820; //   H5
  AA[p].atom[28].icoor[2] = 19.9390; //   H5
  AA[p].atom[29].icoor[0] = 53.5030; //   H6
  AA[p].atom[29].icoor[1] = 18.9580; //   H6
  AA[p].atom[29].icoor[2] = 22.2520; //   H6

  //FIN DE rURA

} // nucleotidos FIN


void change_H_names(ConventionNames opt) {

	// Alanine A A
	switch(opt)
	{
	case iupac:
		strcpy( AA[ALA].atom[5].atom_name, " H  " );
		strcpy( AA[ALA].atom[6].atom_name, " HA " );
		strcpy( AA[ALA].atom[7].atom_name, " HB1" );
		strcpy( AA[ALA].atom[8].atom_name, " HB2" );
		strcpy( AA[ALA].atom[9].atom_name, " HB3" );
		strcpy( AA[ALA].atom[11].atom_name," H2 " );
		strcpy( AA[ALA].atom[12].atom_name," H3 " );
		break;
	case pdb:
		strcpy( AA[ALA].atom[5].atom_name, " H  " );
		strcpy( AA[ALA].atom[6].atom_name, " HA " );
		strcpy( AA[ALA].atom[7].atom_name, "1HB " );
		strcpy( AA[ALA].atom[8].atom_name, "2HB " );
		strcpy( AA[ALA].atom[9].atom_name, "3HB " );
		strcpy( AA[ALA].atom[11].atom_name,"2H  " );
		strcpy( AA[ALA].atom[12].atom_name,"3H  " );
		break;
	case sybil:
		strcpy( AA[ALA].atom[5].atom_name, " H  " );
		strcpy( AA[ALA].atom[6].atom_name, " HA " );
		strcpy( AA[ALA].atom[7].atom_name, " HB1" );
		strcpy( AA[ALA].atom[8].atom_name, " HB2" );
		strcpy( AA[ALA].atom[9].atom_name, " HB3" );
		strcpy( AA[ALA].atom[11].atom_name," H2 " );
		strcpy( AA[ALA].atom[12].atom_name," H3 " );
		break;
	}

	// ASN N

	switch(opt)
	{
	case iupac:
		strcpy( AA[ASN].atom[8].atom_name,  " H  " );
		strcpy( AA[ASN].atom[9].atom_name,  "HD21" );
		strcpy( AA[ASN].atom[10].atom_name, "HD22" );
		strcpy( AA[ASN].atom[11].atom_name, " HA " );
		strcpy( AA[ASN].atom[12].atom_name, " HB2" );
		strcpy( AA[ASN].atom[13].atom_name, " HB3" );
		strcpy( AA[ASN].atom[15].atom_name, " H2 " );
		strcpy( AA[ASN].atom[16].atom_name, " H3 " );
		break;
	case  pdb:
		strcpy( AA[ASN].atom[8].atom_name,  " H  " );
		strcpy( AA[ASN].atom[9].atom_name,  "1HD2" );
		strcpy( AA[ASN].atom[10].atom_name, "2HD2" );
		strcpy( AA[ASN].atom[11].atom_name, " HA " );
		strcpy( AA[ASN].atom[12].atom_name, "1HB " );
		strcpy( AA[ASN].atom[13].atom_name, "2HB " );
		strcpy( AA[ASN].atom[15].atom_name, "2H " );
		strcpy( AA[ASN].atom[16].atom_name, "3H " );
		break;
	case sybil:
		strcpy( AA[ASN].atom[8].atom_name,  " H  " );
		strcpy( AA[ASN].atom[9].atom_name,  "HD21" );
		strcpy( AA[ASN].atom[10].atom_name, "HD22" );
		strcpy( AA[ASN].atom[11].atom_name, " HA " );
		strcpy( AA[ASN].atom[12].atom_name, " HB2" );
		strcpy( AA[ASN].atom[13].atom_name, " HB1" );
		strcpy( AA[ASN].atom[15].atom_name, " H2 " );
		strcpy( AA[ASN].atom[16].atom_name, " H3 " );
		break;
	}


	// CYS C
	switch(opt)
	{
	case iupac:
		strcpy( AA[CYS].atom[6].atom_name,  " H  " );
		strcpy( AA[CYS].atom[7].atom_name,  " HA " );
		strcpy( AA[CYS].atom[8].atom_name,  " HB2" );
		strcpy( AA[CYS].atom[9].atom_name,  " HB3" );
		strcpy( AA[CYS].atom[10].atom_name, " HG " );
		strcpy( AA[CYS].atom[11].atom_name, " OXT" );
		strcpy( AA[CYS].atom[12].atom_name, " H2 " );
		strcpy( AA[CYS].atom[13].atom_name, " H3 " );
		break;
	case  pdb:
		strcpy( AA[CYS].atom[6].atom_name, " H  " );
		strcpy( AA[CYS].atom[7].atom_name, " HA " );
		strcpy( AA[CYS].atom[8].atom_name, "1HB " );
		strcpy( AA[CYS].atom[9].atom_name, "2HB " );
		strcpy( AA[CYS].atom[10].atom_name," HG " );
		strcpy( AA[CYS].atom[11].atom_name," OXT" );
		strcpy( AA[CYS].atom[12].atom_name,"2H  " );
		strcpy( AA[CYS].atom[13].atom_name,"3H  " );
		break;
	case sybil:
		strcpy( AA[CYS].atom[6].atom_name, " H  " );
		strcpy( AA[CYS].atom[7].atom_name, " HA " );
		strcpy( AA[CYS].atom[8].atom_name, " HB2" );
		strcpy( AA[CYS].atom[9].atom_name, " HB1" );
		strcpy( AA[CYS].atom[10].atom_name," HG " );
		strcpy( AA[CYS].atom[11].atom_name," OXT" );
		strcpy( AA[CYS].atom[12].atom_name," H2 " );
		strcpy( AA[CYS].atom[13].atom_name," H3 " );
		break;
	}


	// CYM C NEGATIVE CYS
	switch(opt)
	{
	case iupac:
		strcpy( AA[CYM].atom[6].atom_name,  " H  " );
		strcpy( AA[CYM].atom[7].atom_name,  " HA " );
		strcpy( AA[CYM].atom[8].atom_name,  " HB2" );
		strcpy( AA[CYM].atom[9].atom_name,  " HB3" );
		strcpy( AA[CYM].atom[10].atom_name, " OXT" );
		strcpy( AA[CYM].atom[11].atom_name, " H2 " );
		strcpy( AA[CYM].atom[12].atom_name, " H3 " );
	case  pdb:
		strcpy( AA[CYM].atom[6].atom_name, " H  " );
		strcpy( AA[CYM].atom[7].atom_name, " HA " );
		strcpy( AA[CYM].atom[8].atom_name, "1HB " );
		strcpy( AA[CYM].atom[9].atom_name, "2HB " );
		strcpy( AA[CYM].atom[10].atom_name," OXT" );
		strcpy( AA[CYM].atom[11].atom_name,"2H  " );
		strcpy( AA[CYM].atom[12].atom_name,"3H  " );
		break;
	case sybil:
		strcpy( AA[CYM].atom[6].atom_name, " H  " );
		strcpy( AA[CYM].atom[7].atom_name, " HA " );
		strcpy( AA[CYM].atom[8].atom_name, " HB2" );
		strcpy( AA[CYM].atom[9].atom_name, " HB1" );
		strcpy( AA[CYM].atom[10].atom_name," OXT" );
		strcpy( AA[CYM].atom[11].atom_name," H2 " );
		strcpy( AA[CYM].atom[12].atom_name," H3 " );
		break;
	}

	// CYX C SS-bonded CYS
	switch(opt)
	{
	case iupac:
		strcpy( AA[CYX].atom[6].atom_name,  " H  " );
		strcpy( AA[CYX].atom[7].atom_name,  " HA " );
		strcpy( AA[CYX].atom[8].atom_name,  " HB2" );
		strcpy( AA[CYX].atom[9].atom_name,  " HB3" );
		strcpy( AA[CYX].atom[10].atom_name, " OXT" );
		strcpy( AA[CYX].atom[11].atom_name, " H2 " );
		strcpy( AA[CYX].atom[12].atom_name, " H3 " );
	case  pdb:
		strcpy( AA[CYX].atom[6].atom_name, " H  " );
		strcpy( AA[CYX].atom[7].atom_name, " HA " );
		strcpy( AA[CYX].atom[8].atom_name, "1HB " );
		strcpy( AA[CYX].atom[9].atom_name, "2HB " );
		strcpy( AA[CYX].atom[10].atom_name," OXT" );
		strcpy( AA[CYX].atom[11].atom_name,"2H  " );
		strcpy( AA[CYX].atom[12].atom_name,"3H  " );
		break;
	case sybil:
		strcpy( AA[CYX].atom[6].atom_name, " H  " );
		strcpy( AA[CYX].atom[7].atom_name, " HA " );
		strcpy( AA[CYX].atom[8].atom_name, " HB2" );
		strcpy( AA[CYX].atom[9].atom_name, " HB1" );
		strcpy( AA[CYX].atom[10].atom_name," OXT" );
		strcpy( AA[CYX].atom[11].atom_name," H2 " );
		strcpy( AA[CYX].atom[12].atom_name," H3 " );
		break;
	}

	// ASP D
	switch(opt)
	{
	case iupac:
		strcpy( AA[ASP].atom[8].atom_name,  " H  " );
		strcpy( AA[ASP].atom[9].atom_name,  " HA " );
		strcpy( AA[ASP].atom[10].atom_name, " HB2" );
		strcpy( AA[ASP].atom[11].atom_name, " HB3" );
		strcpy( AA[ASP].atom[13].atom_name, " H2 " );
		strcpy( AA[ASP].atom[14].atom_name, " H3 " );
		break;
	case  pdb:
		strcpy( AA[ASP].atom[8].atom_name,  " H  " );
		strcpy( AA[ASP].atom[9].atom_name,  " HA " );
		strcpy( AA[ASP].atom[10].atom_name, "1HB " );
		strcpy( AA[ASP].atom[11].atom_name, "2HB " );
		strcpy( AA[ASP].atom[13].atom_name, "2H  " );
		strcpy( AA[ASP].atom[14].atom_name, "3H  " );
		break;
	case sybil:
		strcpy( AA[ASP].atom[8].atom_name,  " H  " );
		strcpy( AA[ASP].atom[9].atom_name,  " HA " );
		strcpy( AA[ASP].atom[10].atom_name, " HB2" );
		strcpy( AA[ASP].atom[11].atom_name, " HB1" );
		strcpy( AA[ASP].atom[13].atom_name, " H2 " );
		strcpy( AA[ASP].atom[14].atom_name, " H3 " );
		break;
	}

	// ASH D Neutral ASP
	switch(opt)
	{
	case iupac:
		strcpy( AA[ASH].atom[8].atom_name,  " H  " );
		strcpy( AA[ASH].atom[9].atom_name,  " HA " );
		strcpy( AA[ASH].atom[10].atom_name, " HB2" );
		strcpy( AA[ASH].atom[11].atom_name, " HB3" );
		strcpy( AA[ASH].atom[12].atom_name, " HD2" );
		strcpy( AA[ASH].atom[14].atom_name, " H2 " );
		strcpy( AA[ASH].atom[15].atom_name, " H3 " );
		break;
	case  pdb:
		strcpy( AA[ASH].atom[8].atom_name,  " H  " );
		strcpy( AA[ASH].atom[9].atom_name,  " HA " );
		strcpy( AA[ASH].atom[10].atom_name, "1HB " );
		strcpy( AA[ASH].atom[11].atom_name, "2HB " );
		strcpy( AA[ASH].atom[12].atom_name, " HD2" );
		strcpy( AA[ASH].atom[14].atom_name, "2H  " );
		strcpy( AA[ASH].atom[15].atom_name, "3H  " );
		break;
	case sybil:
		strcpy( AA[ASH].atom[8].atom_name,  " H  " );
		strcpy( AA[ASH].atom[9].atom_name,  " HA " );
		strcpy( AA[ASH].atom[10].atom_name, " HB2" );
		strcpy( AA[ASH].atom[11].atom_name, " HB1" );
		strcpy( AA[ASH].atom[12].atom_name, " HD2" );
		strcpy( AA[ASH].atom[14].atom_name, " H2 " );
		strcpy( AA[ASH].atom[15].atom_name, " H3 " );
		break;
	}

	// GLU E
	switch(opt)
	{
	case iupac:
		strcpy( AA[GLU].atom[9].atom_name,  " H  " );
		strcpy( AA[GLU].atom[10].atom_name, " HA " );
		strcpy( AA[GLU].atom[11].atom_name, " HB2" );
		strcpy( AA[GLU].atom[12].atom_name, " HB3" );
		strcpy( AA[GLU].atom[13].atom_name, " HG2" );
		strcpy( AA[GLU].atom[14].atom_name, " HG3" );
		strcpy( AA[GLU].atom[15].atom_name, " OXT" );
		strcpy( AA[GLU].atom[16].atom_name, " H2 " );
		strcpy( AA[GLU].atom[17].atom_name, " H3 " );
		break;
	case  pdb:
		strcpy( AA[GLU].atom[9].atom_name,  " H  " );
		strcpy( AA[GLU].atom[10].atom_name, " HA " );
		strcpy( AA[GLU].atom[11].atom_name, "1HB " );
		strcpy( AA[GLU].atom[12].atom_name, "2HB " );
		strcpy( AA[GLU].atom[13].atom_name, "1HG " );
		strcpy( AA[GLU].atom[14].atom_name, "2HG " );
		strcpy( AA[GLU].atom[15].atom_name, " OXT" );
		strcpy( AA[GLU].atom[16].atom_name, "2H  " );
		strcpy( AA[GLU].atom[17].atom_name, "3H  " );
		break;
	case sybil:
		strcpy( AA[GLU].atom[9].atom_name,  " H  " );
		strcpy( AA[GLU].atom[10].atom_name, " HA " );
		strcpy( AA[GLU].atom[11].atom_name, " HB2" );
		strcpy( AA[GLU].atom[12].atom_name, " HB1" );
		strcpy( AA[GLU].atom[13].atom_name, " HG2" );
		strcpy( AA[GLU].atom[14].atom_name, " HG1" );
		strcpy( AA[GLU].atom[15].atom_name, " OXT" );
		strcpy( AA[GLU].atom[16].atom_name, " H2 " );
		strcpy( AA[GLU].atom[17].atom_name, " H3 " );
		break;
	}

	// GLH E
	switch(opt)
	{
	case iupac:
		strcpy( AA[GLH].atom[9].atom_name,  " H  " );
		strcpy( AA[GLH].atom[10].atom_name, " HA " );
		strcpy( AA[GLH].atom[11].atom_name, " HB2" );
		strcpy( AA[GLH].atom[12].atom_name, " HB3" );
		strcpy( AA[GLH].atom[13].atom_name, " HG2" );
		strcpy( AA[GLH].atom[14].atom_name, " HG1" );
		strcpy( AA[GLH].atom[15].atom_name, " HE2" );//NEW ADDED
		strcpy( AA[GLH].atom[16].atom_name, " OXT" );
		strcpy( AA[GLH].atom[17].atom_name, " H2 " );
		strcpy( AA[GLH].atom[18].atom_name, " H3 " );
		break;
	case  pdb:
		strcpy( AA[GLH].atom[9].atom_name,  " H  " );
		strcpy( AA[GLH].atom[10].atom_name, " HA " );
		strcpy( AA[GLH].atom[11].atom_name, "1HB " );
		strcpy( AA[GLH].atom[12].atom_name, "2HB " );
		strcpy( AA[GLH].atom[13].atom_name, "1HG " );
		strcpy( AA[GLH].atom[14].atom_name, "2HG " );
		strcpy( AA[GLH].atom[15].atom_name, "2HE " );//NEW ADDED
		strcpy( AA[GLH].atom[16].atom_name, " OXT" );
		strcpy( AA[GLH].atom[17].atom_name, "2H  " );
		strcpy( AA[GLH].atom[18].atom_name, "3H  " );
		break;
	case sybil:
		strcpy( AA[GLH].atom[9].atom_name,  " H  " );
		strcpy( AA[GLH].atom[10].atom_name, " HA " );
		strcpy( AA[GLH].atom[11].atom_name, " HB2" );
		strcpy( AA[GLH].atom[12].atom_name, " HB3" );
		strcpy( AA[GLH].atom[13].atom_name, " HG2" );
		strcpy( AA[GLH].atom[14].atom_name, " HG1" );
		strcpy( AA[GLH].atom[15].atom_name, " HE2" );//NEW ADDED
		strcpy( AA[GLH].atom[16].atom_name, " OXT" );
		strcpy( AA[GLH].atom[17].atom_name, " H2 " );
		strcpy( AA[GLH].atom[18].atom_name, " H3 " );
		break;
	}


	// PHE F
	switch(opt)
	{
	case iupac:
		strcpy( AA[PHE].atom[11].atom_name, " H  " );
		strcpy( AA[PHE].atom[12].atom_name, " HD1" );
		strcpy( AA[PHE].atom[13].atom_name, " HE1" );
		strcpy( AA[PHE].atom[14].atom_name, " HZ " );
		strcpy( AA[PHE].atom[15].atom_name, " HE2" );
		strcpy( AA[PHE].atom[16].atom_name, " HD2" );
		strcpy( AA[PHE].atom[17].atom_name, " HA " );
		strcpy( AA[PHE].atom[18].atom_name, " HB2" );
		strcpy( AA[PHE].atom[19].atom_name, " HB3" );
		strcpy( AA[PHE].atom[21].atom_name, " H2 " );
		strcpy( AA[PHE].atom[22].atom_name, " H3 " );
		break;
	case  pdb:
		strcpy( AA[PHE].atom[11].atom_name, " H  " );
		strcpy( AA[PHE].atom[12].atom_name, "1HD " );
		strcpy( AA[PHE].atom[13].atom_name, "1HE " );
		strcpy( AA[PHE].atom[14].atom_name, " HZ " );
		strcpy( AA[PHE].atom[15].atom_name, "2HE " );
		strcpy( AA[PHE].atom[16].atom_name, "2HD " );
		strcpy( AA[PHE].atom[17].atom_name, " HA " );
		strcpy( AA[PHE].atom[18].atom_name, "1HB " );
		strcpy( AA[PHE].atom[19].atom_name, "2HB " );
		strcpy( AA[PHE].atom[21].atom_name, "2H  " );
		strcpy( AA[PHE].atom[22].atom_name, "3H  " );
		break;
	case sybil:
		strcpy( AA[PHE].atom[11].atom_name, " H  " );
		strcpy( AA[PHE].atom[12].atom_name, " HD1" );
		strcpy( AA[PHE].atom[13].atom_name, " HE1" );
		strcpy( AA[PHE].atom[14].atom_name, " HZ " );
		strcpy( AA[PHE].atom[15].atom_name, " HE2" );
		strcpy( AA[PHE].atom[16].atom_name, " HD2" );
		strcpy( AA[PHE].atom[17].atom_name, " HA " );
		strcpy( AA[PHE].atom[18].atom_name, " HB2" );
		strcpy( AA[PHE].atom[19].atom_name, " HB1" );
		strcpy( AA[PHE].atom[21].atom_name, " H2 " );
		strcpy( AA[PHE].atom[22].atom_name, " H3 " );
		break;
	}

	// GLY G

	switch(opt)
	{
	case iupac:
		strcpy( AA[GLY].atom[4].atom_name, " H  " );
		strcpy( AA[GLY].atom[5].atom_name, " HA2" );
		strcpy( AA[GLY].atom[6].atom_name, " HA3" );
		strcpy( AA[GLY].atom[8].atom_name, " H2 " );
		strcpy( AA[GLY].atom[9].atom_name, " H3 " );
		break;
	case  pdb:
		strcpy( AA[GLY].atom[4].atom_name, " H  " );
		strcpy( AA[GLY].atom[5].atom_name, "1HA " );
		strcpy( AA[GLY].atom[6].atom_name, "2HA " );
		strcpy( AA[GLY].atom[8].atom_name, "2H  " );
		strcpy( AA[GLY].atom[9].atom_name, "3H  " );
		break;
	case sybil:
		strcpy( AA[GLY].atom[4].atom_name, " H  " );
		strcpy( AA[GLY].atom[5].atom_name, " HA2" );
		strcpy( AA[GLY].atom[6].atom_name, " HA1" );
		strcpy( AA[GLY].atom[8].atom_name, " H2 " );
		strcpy( AA[GLY].atom[9].atom_name, " H3 " );
		break;
	}

	// PRO P
	switch(opt)
	{
	case iupac:
		strcpy( AA[PRO].atom[7].atom_name,  " HD2" );
		strcpy( AA[PRO].atom[8].atom_name,  " HD3" );
		strcpy( AA[PRO].atom[9].atom_name,  " HG2" );
		strcpy( AA[PRO].atom[10].atom_name, " HG3" );
		strcpy( AA[PRO].atom[11].atom_name, " HB2" );
		strcpy( AA[PRO].atom[12].atom_name, " HB3" );
		strcpy( AA[PRO].atom[13].atom_name, " HA " );
		strcpy( AA[PRO].atom[15].atom_name, " H2 " );
		strcpy( AA[PRO].atom[16].atom_name, " H3 " );
		break;
	case  pdb:
		strcpy( AA[PRO].atom[7].atom_name,  "1HD " );
		strcpy( AA[PRO].atom[8].atom_name,  "2HD " );
		strcpy( AA[PRO].atom[9].atom_name,  "1HG " );
		strcpy( AA[PRO].atom[10].atom_name, "2HG " );
		strcpy( AA[PRO].atom[11].atom_name, "1HB " );
		strcpy( AA[PRO].atom[12].atom_name, "2HB " );
		strcpy( AA[PRO].atom[13].atom_name, " HA " );
		strcpy( AA[PRO].atom[15].atom_name, "2H  " );
		strcpy( AA[PRO].atom[16].atom_name, "3H  " );
		break;
	case sybil:

		strcpy( AA[PRO].atom[7].atom_name,  " HD2" );
		strcpy( AA[PRO].atom[8].atom_name,  " HD1" );
		strcpy( AA[PRO].atom[9].atom_name,  " HG2" );
		strcpy( AA[PRO].atom[10].atom_name, " HG1" );
		strcpy( AA[PRO].atom[11].atom_name, " HB2" );
		strcpy( AA[PRO].atom[12].atom_name, " HB1" );
		strcpy( AA[PRO].atom[13].atom_name, " HA " );
		strcpy( AA[PRO].atom[15].atom_name, " H2 " );
		strcpy( AA[PRO].atom[16].atom_name, " H3 " );


		break;
	}




	//  VAL V

	switch(opt)
	{
	case iupac:
		strcpy( AA[VAL].atom[7].atom_name, " H  " );
		strcpy( AA[VAL].atom[8].atom_name, " HA " );
		strcpy( AA[VAL].atom[9].atom_name, " HB " );
		strcpy( AA[VAL].atom[10].atom_name, "HG11" );
		strcpy( AA[VAL].atom[11].atom_name, "HG12" );
		strcpy( AA[VAL].atom[12].atom_name, "HG13" );
		strcpy( AA[VAL].atom[13].atom_name, "HG21" );
		strcpy( AA[VAL].atom[14].atom_name, "HG22" );
		strcpy( AA[VAL].atom[15].atom_name, "HG23" );
		strcpy( AA[VAL].atom[17].atom_name, " H2 " );
		strcpy( AA[VAL].atom[18].atom_name, " H3 " );
		break;
	case  pdb:
		strcpy( AA[VAL].atom[7].atom_name, " H  " );
		strcpy( AA[VAL].atom[8].atom_name, " HA " );
		strcpy( AA[VAL].atom[9].atom_name, " HB " );
		strcpy( AA[VAL].atom[10].atom_name, "1HG1" );
		strcpy( AA[VAL].atom[11].atom_name, "2HG1" );
		strcpy( AA[VAL].atom[12].atom_name, "3HG1" );
		strcpy( AA[VAL].atom[13].atom_name, "1HG2" );
		strcpy( AA[VAL].atom[14].atom_name, "2HG2" );
		strcpy( AA[VAL].atom[15].atom_name, "3HG2" );

		strcpy( AA[VAL].atom[17].atom_name, "2H  " );
		strcpy( AA[VAL].atom[18].atom_name, "3H  " );
		break;
	case sybil:
		strcpy( AA[VAL].atom[7].atom_name, " H  " );
		strcpy( AA[VAL].atom[8].atom_name, " HA " );
		strcpy( AA[VAL].atom[9].atom_name, " HB " );
		strcpy( AA[VAL].atom[10].atom_name, "HG11" );
		strcpy( AA[VAL].atom[11].atom_name, "HG12" );
		strcpy( AA[VAL].atom[12].atom_name, "HG13" );
		strcpy( AA[VAL].atom[13].atom_name, "HG21" );
		strcpy( AA[VAL].atom[14].atom_name, "HG22" );
		strcpy( AA[VAL].atom[15].atom_name, "HG23" );
		strcpy( AA[VAL].atom[17].atom_name, " H2 " );
		strcpy( AA[VAL].atom[18].atom_name, " H3 " );
		break;
	}





	// ILE I

	switch(opt)
	{
	case iupac:
		strcpy( AA[ILE].atom[8].atom_name,  " H  " );
		strcpy( AA[ILE].atom[9].atom_name,  " HA " );
		strcpy( AA[ILE].atom[10].atom_name, " HB " );
		strcpy( AA[ILE].atom[11].atom_name, "HG21" );
		strcpy( AA[ILE].atom[12].atom_name, "HG22" );
		strcpy( AA[ILE].atom[13].atom_name, "HG23" );
		strcpy( AA[ILE].atom[14].atom_name, "HG12" );
		strcpy( AA[ILE].atom[15].atom_name, "HG13" );
		strcpy( AA[ILE].atom[16].atom_name, "HD11" );
		strcpy( AA[ILE].atom[17].atom_name, "HD12" );
		strcpy( AA[ILE].atom[18].atom_name, "HD13" );
		strcpy( AA[ILE].atom[20].atom_name, " H2 " );
		strcpy( AA[ILE].atom[21].atom_name, " H3 " );
		break;
	case  pdb:
		strcpy( AA[ILE].atom[8].atom_name, " H  " );
		strcpy( AA[ILE].atom[9].atom_name, " HA " );
		strcpy( AA[ILE].atom[10].atom_name, " HB " );
		strcpy( AA[ILE].atom[11].atom_name, "1HG2" );
		strcpy( AA[ILE].atom[12].atom_name, "2HG2" );
		strcpy( AA[ILE].atom[13].atom_name, "3HG2" );
		strcpy( AA[ILE].atom[14].atom_name, "1HG1" );
		strcpy( AA[ILE].atom[15].atom_name, "2HG1" );
		strcpy( AA[ILE].atom[16].atom_name, "1HD1" );
		strcpy( AA[ILE].atom[17].atom_name, "2HD1" );
		strcpy( AA[ILE].atom[18].atom_name, "3HD1" );
		strcpy( AA[ILE].atom[20].atom_name, "2H  " );
		strcpy( AA[ILE].atom[21].atom_name, "3H  " );
		break;
	case sybil:
		strcpy( AA[ILE].atom[8].atom_name,  " H  " );
		strcpy( AA[ILE].atom[9].atom_name,  " HA " );
		strcpy( AA[ILE].atom[10].atom_name, " HB " );
		strcpy( AA[ILE].atom[11].atom_name, "HG21" );
		strcpy( AA[ILE].atom[12].atom_name, "HG22" );
		strcpy( AA[ILE].atom[13].atom_name, "HG23" );
		strcpy( AA[ILE].atom[14].atom_name, "HG12" );
		strcpy( AA[ILE].atom[15].atom_name, "HG11" );
		strcpy( AA[ILE].atom[16].atom_name, "HD11" );
		strcpy( AA[ILE].atom[17].atom_name, "HD12" );
		strcpy( AA[ILE].atom[18].atom_name, "HD13" );
		strcpy( AA[ILE].atom[20].atom_name, " H2 " );
		strcpy( AA[ILE].atom[21].atom_name, " H3 " );
		break;
	}

	// HIS  H;
	switch(opt)
	{
	case iupac:
		strcpy( AA[HIS].atom[10].atom_name," H  " );
		strcpy( AA[HIS].atom[11].atom_name," HE2" );// <--- The proton
		strcpy( AA[HIS].atom[12].atom_name," HA " );
		strcpy( AA[HIS].atom[13].atom_name," HB2" );
		strcpy( AA[HIS].atom[14].atom_name," HB3" );
		strcpy( AA[HIS].atom[15].atom_name," HE1" );
		strcpy( AA[HIS].atom[16].atom_name," HD2" );
		strcpy( AA[HIS].atom[18].atom_name," H2 " );
		strcpy( AA[HIS].atom[19].atom_name," H3 " );
		break;
	case  pdb:
		strcpy( AA[HIS].atom[10].atom_name," H  " );
		strcpy( AA[HIS].atom[11].atom_name," HE2" );// <--- The proton
		strcpy( AA[HIS].atom[12].atom_name," HA " );
		strcpy( AA[HIS].atom[13].atom_name,"1HB " );
		strcpy( AA[HIS].atom[14].atom_name,"2HB " );
		strcpy( AA[HIS].atom[15].atom_name," HE1" );
		strcpy( AA[HIS].atom[16].atom_name," HD2" );
		strcpy( AA[HIS].atom[18].atom_name,"2H  " );
		strcpy( AA[HIS].atom[19].atom_name,"3H  " );
		break;
	case sybil:
		strcpy( AA[HIS].atom[10].atom_name," H  " );
		strcpy( AA[HIS].atom[11].atom_name," HE2" );// <--- The proton
		strcpy( AA[HIS].atom[12].atom_name," HA " );
		strcpy( AA[HIS].atom[13].atom_name," HB2" );
		strcpy( AA[HIS].atom[14].atom_name," HB3" );
		strcpy( AA[HIS].atom[15].atom_name," HE1" );
		strcpy( AA[HIS].atom[16].atom_name," HD2" );
		strcpy( AA[HIS].atom[18].atom_name," H2 " );
		strcpy( AA[HIS].atom[19].atom_name," H3 " );
		break;
	}

	// HIP H
	switch(opt)
	{
	case iupac:
		strcpy( AA[HIP].atom[10].atom_name," H  " );
		strcpy( AA[HIP].atom[11].atom_name," HE2" );// <--- The proton
		strcpy( AA[HIP].atom[12].atom_name," HD1" );// <--- The proton
		strcpy( AA[HIP].atom[13].atom_name," HA " );
		strcpy( AA[HIP].atom[14].atom_name," HB2" );
		strcpy( AA[HIP].atom[15].atom_name," HB3" );
		strcpy( AA[HIP].atom[16].atom_name," HE1" );
		strcpy( AA[HIP].atom[17].atom_name," HD2" );
		strcpy( AA[HIP].atom[19].atom_name," H2 " );
		strcpy( AA[HIP].atom[20].atom_name," H3 " );
		break;
	case  pdb:
		strcpy( AA[HIP].atom[10].atom_name," H  " );
		strcpy( AA[HIP].atom[11].atom_name," HE2" );// <--- The proton
		strcpy( AA[HIP].atom[12].atom_name," HD1" );// <--- The proton.
		strcpy( AA[HIP].atom[13].atom_name," HA " );
		strcpy( AA[HIP].atom[14].atom_name,"1HB " );
		strcpy( AA[HIP].atom[15].atom_name,"2HB " );
		strcpy( AA[HIP].atom[16].atom_name," HE1" );
		strcpy( AA[HIP].atom[17].atom_name," HD2" );
		strcpy( AA[HIP].atom[19].atom_name,"2H  " );
		strcpy( AA[HIP].atom[20].atom_name,"3H  " );
		break;
	case sybil:
		strcpy( AA[HIP].atom[10].atom_name," H  " );
		strcpy( AA[HIP].atom[11].atom_name," HE2" );// <--- The proton
		strcpy( AA[HIP].atom[12].atom_name," HD1" );// <--- The proton.
		strcpy( AA[HIP].atom[13].atom_name," HA " );
		strcpy( AA[HIP].atom[14].atom_name," HB2" );
		strcpy( AA[HIP].atom[15].atom_name," HB3" );
		strcpy( AA[HIP].atom[16].atom_name," HE1" );
		strcpy( AA[HIP].atom[17].atom_name," HD2" );
		strcpy( AA[HIP].atom[19].atom_name," H2 " );
		strcpy( AA[HIP].atom[20].atom_name," H3 " );
		break;
	}

	// HID
	switch(opt)
	{
	case iupac:
		strcpy( AA[HID].atom[10].atom_name," H  " );
		strcpy( AA[HID].atom[11].atom_name," HD1" );// <--- The proton
		strcpy( AA[HID].atom[12].atom_name," HA " );
		strcpy( AA[HID].atom[13].atom_name," HB2" );
		strcpy( AA[HID].atom[14].atom_name," HB3" );
		strcpy( AA[HID].atom[15].atom_name," HE1" );
		strcpy( AA[HID].atom[16].atom_name," HD2" );
		strcpy( AA[HID].atom[18].atom_name," H2 " );
		strcpy( AA[HID].atom[19].atom_name," H3 " );
		break;
	case  pdb:
		strcpy( AA[HID].atom[10].atom_name," H  " );
		strcpy( AA[HID].atom[11].atom_name," HD1" );// <--- The proton
		strcpy( AA[HID].atom[12].atom_name," HA " );
		strcpy( AA[HID].atom[13].atom_name,"1HB " );
		strcpy( AA[HID].atom[14].atom_name,"2HB " );
		strcpy( AA[HID].atom[15].atom_name," HE1" );
		strcpy( AA[HID].atom[16].atom_name," HD2" );
		strcpy( AA[HID].atom[18].atom_name,"2H  " );
		strcpy( AA[HID].atom[19].atom_name,"3H  " );
		break;
	case sybil:
		strcpy( AA[HID].atom[10].atom_name," H  " );
		strcpy( AA[HID].atom[11].atom_name," HD1" );// <--- The proton
		strcpy( AA[HID].atom[12].atom_name," HA " );
		strcpy( AA[HID].atom[13].atom_name," HB2" );
		strcpy( AA[HID].atom[14].atom_name," HB3" );
		strcpy( AA[HID].atom[15].atom_name," HE1" );
		strcpy( AA[HID].atom[16].atom_name," HD2" );
		strcpy( AA[HID].atom[18].atom_name," H2 " );
		strcpy( AA[HID].atom[19].atom_name," H3 " );
		break;
	}




	// HIE
	switch(opt)
	{
	case iupac:
		strcpy( AA[HIE].atom[10].atom_name," H  " );
		strcpy( AA[HIE].atom[11].atom_name," HE2" );// <--- The proton
		strcpy( AA[HIE].atom[12].atom_name," HA " );
		strcpy( AA[HIE].atom[13].atom_name," HB2" );
		strcpy( AA[HIE].atom[14].atom_name," HB3" );
		strcpy( AA[HIE].atom[15].atom_name," HE1" );
		strcpy( AA[HIE].atom[16].atom_name," HD2" );
		strcpy( AA[HIE].atom[18].atom_name," H2 " );
		strcpy( AA[HIE].atom[19].atom_name," H3 " );
		break;
	case  pdb:
		strcpy( AA[HIE].atom[10].atom_name," H  " );
		strcpy( AA[HIE].atom[11].atom_name," HE2" );// <--- The proton
		strcpy( AA[HIE].atom[12].atom_name," HA " );
		strcpy( AA[HIE].atom[13].atom_name,"1HB " );
		strcpy( AA[HIE].atom[14].atom_name,"2HB " );
		strcpy( AA[HIE].atom[15].atom_name," HE1" );
		strcpy( AA[HIE].atom[16].atom_name," HD2" );
		strcpy( AA[HIE].atom[18].atom_name,"2H  " );
		strcpy( AA[HIE].atom[19].atom_name,"3H  " );
		break;
	case sybil:
		strcpy( AA[HIE].atom[10].atom_name," H  " );
		strcpy( AA[HIE].atom[11].atom_name," HE2" );// <--- The proton
		strcpy( AA[HIE].atom[12].atom_name," HA " );
		strcpy( AA[HIE].atom[13].atom_name," HB2" );
		strcpy( AA[HIE].atom[14].atom_name," HB3" );
		strcpy( AA[HIE].atom[15].atom_name," HE1" );
		strcpy( AA[HIE].atom[16].atom_name," HD2" );
		strcpy( AA[HIE].atom[18].atom_name," H2 " );
		strcpy( AA[HIE].atom[19].atom_name," H3 " );
		break;
	}


	// LEU L
	switch(opt)
	{
	case iupac:
		strcpy( AA[LEU].atom[8].atom_name,  " H  " );
		strcpy( AA[LEU].atom[9].atom_name,  " HA " );
		strcpy( AA[LEU].atom[10].atom_name, " HB2" );
		strcpy( AA[LEU].atom[11].atom_name, " HB3" );
		strcpy( AA[LEU].atom[12].atom_name, " HG " );
		strcpy( AA[LEU].atom[13].atom_name, "HD11" );
		strcpy( AA[LEU].atom[14].atom_name, "HD12" );
		strcpy( AA[LEU].atom[15].atom_name, "HD13" );
		strcpy( AA[LEU].atom[16].atom_name, "HD21" );
		strcpy( AA[LEU].atom[17].atom_name, "HD22" );
		strcpy( AA[LEU].atom[18].atom_name, "HD23" );
		strcpy( AA[LEU].atom[19].atom_name, " OXT" );
		strcpy( AA[LEU].atom[20].atom_name, " H2 " );
		strcpy( AA[LEU].atom[21].atom_name, " H3 " );
		break;
	case  pdb:
		strcpy( AA[LEU].atom[8].atom_name,  " H  " );
		strcpy( AA[LEU].atom[9].atom_name,  " HA " );
		strcpy( AA[LEU].atom[10].atom_name, " HB2" );
		strcpy( AA[LEU].atom[11].atom_name, " HB3" );
		strcpy( AA[LEU].atom[12].atom_name, " HG " );
		strcpy( AA[LEU].atom[13].atom_name, "1HD1" );
		strcpy( AA[LEU].atom[14].atom_name, "2HD1" );
		strcpy( AA[LEU].atom[15].atom_name, "3HD1" );
		strcpy( AA[LEU].atom[16].atom_name, "1HD2" );
		strcpy( AA[LEU].atom[17].atom_name, "2HD2" );
		strcpy( AA[LEU].atom[18].atom_name, "3HD2" );
		strcpy( AA[LEU].atom[19].atom_name, " OXT" );
		strcpy( AA[LEU].atom[20].atom_name, "2H  " );
		strcpy( AA[LEU].atom[21].atom_name, "3H  " );
		break;
	case sybil:
		strcpy( AA[LEU].atom[8].atom_name,  " H  " );
		strcpy( AA[LEU].atom[9].atom_name,  " HA " );
		strcpy( AA[LEU].atom[10].atom_name, " HB2" );
		strcpy( AA[LEU].atom[11].atom_name, " HB1" );
		strcpy( AA[LEU].atom[12].atom_name, " HG " );
		strcpy( AA[LEU].atom[13].atom_name, "HD11" );
		strcpy( AA[LEU].atom[14].atom_name, "HD12" );
		strcpy( AA[LEU].atom[15].atom_name, "HD13" );
		strcpy( AA[LEU].atom[16].atom_name, "HD21" );
		strcpy( AA[LEU].atom[17].atom_name, "HD22" );
		strcpy( AA[LEU].atom[18].atom_name, "HD23" );
		strcpy( AA[LEU].atom[19].atom_name, " OXT" );
		strcpy( AA[LEU].atom[20].atom_name, " H2 " );
		strcpy( AA[LEU].atom[21].atom_name, " H3 " );
		break;
	}


	// MET M
	switch(opt)
	{
	case iupac:
		strcpy( AA[MET].atom[8].atom_name,  " H  " );
		strcpy( AA[MET].atom[9].atom_name,  " HA " );
		strcpy( AA[MET].atom[10].atom_name, " HB2" );
		strcpy( AA[MET].atom[11].atom_name, " HB3" );
		strcpy( AA[MET].atom[12].atom_name, " HG2" );
		strcpy( AA[MET].atom[13].atom_name, " HG3" );
		strcpy( AA[MET].atom[14].atom_name, " HE1" );
		strcpy( AA[MET].atom[15].atom_name, " HE2" );
		strcpy( AA[MET].atom[16].atom_name, " HE3" );
		strcpy( AA[MET].atom[17].atom_name, " OXT" );
		strcpy( AA[MET].atom[18].atom_name, " H2 " );
		strcpy( AA[MET].atom[19].atom_name, " H3 " );
		break;
	case  pdb:
		strcpy( AA[MET].atom[8].atom_name,  " H  " );
		strcpy( AA[MET].atom[9].atom_name,  " HA " );
		strcpy( AA[MET].atom[10].atom_name, "1HB " );
		strcpy( AA[MET].atom[11].atom_name, "2HB " );
		strcpy( AA[MET].atom[12].atom_name, "1HG " );
		strcpy( AA[MET].atom[13].atom_name, "2HG " );
		strcpy( AA[MET].atom[14].atom_name, "1HE " );
		strcpy( AA[MET].atom[15].atom_name, "2HE " );
		strcpy( AA[MET].atom[16].atom_name, "3HE " );
		strcpy( AA[MET].atom[17].atom_name, " OXT" );
		strcpy( AA[MET].atom[18].atom_name, "2H  " );
		strcpy( AA[MET].atom[19].atom_name, "3H  " );
		break;
	case sybil:
		strcpy( AA[MET].atom[8].atom_name,  " H  " );
		strcpy( AA[MET].atom[9].atom_name,  " HA " );
		strcpy( AA[MET].atom[10].atom_name, " HB2" );
		strcpy( AA[MET].atom[11].atom_name, " HB1" );
		strcpy( AA[MET].atom[12].atom_name, " HG2" );
		strcpy( AA[MET].atom[13].atom_name, " HG1" );
		strcpy( AA[MET].atom[14].atom_name, " HE1" );
		strcpy( AA[MET].atom[15].atom_name, " HE2" );
		strcpy( AA[MET].atom[16].atom_name, " HE3" );
		strcpy( AA[MET].atom[17].atom_name, " OXT" );
		strcpy( AA[MET].atom[18].atom_name, " H2  " );
		strcpy( AA[MET].atom[19].atom_name, " H3  " );
		break;
	}
	// MSE M
	switch(opt)
	{
	case iupac:
		strcpy( AA[MSE].atom[8].atom_name,  " H  " );
		strcpy( AA[MSE].atom[9].atom_name,  " HA " );
		strcpy( AA[MSE].atom[10].atom_name, " HB2" );
		strcpy( AA[MSE].atom[11].atom_name, " HB3" );
		strcpy( AA[MSE].atom[12].atom_name, " HG2" );
		strcpy( AA[MSE].atom[13].atom_name, " HG3" );
		strcpy( AA[MSE].atom[14].atom_name, " HE1" );
		strcpy( AA[MSE].atom[15].atom_name, " HE2" );
		strcpy( AA[MSE].atom[16].atom_name, " HE3" );
		strcpy( AA[MSE].atom[17].atom_name, " OXT" );
		strcpy( AA[MSE].atom[18].atom_name, " H2 " );
		strcpy( AA[MSE].atom[19].atom_name, " H3 " );
		break;
	case  pdb:
		strcpy( AA[MSE].atom[8].atom_name,  " H  " );
		strcpy( AA[MSE].atom[9].atom_name,  " HA " );
		strcpy( AA[MSE].atom[10].atom_name, "1HB " );
		strcpy( AA[MSE].atom[11].atom_name, "2HB " );
		strcpy( AA[MSE].atom[12].atom_name, "1HG " );
		strcpy( AA[MSE].atom[13].atom_name, "2HG " );
		strcpy( AA[MSE].atom[14].atom_name, "1HE " );
		strcpy( AA[MSE].atom[15].atom_name, "2HE " );
		strcpy( AA[MSE].atom[16].atom_name, "3HE " );
		strcpy( AA[MSE].atom[17].atom_name, " OXT" );
		strcpy( AA[MSE].atom[18].atom_name, "2H  " );
		strcpy( AA[MSE].atom[19].atom_name, "3H  " );
		break;
	case sybil:
		strcpy( AA[MSE].atom[8].atom_name,  " H  " );
		strcpy( AA[MSE].atom[9].atom_name,  " HA " );
		strcpy( AA[MSE].atom[10].atom_name, " HB2" );
		strcpy( AA[MSE].atom[11].atom_name, " HB1" );
		strcpy( AA[MSE].atom[12].atom_name, " HG2" );
		strcpy( AA[MSE].atom[13].atom_name, " HG1" );
		strcpy( AA[MSE].atom[14].atom_name, " HE1" );
		strcpy( AA[MSE].atom[15].atom_name, " HE2" );
		strcpy( AA[MSE].atom[16].atom_name, " HE3" );
		strcpy( AA[MSE].atom[17].atom_name, " OXT" );
		strcpy( AA[MSE].atom[18].atom_name, " H2 " );
		strcpy( AA[MSE].atom[19].atom_name, " H3 " );
		break;
	}


	// GLN Q
	switch(opt)
	{
	case iupac:
		strcpy( AA[GLN].atom[9].atom_name,  " H  " );
		strcpy( AA[GLN].atom[10].atom_name, "HE21" );
		strcpy( AA[GLN].atom[11].atom_name, "HE22" );
		strcpy( AA[GLN].atom[12].atom_name, " HA " );
		strcpy( AA[GLN].atom[13].atom_name, " HB2" );
		strcpy( AA[GLN].atom[14].atom_name, " HB3" );
		strcpy( AA[GLN].atom[15].atom_name, " HG2" );
		strcpy( AA[GLN].atom[16].atom_name, " HG3" );
		strcpy( AA[GLN].atom[17].atom_name, " OXT" );
		strcpy( AA[GLN].atom[18].atom_name, " H2 " );
		strcpy( AA[GLN].atom[19].atom_name, " H3 " );
		break;
	case  pdb:
		strcpy( AA[GLN].atom[9].atom_name,  " H  " );
		strcpy( AA[GLN].atom[10].atom_name, "1HE2" );
		strcpy( AA[GLN].atom[11].atom_name, "2HE2" );
		strcpy( AA[GLN].atom[12].atom_name, " HA " );
		strcpy( AA[GLN].atom[13].atom_name, "1HB " );
		strcpy( AA[GLN].atom[14].atom_name, "2HB " );
		strcpy( AA[GLN].atom[15].atom_name, "1HG " );
		strcpy( AA[GLN].atom[16].atom_name, "2HG " );
		strcpy( AA[GLN].atom[17].atom_name, " OXT" );
		strcpy( AA[GLN].atom[18].atom_name, "2H  " );
		strcpy( AA[GLN].atom[19].atom_name, "3H  " );
		break;
	case sybil:
		strcpy( AA[GLN].atom[9].atom_name,  " H  " );
		strcpy( AA[GLN].atom[10].atom_name, "HE21" );
		strcpy( AA[GLN].atom[11].atom_name, "HE22" );
		strcpy( AA[GLN].atom[12].atom_name, " HA " );
		strcpy( AA[GLN].atom[13].atom_name, " HB2" );
		strcpy( AA[GLN].atom[14].atom_name, " HB1" );
		strcpy( AA[GLN].atom[15].atom_name, " HG2" );
		strcpy( AA[GLN].atom[16].atom_name, " HG1" );
		strcpy( AA[GLN].atom[17].atom_name, " OXT" );
		strcpy( AA[GLN].atom[18].atom_name, " H2 " );
		strcpy( AA[GLN].atom[19].atom_name, " H3 " );
		break;
	}


	// ARG R
	switch(opt)
	{
	case iupac:
		strcpy( AA[ARG].atom[11].atom_name, " H  " );
		strcpy( AA[ARG].atom[12].atom_name, "HH11" );
		strcpy( AA[ARG].atom[13].atom_name, "HH12" );
		strcpy( AA[ARG].atom[14].atom_name, "HH21" );
		strcpy( AA[ARG].atom[15].atom_name, "HH22" );
		strcpy( AA[ARG].atom[16].atom_name, " HE " );
		strcpy( AA[ARG].atom[17].atom_name, " HA " );
		strcpy( AA[ARG].atom[18].atom_name, " HB2" );
		strcpy( AA[ARG].atom[19].atom_name, " HB3" );
		strcpy( AA[ARG].atom[20].atom_name, " HG2" );
		strcpy( AA[ARG].atom[21].atom_name, " HG3" );
		strcpy( AA[ARG].atom[22].atom_name, " HD2" );
		strcpy( AA[ARG].atom[23].atom_name, " HD3" );
		strcpy( AA[ARG].atom[24].atom_name, " OXT" );
		strcpy( AA[ARG].atom[25].atom_name, " H2 " );
		strcpy( AA[ARG].atom[26].atom_name, " H3 " );
		break;
	case  pdb:
		strcpy( AA[ARG].atom[11].atom_name, " H  " );
		strcpy( AA[ARG].atom[12].atom_name, "1HH1" );
		strcpy( AA[ARG].atom[13].atom_name, "2HH1" );
		strcpy( AA[ARG].atom[14].atom_name, "1HH2" );
		strcpy( AA[ARG].atom[15].atom_name, "2HH2" );
		strcpy( AA[ARG].atom[16].atom_name, " HE " );
		strcpy( AA[ARG].atom[17].atom_name, " HA " );
		strcpy( AA[ARG].atom[18].atom_name, "1HB " );
		strcpy( AA[ARG].atom[19].atom_name, "2HB " );
		strcpy( AA[ARG].atom[20].atom_name, "1HG " );
		strcpy( AA[ARG].atom[21].atom_name, "2HG " );
		strcpy( AA[ARG].atom[22].atom_name, "1HD " );
		strcpy( AA[ARG].atom[23].atom_name, "2HD " );
		strcpy( AA[ARG].atom[24].atom_name, " OXT" );
		strcpy( AA[ARG].atom[25].atom_name, "2H  " );
		strcpy( AA[ARG].atom[26].atom_name, "3H  " );
		break;
	case sybil:
		strcpy( AA[ARG].atom[11].atom_name, " H  " );
		strcpy( AA[ARG].atom[12].atom_name, "HH11" );
		strcpy( AA[ARG].atom[13].atom_name, "HH12" );
		strcpy( AA[ARG].atom[14].atom_name, "HH21" );
		strcpy( AA[ARG].atom[15].atom_name, "HH22" );
		strcpy( AA[ARG].atom[16].atom_name, " HE " );
		strcpy( AA[ARG].atom[17].atom_name, " HA " );
		strcpy( AA[ARG].atom[18].atom_name, " HB2" );
		strcpy( AA[ARG].atom[19].atom_name, " HB1" );
		strcpy( AA[ARG].atom[20].atom_name, " HG2" );
		strcpy( AA[ARG].atom[21].atom_name, " HG1" );
		strcpy( AA[ARG].atom[22].atom_name, " HD2" );
		strcpy( AA[ARG].atom[23].atom_name, " HD1" );
		strcpy( AA[ARG].atom[24].atom_name, " OXT" );
		strcpy( AA[ARG].atom[25].atom_name, " H2 " );
		strcpy( AA[ARG].atom[26].atom_name, " H3 " );
		break;
	}

	// SER S
	switch(opt)
	{
	case iupac:
		strcpy( AA[SER].atom[6].atom_name,  " H  " );
		strcpy( AA[SER].atom[7].atom_name,  " HG " );
		strcpy( AA[SER].atom[8].atom_name,  " HA " );
		strcpy( AA[SER].atom[9].atom_name,  " HB2" );
		strcpy( AA[SER].atom[10].atom_name, " HB3" );
		strcpy( AA[SER].atom[12].atom_name, " H2 " );
		strcpy( AA[SER].atom[13].atom_name, " H3 " );
		break;
	case  pdb:
		strcpy( AA[SER].atom[6].atom_name,  " H  " );
		strcpy( AA[SER].atom[7].atom_name,  " HG " );
		strcpy( AA[SER].atom[8].atom_name,  " HA " );
		strcpy( AA[SER].atom[9].atom_name,  "2HB " );
		strcpy( AA[SER].atom[10].atom_name, "1HB " );
		strcpy( AA[SER].atom[12].atom_name, "2H  " );
		strcpy( AA[SER].atom[13].atom_name, "3H  " );
		break;
	case sybil:
		strcpy( AA[SER].atom[6].atom_name,  " H  " );
		strcpy( AA[SER].atom[7].atom_name,  " HG " );
		strcpy( AA[SER].atom[8].atom_name,  " HA " );
		strcpy( AA[SER].atom[9].atom_name,  " HB2" );
		strcpy( AA[SER].atom[10].atom_name, " HB1" );
		strcpy( AA[SER].atom[12].atom_name, " H2 " );
		strcpy( AA[SER].atom[13].atom_name, " H3 " );
		break;
	}


	// THR T
	switch(opt)
	{
	case iupac:
		strcpy( AA[THR].atom[7].atom_name, " H  " );
		strcpy( AA[THR].atom[8].atom_name, " HG1" );
		strcpy( AA[THR].atom[9].atom_name, " HA " );
		strcpy( AA[THR].atom[10].atom_name," HB " );
		strcpy( AA[THR].atom[11].atom_name,"HG21" );
		strcpy( AA[THR].atom[12].atom_name,"HG22" );
		strcpy( AA[THR].atom[13].atom_name,"HG23" );
		strcpy( AA[THR].atom[15].atom_name," H2 " );
		strcpy( AA[THR].atom[16].atom_name," H3 " );
		break;
	case  pdb:
		strcpy( AA[THR].atom[7].atom_name, " H  " );
		strcpy( AA[THR].atom[8].atom_name, " HG1" );
		strcpy( AA[THR].atom[9].atom_name, " HA " );
		strcpy( AA[THR].atom[10].atom_name," HB " );
		strcpy( AA[THR].atom[11].atom_name,"HG21" );
		strcpy( AA[THR].atom[12].atom_name,"HG22" );
		strcpy( AA[THR].atom[13].atom_name,"HG23" );
		strcpy( AA[THR].atom[15].atom_name,"2H  " );
		strcpy( AA[THR].atom[16].atom_name,"3H  " );
		break;
	case sybil:
		strcpy( AA[THR].atom[7].atom_name, " H  " );
		strcpy( AA[THR].atom[8].atom_name, " HG1" );
		strcpy( AA[THR].atom[9].atom_name, " HA " );
		strcpy( AA[THR].atom[10].atom_name," HB " );
		strcpy( AA[THR].atom[11].atom_name,"HG21" );
		strcpy( AA[THR].atom[12].atom_name,"HG22" );
		strcpy( AA[THR].atom[13].atom_name,"HG23" );
		strcpy( AA[THR].atom[15].atom_name," H2 " );
		strcpy( AA[THR].atom[16].atom_name," H3 " );
		break;
	}

	// TRP W
	switch(opt)
	{
	case iupac:
		strcpy( AA[TRP].atom[14].atom_name, " H  " );
		strcpy( AA[TRP].atom[15].atom_name, " HE1" );
		strcpy( AA[TRP].atom[16].atom_name, " HD1" );
		strcpy( AA[TRP].atom[17].atom_name, " HZ2" );
		strcpy( AA[TRP].atom[18].atom_name, " HH2" );
		strcpy( AA[TRP].atom[19].atom_name, " HZ3" );
		strcpy( AA[TRP].atom[20].atom_name, " HE3" );
		strcpy( AA[TRP].atom[21].atom_name, " HA " );
		strcpy( AA[TRP].atom[22].atom_name, " HB2" );
		strcpy( AA[TRP].atom[23].atom_name, " HB3" );
		strcpy( AA[TRP].atom[25].atom_name, " H2 " );
		strcpy( AA[TRP].atom[26].atom_name, " H3 " );
		break;
	case  pdb:
		strcpy( AA[TRP].atom[14].atom_name, " H  " );
		strcpy( AA[TRP].atom[15].atom_name, " HE1" );
		strcpy( AA[TRP].atom[16].atom_name, " HD1" );
		strcpy( AA[TRP].atom[17].atom_name, " HZ2" );
		strcpy( AA[TRP].atom[18].atom_name, " HH2" );
		strcpy( AA[TRP].atom[19].atom_name, " HZ3" );
		strcpy( AA[TRP].atom[20].atom_name, " HE3" );
		strcpy( AA[TRP].atom[21].atom_name, " HA " );
		strcpy( AA[TRP].atom[22].atom_name, "1HB " );
		strcpy( AA[TRP].atom[23].atom_name, "2HB " );
		strcpy( AA[TRP].atom[25].atom_name, "2H  " );
		strcpy( AA[TRP].atom[26].atom_name, "3H  " );
		break;
	case sybil:
		strcpy( AA[TRP].atom[14].atom_name, " H  " );
		strcpy( AA[TRP].atom[15].atom_name, " HE1" );
		strcpy( AA[TRP].atom[16].atom_name, " HD1" );
		strcpy( AA[TRP].atom[17].atom_name, " HZ2" );
		strcpy( AA[TRP].atom[18].atom_name, " HH2" );
		strcpy( AA[TRP].atom[19].atom_name, " HZ3" );
		strcpy( AA[TRP].atom[20].atom_name, " HE3" );
		strcpy( AA[TRP].atom[21].atom_name, " HA " );
		strcpy( AA[TRP].atom[22].atom_name, " HB2" );
		strcpy( AA[TRP].atom[23].atom_name, " HB1" );
		strcpy( AA[TRP].atom[25].atom_name, " H2 " );
		strcpy( AA[TRP].atom[26].atom_name, " H3 " );
		break;
	}

	// TYR Y
	switch(opt)
	{
	case iupac:
		strcpy( AA[TYR].atom[12].atom_name, " H  " );
		strcpy( AA[TYR].atom[13].atom_name, " HH " );
		strcpy( AA[TYR].atom[14].atom_name, " HD1" );
		strcpy( AA[TYR].atom[15].atom_name, " HE1" );
		strcpy( AA[TYR].atom[16].atom_name, " HE2" );
		strcpy( AA[TYR].atom[17].atom_name, " HD2" );
		strcpy( AA[TYR].atom[18].atom_name, " HA " );
		strcpy( AA[TYR].atom[19].atom_name, " HB2" );
		strcpy( AA[TYR].atom[20].atom_name, " HB3" );
		strcpy( AA[TYR].atom[22].atom_name, " H2 " );
		strcpy( AA[TYR].atom[23].atom_name, " H3 " );
		break;
	case  pdb:
		strcpy( AA[TYR].atom[12].atom_name, " H  " );
		strcpy( AA[TYR].atom[13].atom_name, " HH " );
		strcpy( AA[TYR].atom[14].atom_name, " HD1" );
		strcpy( AA[TYR].atom[15].atom_name, " HE1" );
		strcpy( AA[TYR].atom[16].atom_name, " HE2" );
		strcpy( AA[TYR].atom[17].atom_name, " HD2" );
		strcpy( AA[TYR].atom[18].atom_name, " HA " );
		strcpy( AA[TYR].atom[19].atom_name, "1HB " );
		strcpy( AA[TYR].atom[20].atom_name, "2HB " );
		strcpy( AA[TYR].atom[22].atom_name, "2H  " );
		strcpy( AA[TYR].atom[23].atom_name, "3H  " );
		break;
	case sybil:
		strcpy( AA[TYR].atom[12].atom_name, " H  " );
		strcpy( AA[TYR].atom[13].atom_name, " HH " );
		strcpy( AA[TYR].atom[14].atom_name, " HD1" );
		strcpy( AA[TYR].atom[15].atom_name, " HE1" );
		strcpy( AA[TYR].atom[16].atom_name, " HE2" );
		strcpy( AA[TYR].atom[17].atom_name, " HD2" );
		strcpy( AA[TYR].atom[18].atom_name, " HA " );
		strcpy( AA[TYR].atom[19].atom_name, " HB2" );
		strcpy( AA[TYR].atom[20].atom_name, " HB1" );
		strcpy( AA[TYR].atom[22].atom_name, " H2 " );
		strcpy( AA[TYR].atom[23].atom_name, " H3 " );
		break;
	}

	// TYM Y Negative TYR
	switch(opt)
	{
	case iupac:
		strcpy( AA[TYM].atom[12].atom_name, " H  " );
		strcpy( AA[TYM].atom[13].atom_name, " HD1" );
		strcpy( AA[TYM].atom[14].atom_name, " HE1" );
		strcpy( AA[TYM].atom[15].atom_name, " HE2" );
		strcpy( AA[TYM].atom[16].atom_name, " HD2" );
		strcpy( AA[TYM].atom[17].atom_name, " HA " );
		strcpy( AA[TYM].atom[18].atom_name, " HB2" );
		strcpy( AA[TYM].atom[19].atom_name, " HB3" );
		strcpy( AA[TYM].atom[21].atom_name, " H2 " );
		strcpy( AA[TYM].atom[22].atom_name, " H3 " );
		break;
	case  pdb:
		strcpy( AA[TYM].atom[12].atom_name, " H  " );
		strcpy( AA[TYM].atom[13].atom_name, " HD1" );
		strcpy( AA[TYM].atom[14].atom_name, " HE1" );
		strcpy( AA[TYM].atom[15].atom_name, " HE2" );
		strcpy( AA[TYM].atom[16].atom_name, " HD2" );
		strcpy( AA[TYM].atom[17].atom_name, " HA " );
		strcpy( AA[TYM].atom[18].atom_name, "1HB " );
		strcpy( AA[TYM].atom[19].atom_name, "2HB " );
		strcpy( AA[TYM].atom[21].atom_name, "2H  " );
		strcpy( AA[TYM].atom[22].atom_name, "3H  " );
		break;
	case sybil:
		strcpy( AA[TYM].atom[12].atom_name, " H  " );
		strcpy( AA[TYM].atom[13].atom_name, " HD1" );
		strcpy( AA[TYM].atom[14].atom_name, " HE1" );
		strcpy( AA[TYM].atom[15].atom_name, " HE2" );
		strcpy( AA[TYM].atom[16].atom_name, " HD2" );
		strcpy( AA[TYM].atom[17].atom_name, " HA " );
		strcpy( AA[TYM].atom[18].atom_name, " HB2" );
		strcpy( AA[TYM].atom[19].atom_name, " HB1" );
		strcpy( AA[TYM].atom[21].atom_name, " H2 " );
		strcpy( AA[TYM].atom[22].atom_name, " H3 " );
		break;
	}

	// LYS K
	switch(opt)
	{
	case iupac:
		strcpy( AA[LYS].atom[9].atom_name,  " H  " );
		strcpy( AA[LYS].atom[10].atom_name, " HZ1" );
		strcpy( AA[LYS].atom[11].atom_name, " HZ2" );
		strcpy( AA[LYS].atom[12].atom_name, " HZ3" );
		strcpy( AA[LYS].atom[13].atom_name, " HA " );
		strcpy( AA[LYS].atom[14].atom_name, " HB2" );
		strcpy( AA[LYS].atom[15].atom_name, " HB3" );
		strcpy( AA[LYS].atom[16].atom_name, " HG2" );
		strcpy( AA[LYS].atom[17].atom_name, " HG3" );
		strcpy( AA[LYS].atom[18].atom_name, " HD2" );
		strcpy( AA[LYS].atom[19].atom_name, " HD3" );
		strcpy( AA[LYS].atom[20].atom_name, " HE2" );
		strcpy( AA[LYS].atom[21].atom_name, " HE3" );
		strcpy( AA[LYS].atom[22].atom_name, " OXT" );
		strcpy( AA[LYS].atom[23].atom_name, " H2 " );
		strcpy( AA[LYS].atom[24].atom_name, " H3 " );
		break;
	case  pdb:
		strcpy( AA[LYS].atom[9].atom_name,  " H  " );
		strcpy( AA[LYS].atom[10].atom_name, "1HZ " );
		strcpy( AA[LYS].atom[11].atom_name, "2HZ " );
		strcpy( AA[LYS].atom[12].atom_name, "3HZ " );
		strcpy( AA[LYS].atom[13].atom_name, " HA " );
		strcpy( AA[LYS].atom[14].atom_name, "1HB " );
		strcpy( AA[LYS].atom[15].atom_name, "2HB " );
		strcpy( AA[LYS].atom[16].atom_name, "1HG " );
		strcpy( AA[LYS].atom[17].atom_name, "2HG " );
		strcpy( AA[LYS].atom[18].atom_name, "1HD " );
		strcpy( AA[LYS].atom[19].atom_name, "2HD " );
		strcpy( AA[LYS].atom[20].atom_name, "1HE " );
		strcpy( AA[LYS].atom[21].atom_name, "2HE " );
		strcpy( AA[LYS].atom[22].atom_name, " OXT" );
		strcpy( AA[LYS].atom[23].atom_name, "2H  " );
		strcpy( AA[LYS].atom[24].atom_name, "3H  " );
		break;
	case sybil:
		strcpy( AA[LYS].atom[9].atom_name,  " H  " );
		strcpy( AA[LYS].atom[10].atom_name, " HZ1" );
		strcpy( AA[LYS].atom[11].atom_name, " HZ2" );
		strcpy( AA[LYS].atom[12].atom_name, " HZ3" );
		strcpy( AA[LYS].atom[13].atom_name, " HA " );
		strcpy( AA[LYS].atom[14].atom_name, " HB2" );
		strcpy( AA[LYS].atom[15].atom_name, " HB1" );
		strcpy( AA[LYS].atom[16].atom_name, " HG2" );
		strcpy( AA[LYS].atom[17].atom_name, " HG1" );
		strcpy( AA[LYS].atom[18].atom_name, " HD2" );
		strcpy( AA[LYS].atom[19].atom_name, " HD1" );
		strcpy( AA[LYS].atom[20].atom_name, " HE2" );
		strcpy( AA[LYS].atom[21].atom_name, " HE1" );
		strcpy( AA[LYS].atom[22].atom_name, " OXT" );
		strcpy( AA[LYS].atom[23].atom_name, " H2 " );
		strcpy( AA[LYS].atom[24].atom_name, " H3 " );
		break;
	}

	// LYN K Neutral LYS
	switch(opt)
	{
	case iupac:
		strcpy( AA[LYN].atom[9].atom_name,  " H  " );
		strcpy( AA[LYN].atom[10].atom_name, " HZ1" );
		strcpy( AA[LYN].atom[11].atom_name, " HZ2" );
		strcpy( AA[LYN].atom[12].atom_name, " HA " );
		strcpy( AA[LYN].atom[13].atom_name, " HB2" );
		strcpy( AA[LYN].atom[14].atom_name, " HB3" );
		strcpy( AA[LYN].atom[15].atom_name, " HG2" );
		strcpy( AA[LYN].atom[16].atom_name, " HG3" );
		strcpy( AA[LYN].atom[17].atom_name, " HD2" );
		strcpy( AA[LYN].atom[18].atom_name, " HD3" );
		strcpy( AA[LYN].atom[19].atom_name, " HE2" );
		strcpy( AA[LYN].atom[20].atom_name, " HE3" );
		strcpy( AA[LYN].atom[21].atom_name, " OXT" );
		strcpy( AA[LYN].atom[22].atom_name, " H2 " );
		strcpy( AA[LYN].atom[23].atom_name, " H3 " );
		break;
	case  pdb:
		strcpy( AA[LYN].atom[9].atom_name,  " H  " );
		strcpy( AA[LYN].atom[10].atom_name, "1HZ " );
		strcpy( AA[LYN].atom[11].atom_name, "2HZ " );
		strcpy( AA[LYN].atom[12].atom_name, " HA " );
		strcpy( AA[LYN].atom[13].atom_name, "1HB " );
		strcpy( AA[LYN].atom[14].atom_name, "2HB " );
		strcpy( AA[LYN].atom[15].atom_name, "1HG " );
		strcpy( AA[LYN].atom[16].atom_name, "2HG " );
		strcpy( AA[LYN].atom[17].atom_name, "1HD " );
		strcpy( AA[LYN].atom[18].atom_name, "2HD " );
		strcpy( AA[LYN].atom[19].atom_name, "1HE " );
		strcpy( AA[LYN].atom[20].atom_name, "2HE " );
		strcpy( AA[LYN].atom[21].atom_name, " OXT" );
		strcpy( AA[LYN].atom[22].atom_name, "2H  " );
		strcpy( AA[LYN].atom[23].atom_name, "3H  " );
		break;
	case sybil:
		strcpy( AA[LYN].atom[9].atom_name,  " H  " );
		strcpy( AA[LYN].atom[10].atom_name, " HZ1" );
		strcpy( AA[LYN].atom[11].atom_name, " HZ2" );
		strcpy( AA[LYN].atom[12].atom_name, " HA " );
		strcpy( AA[LYN].atom[13].atom_name, " HB2" );
		strcpy( AA[LYN].atom[14].atom_name, " HB1" );
		strcpy( AA[LYN].atom[15].atom_name, " HG2" );
		strcpy( AA[LYN].atom[16].atom_name, " HG1" );
		strcpy( AA[LYN].atom[17].atom_name, " HD2" );
		strcpy( AA[LYN].atom[18].atom_name, " HD1" );
		strcpy( AA[LYN].atom[19].atom_name, " HE2" );
		strcpy( AA[LYN].atom[20].atom_name, " HE1" );
		strcpy( AA[LYN].atom[21].atom_name, " OXT" );
		strcpy( AA[LYN].atom[22].atom_name, " H2 " );
		strcpy( AA[LYN].atom[23].atom_name, " H3 " );
		break;
	}

}

// number to convention name
Convention int2cff (int ff_convention)
{
	Convention conv_ff;

	switch(ff_convention)
		{
		case 0:
			// fprintf(stdout,"Force Field: CHARMM\n");
			conv_ff=Rosseta;
			break;
		case 1:
			// fprintf(stdout,"Force Field: ICM\n");
			conv_ff=ICM;
			break;
		case 2:
			// fprintf(stdout,"Force Field: EEF1\n");
			conv_ff=EEF1;
			break;
		case 3:
			// fprintf(stdout,"Force Field: Sybil\n");
			conv_ff=Sybil;
			break;
		default:
			fprintf(stdout,"Invalid convention set to CHARMM\n");
			conv_ff=Rosseta;
			break;
		}


	return conv_ff;
}

// number to convention name
ConventionNames int2cName (int name_convention) {

	ConventionNames conv;

	switch(name_convention)
		{
		case 0:
			// fprintf(stdout,"Name H convention: iupac\n");
			conv=iupac;
			break;
		case 1:
			// fprintf(stdout,"Name H convention: pdb\n");
			conv=pdb;
			break;
		case 2:
			// fprintf(stdout,"Name H convention: Sybil\n");
			conv=sybil;
			break;
		default:
			fprintf(stdout,"Invalid H name convention set to PDB 3.0\n");
			conv=pdb;
			break;
		}

	return conv;



}


void init_aminoacids(Convention opt, ConventionNames opt2)
{
if(!aminoacids_is_init)
{
 aminoacids_is_init=true;


 // FF convetion names...
 switch(opt)
 {

  case Rosseta:
    atom_types=atom_types_Rosseta;
    num_atom_type=54;
    break;

  case ICM:
    atom_types=atom_types_ICM;
    num_atom_type=54;
    break;

  case EEF1:
    atom_types=atom_types_EEF1;
    num_atom_type=31;
    break;
  case Sybil:
      atom_types=atom_types_Sybil;
      num_atom_type=45;
    break;

}


  init_ALA(opt);  init_CYS(opt);  init_ASP(opt); init_PHE(opt);  init_GLY(opt);  init_HIS(opt);
  init_ILE(opt);  init_GLU(opt);  init_LYS(opt); init_LEU(opt);  init_MET(opt);  init_ASN(opt);
  init_PRO(opt);  init_GLN(opt);  init_ARG(opt); init_SER(opt);  init_THR(opt);  init_VAL(opt);
  init_TRP(opt);  init_TYR(opt);  init_ASH(opt); init_CYX(opt);  init_CYM(opt);  init_GLH(opt);
  init_HIP(opt);  init_HID(opt);  init_HIE(opt); init_LYN(opt);	 init_TYM(opt);  init_MSE(opt);
  init_NtE(opt); init_CtE(opt);
  nucleotids();
  // MON: Do not alter this order. HIS, HIP, HID, and HIE must have ids 6, 24, 25, and 26, respectively (required for KORPE).

	for(int i=0; i<N_AMINO; i++)
	{
		AA[i].nelec = 0.0;
		for(int k=0; k<AA[i].nheavyatoms; k++)
		{
			int atn = atom_types[AA[i].atom[k].fullatom_type-1].at;
			AA[i].nelec += Table_Elements::getElement(atn)->number;
		}
	}

  // change if 1-2 (iupac) instead 2-3 H (pdb format)
//
//  if(opt2==pdb) {
////    aa_iupac2pdb();
//    fprintf(stdout,"->Name convention: pdb\n");
//  }
//
//  else if (opt2==iupac) {
//  //   aa_iupac2pqr();
//    fprintf(stdout,"->Name convention: iupac \n");
//  }
//
//  else if (opt2==sybil) {
//    //   aa_iupac2pqr();
//      fprintf(stdout,"->Name convention: sybil \n");
//    }


  change_H_names(opt2);


}
}


//Inicializacion de tabla de fuerzas de VDW para Carbono
void init_vdw_c(float probeRad,float probeEmax,float r, float f1, float f2)
{
    int i,j;
    float dist;
    float eps, sig, rat2, rat6;

    //FILE *f=fopen("tablaVDW.txt","wt");


    for(i=0;i<num_atom_type;i++)
    {
      //printf("atom_types=%f\n",atom_types[17].vdwd);
      //exit(1);
      //eps=sqrt(atom_types[i].vdwd*atom_types[18].vdwd);
      //eps=sqrt(0.14*0.14);
      eps=sqrt(atom_types[i].vdwd*0.14);
      //sig=atom_types[i].vdw+atom_types[18].vdw;

			sig=atom_types[i].vdw*r+probeRad;

      //fprintf(f,"N_AMINO=%3d:\t",i);
      for(j=0;j<MAX_VDW;j++)
      {
        dist=j*VDW_STEP;
        if(dist==0)
          VDW_C[i][j]=probeEmax;
        else
        {
          rat2=(sig*sig)/(dist*dist);
          rat6=rat2*rat2*rat2;
          VDW_C[i][j]=(eps*(f1*rat6*rat6-f2*2*rat6));
          if(VDW_C[i][j]>0)
            VDW_C[i][j]=(VDW_C[i][j]*probeEmax)/(VDW_C[i][j]+probeEmax);
        }

      }
    }
/*
    for(j=0;j<MAX_VDW;j++)
    {
      fprintf(f,"\n");
      for(i=0;i<N_AMINO;i++)
      {
        fprintf(f,"%10.6f\t",VDW_C[i][j]);
      }
    }


  fclose(f);*/
}




void init_vdw_h(void)
{
    int i,j;
    float dist;
    float eps, sig, rat2, rat6;

    for(i=0;i<num_atom_type;i++)
    {
      eps=sqrt(atom_types[i].vdwd*atom_types[21].vdwd);
      sig=atom_types[i].vdw+atom_types[21].vdw;
      for(j=0;j<MAX_VDW;j++)
      {
        dist=j*VDW_STEP;
        if(dist==0)
          VDW_H[i][j]=EMAX;
        else
        {
          rat2=(sig*sig)/(dist*dist);
          rat6=rat2*rat2*rat2;
          VDW_H[i][j]=(eps*(rat6*rat6-2*rat6));
          if(VDW_H[i][j]>0)
            VDW_H[i][j]=(VDW_H[i][j]*EMAX)/(VDW_H[i][j]+EMAX);
        }

      }
    }
}

void init_ele(void)
{
    int j;
    float dist;
    float beta=86.9525;//(epsilon0-A);
    float lambda_b=(-0.003627)*beta;
    float er;
    //float diel=-0.1146;
    float approx_zero=1.0E-6;
      for(j=0;j<MAX_ELE;j++)
      {
        dist=j*ELE_STEP;
        if(dist<0.5)
        {
          dist=0.5;
        }

        er=(-8.5525+(beta/(1.0+7.7839*exp(lambda_b*dist))));
        if(er<approx_zero)
          er=1.0;

        //ELE[j]=(332.0/(dist*dist*4.0));
        ELE[j]=(332.0*ELEC_MAX)/((dist*dist*4.0)+ELEC_MAX);
       // ELE[j]=(332.0/(dist*er));
      }
}

//version PABLO

void change_charges_withoutH()
{

// ALA

    AA[ALA].atom[0].charge = -0.16; //    N
    AA[ALA].atom[1].charge = 0.16; //    CA
    AA[ALA].atom[2].charge = 0.51; //    C
    AA[ALA].atom[3].charge = -0.51; //    O
    AA[ALA].atom[4].charge = 0.00; //    CB
    AA[ALA].atom[5].charge = 0.00; //    H
    AA[ALA].atom[6].charge = 0.00; //    HA
    AA[ALA].atom[7].charge = 0.00; //   1HB
    AA[ALA].atom[8].charge = 0.00; //   2HB
    AA[ALA].atom[9].charge = 0.00; //   3HB

    AA[ALA].atom[10].charge = -0.67; //   0XT
    AA[ALA].atom[11].charge = 0.00; //  1H
    AA[ALA].atom[12].charge = 0.00; //  2H

// CYS
    AA[CYS].atom[0].charge = -0.16; //   N
    AA[CYS].atom[1].charge = 0.16; //   CA
    AA[CYS].atom[2].charge = 0.51; //   C
    AA[CYS].atom[3].charge = -0.51; //   O
    AA[CYS].atom[4].charge =  0.07; //   CB
    AA[CYS].atom[5].charge = -0.07; //   SG
    AA[CYS].atom[6].charge = 0.00; //   H
    AA[CYS].atom[7].charge = 0.00; //   HA
    AA[CYS].atom[8].charge = 0.00; //   2HB
    AA[CYS].atom[9].charge = 0.00; //   3HB
    AA[CYS].atom[10].charge = 0.00; //   HG

    AA[CYS].atom[11].charge = -0.67; //   0XT
    AA[CYS].atom[12].charge = 0.00; //  1H
    AA[CYS].atom[13].charge = 0.00; //  2H

// ASP

     AA[ASP].atom[0].charge = -0.16; //    N
     AA[ASP].atom[1].charge = 0.16; //    CA
     AA[ASP].atom[2].charge = 0.51; //     C
     AA[ASP].atom[3].charge = -0.51; //    O
     AA[ASP].atom[4].charge = -0.10; //    CB
     AA[ASP].atom[5].charge = 0.6; //    CG
     AA[ASP].atom[6].charge = -0.76; //    OD1
     AA[ASP].atom[7].charge = -0.76; //    OD2
     AA[ASP].atom[8].charge = 0.00; //    H  <--- Nt(0.33)
     AA[ASP].atom[9].charge = 0.00; //    HA  <--- Nt(0.10)
     AA[ASP].atom[10].charge = 0.00; //    2HB
     AA[ASP].atom[11].charge = 0.00; //    3HB

     AA[ASP].atom[12].charge = -0.67; //   0XT
     AA[ASP].atom[13].charge = 0.00; //  1H
     AA[ASP].atom[14].charge = 0.00; //  2H

// GLU

      AA[GLU].atom[0].charge = -0.16; //    N
      AA[GLU].atom[1].charge = 0.16; //    CA
      AA[GLU].atom[2].charge = 0.51; //    C
      AA[GLU].atom[3].charge = -0.51; //    O
      AA[GLU].atom[4].charge =  0.00; //    CB
      AA[GLU].atom[5].charge = -0.18; //    CG
      AA[GLU].atom[6].charge = 0.62; //    CD
      AA[GLU].atom[7].charge = -0.76; //    OE1
      AA[GLU].atom[8].charge = -0.76; //    OE2
      AA[GLU].atom[9].charge = 0.00; //    H
      AA[GLU].atom[10].charge = 0.00; //    HA
      AA[GLU].atom[11].charge = 0.00; //    2HB
      AA[GLU].atom[12].charge = 0.00; //    3HB
      AA[GLU].atom[13].charge = 0.00; //    2HG
      AA[GLU].atom[14].charge = 0.00; //    3HG

      AA[GLU].atom[15].charge = -0.67; //   0XT
      AA[GLU].atom[16].charge = 0.00; //  1H
      AA[GLU].atom[17].charge = 0.00; //  2H

// PHE

      AA[PHE].atom[0].charge = -0.16; //    N
      AA[PHE].atom[1].charge = 0.16; //    CA  <--- Nt(0.21)
      AA[PHE].atom[2].charge = 0.51; //    C
      AA[PHE].atom[3].charge = -0.51; //    O
      AA[PHE].atom[4].charge = 0.00; //    CB
      AA[PHE].atom[5].charge = 0.00; //    CG
      AA[PHE].atom[6].charge = 0.00; //    CD1
      AA[PHE].atom[7].charge = 0.00; //    CD2
      AA[PHE].atom[8].charge = 0.00; //    CE1
      AA[PHE].atom[9].charge = 0.00; //    CE2
      AA[PHE].atom[10].charge = 0.00; //    CZ
      AA[PHE].atom[11].charge = 0.00; //    H
      AA[PHE].atom[12].charge = 0.00; //    HD1
      AA[PHE].atom[13].charge = 0.00; //    HE1
      AA[PHE].atom[14].charge = 0.00; //    HZ
      AA[PHE].atom[15].charge = 0.00; //    HE2
      AA[PHE].atom[16].charge = 0.00; //    HD2
      AA[PHE].atom[17].charge = 0.00; //    HA
      AA[PHE].atom[18].charge = 0.00; //    2HB
      AA[PHE].atom[19].charge = 0.00; //    3HB

      AA[PHE].atom[20].charge = -0.67; //   0XT
      AA[PHE].atom[21].charge = 0.00; //  1H
      AA[PHE].atom[22].charge = 0.00; //  2H

 // GLY

      AA[GLY].atom[0].charge = -0.16; //    N <---Nt(-0.30)
      AA[GLY].atom[1].charge =  0.16; //    CA <---Nt(0.13)
      AA[GLY].atom[2].charge = 0.51; //    C <---Ct(0.34)
      AA[GLY].atom[3].charge = -0.51; //    O <---Ct(0.67)
      AA[GLY].atom[4].charge = 0.00; //    H <---Nt(0.33)
      AA[GLY].atom[5].charge = 0.00; //    2HA
      AA[GLY].atom[6].charge = 0.00; //    3HA

// HYS

      AA[HIE].atom[0].charge = -0.16; //    N
      AA[HIE].atom[1].charge = 0.16; //    CA
      AA[HIE].atom[2].charge = 0.51; //    C
      AA[HIE].atom[3].charge = -0.51; //    O
      AA[HIE].atom[4].charge =  0.00; //    CB
      AA[HIE].atom[5].charge =  0.10; //    CG
      AA[HIE].atom[6].charge = -0.40; //    ND1
      AA[HIE].atom[7].charge =  0.10; //    CD2
      AA[HIE].atom[8].charge =  0.30; //    CE1
      AA[HIE].atom[9].charge = -0.10; //    NE2 --> HE1
      AA[HIE].atom[10].charge = 0.00; //    H
      AA[HIE].atom[11].charge = 0.00; //    HE2
      AA[HIE].atom[12].charge = 0.00; //    HA
      AA[HIE].atom[13].charge = 0.00; //    2HB
      AA[HIE].atom[14].charge = 0.00; //    3HB
      AA[HIE].atom[15].charge = 0.00; //    HE1 --->AQUI
      AA[HIE].atom[16].charge = 0.00; //    HD2 --- (no existe...en esto)

      AA[HIE].atom[17].charge = -0.67; //   0XT
      AA[HIE].atom[18].charge = 0.00; //  1H
      AA[HIE].atom[19].charge = 0.00; //  2H


// ILE

      AA[ILE].atom[0].charge = -0.16; //    N
      AA[ILE].atom[1].charge = 0.16; //    CA
      AA[ILE].atom[2].charge = 0.51; //    C
      AA[ILE].atom[3].charge = -0.51; //    O
      AA[ILE].atom[4].charge =  0.00; //    CB
      AA[ILE].atom[5].charge =  0.00; //    CG1
      AA[ILE].atom[6].charge =  0.00; //    CG2
      AA[ILE].atom[7].charge =  0.00; //    CD1
      AA[ILE].atom[8].charge = 0.00; //    H
      AA[ILE].atom[9].charge = 0.00; //    HA
      AA[ILE].atom[10].charge = 0.00; //    HB
      AA[ILE].atom[11].charge = 0.00; //    1HG2
      AA[ILE].atom[12].charge = 0.00; //    2HG2
      AA[ILE].atom[13].charge = 0.00; //    3HG2
      AA[ILE].atom[14].charge = 0.00; //    2HG1
      AA[ILE].atom[15].charge = 0.00; //    3HG1
      AA[ILE].atom[16].charge = 0.00; //    1HD1
      AA[ILE].atom[17].charge = 0.00; //    2HD1
      AA[ILE].atom[18].charge = 0.00; //    3HD1

      AA[ILE].atom[19].charge = -0.67; //   0XT
      AA[ILE].atom[20].charge = 0.00; //  1H
      AA[ILE].atom[21].charge = 0.00; //  2H

// LYS
			//correcto
      AA[LYS].atom[0].charge = -0.16; //    N
      AA[LYS].atom[1].charge = 0.16; //    CA
      AA[LYS].atom[2].charge = 0.51; //    C
      AA[LYS].atom[3].charge = -0.51; //    O
      AA[LYS].atom[4].charge =  0.00; //    CB
      AA[LYS].atom[5].charge =  0.00; //    CG
      AA[LYS].atom[6].charge =  0.00; //    CD
      AA[LYS].atom[7].charge =  0.31; //    CE
      AA[LYS].atom[8].charge =  0.69; //    NZ
      AA[LYS].atom[9].charge =  0.00; //    H
      AA[LYS].atom[10].charge = 0.00; //    1HZ
      AA[LYS].atom[11].charge = 0.00; //    2HZ
      AA[LYS].atom[12].charge = 0.00; //    3HZ
      AA[LYS].atom[13].charge = 0.00; //    HA
      AA[LYS].atom[14].charge = 0.00; //    2HB
      AA[LYS].atom[15].charge = 0.00; //    3HB
      AA[LYS].atom[16].charge = 0.00; //    2HG
      AA[LYS].atom[17].charge = 0.00; //    3HG
      AA[LYS].atom[18].charge = 0.00; //    2HD
      AA[LYS].atom[19].charge = 0.00; //    3HD
      AA[LYS].atom[20].charge = 0.00; //    2HE
      AA[LYS].atom[21].charge = 0.00; //    3HE

//  LEU

      AA[LEU].atom[0].charge = -0.16; //    N
      AA[LEU].atom[1].charge =  0.16; //    CA
      AA[LEU].atom[2].charge =  0.51; //    C
      AA[LEU].atom[3].charge = -0.51; //    O
      AA[LEU].atom[4].charge =  0.00; //    CB
      AA[LEU].atom[5].charge =  0.00; //    CG
      AA[LEU].atom[6].charge =  0.00; //    CD1
      AA[LEU].atom[7].charge =  0.00; //    CD2
      AA[LEU].atom[8].charge = 0.00; //    H
      AA[LEU].atom[9].charge = 0.00; //    HA
      AA[LEU].atom[10].charge = 0.00; //    2HB
      AA[LEU].atom[11].charge = 0.00; //    3HB
      AA[LEU].atom[12].charge = 0.00; //    HG
      AA[LEU].atom[13].charge = 0.00; //    1HD1
      AA[LEU].atom[14].charge = 0.00; //    2HD1
      AA[LEU].atom[15].charge = 0.00; //    3HD1
      AA[LEU].atom[16].charge = 0.00; //    1HD2
      AA[LEU].atom[17].charge = 0.00; //    2HD2
      AA[LEU].atom[18].charge = 0.00; //    3HD2

      AA[LEU].atom[19].charge = -0.67; //   0XT
      AA[LEU].atom[20].charge = 0.00; //  1H
      AA[LEU].atom[21].charge = 0.00; //  2H

 //  MET

      AA[MET].atom[0].charge = -0.16; //    N
      AA[MET].atom[1].charge = 0.16; //    CA
      AA[MET].atom[2].charge = 0.51; //    C
      AA[MET].atom[3].charge = -0.51; //    O
      AA[MET].atom[4].charge =  0.00; //    CB
      AA[MET].atom[5].charge =  0.04; //    CG
      AA[MET].atom[6].charge = -0.09; //    SD
      AA[MET].atom[7].charge =  0.05; //    CE
      AA[MET].atom[8].charge = 0.00; //    H
      AA[MET].atom[9].charge = 0.00; //    HA
      AA[MET].atom[10].charge = 0.00; //   2HB
      AA[MET].atom[11].charge = 0.00; //   3HB
      AA[MET].atom[12].charge = 0.00; //   2HG
      AA[MET].atom[13].charge = 0.00; //   3HG
      AA[MET].atom[14].charge = 0.00; //   1HE
      AA[MET].atom[15].charge = 0.00; //   2HE
      AA[MET].atom[16].charge = 0.00; //   3HE

      AA[MET].atom[17].charge = -0.67; //   0XT
      AA[MET].atom[18].charge = 0.00; //  1H
      AA[MET].atom[19].charge = 0.00; //  2H



//   ASN

      AA[ASN].atom[0].charge = -0.16; //    N
      AA[ASN].atom[1].charge = 0.16; //    CA
      AA[ASN].atom[2].charge = 0.51; //    C
      AA[ASN].atom[3].charge = -0.51; //    O
      AA[ASN].atom[4].charge =  0.00; //    CB
      AA[ASN].atom[5].charge = 0.55; //    CG
      AA[ASN].atom[6].charge = -0.55; //    OD1
      AA[ASN].atom[7].charge =  0.00; //    ND2
      AA[ASN].atom[8].charge = 0.00; //    H
      AA[ASN].atom[9].charge = 0.00; //    1HD2
      AA[ASN].atom[10].charge = 0.00; //    2HD2
      AA[ASN].atom[11].charge = 0.00; //    HA
      AA[ASN].atom[12].charge = 0.00; //    2HB
      AA[ASN].atom[13].charge = 0.00; //    3HB

      AA[ASN].atom[14].charge = -0.67; //   0XT
      AA[ASN].atom[15].charge = 0.00; //  1H
      AA[ASN].atom[16].charge = 0.00; //  2H

// PRO

      AA[PRO].atom[0].charge = -0.29; //    N
      AA[PRO].atom[1].charge = 0.11; //    CA
      AA[PRO].atom[2].charge = 0.51; //    C
      AA[PRO].atom[3].charge = -0.51; //    O
      AA[PRO].atom[4].charge = 0.00; //    CB
      AA[PRO].atom[5].charge = 0.00; //    CG
      AA[PRO].atom[6].charge = 0.00; //    CD
      AA[PRO].atom[7].charge = 0.00; //   2HD
      AA[PRO].atom[8].charge = 0.00; //   3HD
      AA[PRO].atom[9].charge = 0.00; //   2HG
      AA[PRO].atom[10].charge = 0.00; //   3HG
      AA[PRO].atom[11].charge = 0.00; //   2HB
      AA[PRO].atom[12].charge = 0.00; //   3HB
      AA[PRO].atom[13].charge = 0.00; //    HA

      AA[PRO].atom[14].charge = -0.67; //   0XT
      AA[PRO].atom[15].charge = 0.00; //  1H
      AA[PRO].atom[16].charge = 0.00; //  2H

// GLN

      AA[GLN].atom[0].charge = -0.16; //    N
      AA[GLN].atom[1].charge = 0.16; //    CA
      AA[GLN].atom[2].charge = 0.51; //    C
      AA[GLN].atom[3].charge = -0.51; //    O
      AA[GLN].atom[4].charge = 0.00; //    CB
      AA[GLN].atom[5].charge = 0.00; //    CG
      AA[GLN].atom[6].charge = 0.55; //    CD
      AA[GLN].atom[7].charge = -0.55; //    OE1
      AA[GLN].atom[8].charge =  0.00; //    NE2
      AA[GLN].atom[9].charge = 0.00; //    H
      AA[GLN].atom[10].charge = 0.00; //   1HE2
      AA[GLN].atom[11].charge = 0.00; //   2HE2
      AA[GLN].atom[12].charge = 0.00; //    HA
      AA[GLN].atom[13].charge = 0.00; //   2HB
      AA[GLN].atom[14].charge = 0.00; //   3HB
      AA[GLN].atom[15].charge = 0.00; //   2HG
      AA[GLN].atom[16].charge = 0.00; //   3HG

      AA[GLN].atom[17].charge = -0.67; //   0XT
      AA[GLN].atom[18].charge = 0.00; //  1H
      AA[GLN].atom[19].charge = 0.00; //  2H

// ARG

      AA[ARG].atom[0].charge = -0.16; //    N
      AA[ARG].atom[1].charge = 0.16; //    CA
      AA[ARG].atom[2].charge = 0.51; //    C
      AA[ARG].atom[3].charge = -0.51; //    O
      AA[ARG].atom[4].charge = 0.00; //    CB
      AA[ARG].atom[5].charge = 0.00; //    CG
      AA[ARG].atom[6].charge = 0.05; //    CD
      AA[ARG].atom[7].charge = 0.10; //    NE
      AA[ARG].atom[8].charge = 0.30; //    CZ
      AA[ARG].atom[9].charge =  0.25; //    NH1
      AA[ARG].atom[10].charge = 0.25; //    NH2
      AA[ARG].atom[11].charge = 0.00; //    H
      AA[ARG].atom[12].charge = 0.00; //   1HH1
      AA[ARG].atom[13].charge = 0.00; //   2HH1
      AA[ARG].atom[14].charge = 0.00; //   1HH2
      AA[ARG].atom[15].charge = 0.00; //   2HH2
      AA[ARG].atom[16].charge = 0.00; //    HE
      AA[ARG].atom[17].charge = 0.00; //    HA
      AA[ARG].atom[18].charge = 0.00; //   2HB
      AA[ARG].atom[19].charge = 0.00; //   3HB
      AA[ARG].atom[20].charge = 0.00; //   2HG
      AA[ARG].atom[21].charge = 0.00; //   3HG
      AA[ARG].atom[22].charge = 0.00; //   2HD
      AA[ARG].atom[23].charge = 0.00; //   3HD

// SER

      AA[SER].atom[0].charge = -0.16; //    N
      AA[SER].atom[1].charge = 0.16; //    CA
      AA[SER].atom[2].charge = 0.51; //    C
      AA[SER].atom[3].charge = -0.51; //    O
      AA[SER].atom[4].charge = 0.13; //    CB
      AA[SER].atom[5].charge = -0.13; //    OG
      AA[SER].atom[6].charge = 0.00; //    H
      AA[SER].atom[7].charge = 0.00; //    HG
      AA[SER].atom[8].charge = 0.00; //    HA
      AA[SER].atom[9].charge = 0.00; //   2HB
      AA[SER].atom[10].charge = 0.00; //   3HB

      AA[SER].atom[11].charge = -0.67; //   0XT
      AA[SER].atom[12].charge = 0.00; //  1H
      AA[SER].atom[13].charge = 0.00; //  2H

// THR

      AA[THR].atom[0].charge = -0.16; //    N
      AA[THR].atom[1].charge = 0.16; //    CA
      AA[THR].atom[2].charge = 0.51; //    C
      AA[THR].atom[3].charge = -0.51; //    O
      AA[THR].atom[4].charge = 0.13; //    CB
      AA[THR].atom[5].charge = -0.13; //    OG1
      AA[THR].atom[6].charge = -0.00; //    CG2
      AA[THR].atom[7].charge = 0.00; //    H
      AA[THR].atom[8].charge = 0.00; //    HG1
      AA[THR].atom[9].charge = 0.00; //    HA
      AA[THR].atom[10].charge = 0.00; //    HB
      AA[THR].atom[11].charge = 0.00; //    1HG2
      AA[THR].atom[12].charge = 0.00; //    2HG2
      AA[THR].atom[13].charge = 0.00; //    3HG2

      AA[THR].atom[14].charge = -0.67; //   0XT
      AA[THR].atom[15].charge = 0.00; //  1H
      AA[THR].atom[16].charge = 0.00; //  2H

 // VAL

      AA[VAL].atom[0].charge = -0.16; //   N
      AA[VAL].atom[1].charge = 0.16; //   CA
      AA[VAL].atom[2].charge = 0.51; //   C
      AA[VAL].atom[3].charge = -0.51; //   O
      AA[VAL].atom[4].charge = -0.00; //   CB
      AA[VAL].atom[5].charge = -0.00; //   CG1
      AA[VAL].atom[6].charge = -0.00; //   CG2
      AA[VAL].atom[7].charge = 0.00; //   H
      AA[VAL].atom[8].charge = 0.00; //   HA
      AA[VAL].atom[9].charge = 0.00; //   HB
      AA[VAL].atom[10].charge = 0.00; //   1HG1
      AA[VAL].atom[11].charge = 0.00; //   2HG1
      AA[VAL].atom[12].charge = 0.00; //   3HG1
      AA[VAL].atom[13].charge = 0.00; //   1HG2
      AA[VAL].atom[14].charge = 0.00; //   2HG2
      AA[VAL].atom[15].charge = 0.00; //   3HG2

      AA[VAL].atom[16].charge = -0.67; //   0XT
      AA[VAL].atom[17].charge = 0.00; //  1H
      AA[VAL].atom[18].charge = 0.00; //  2H

// TRP

      AA[TRP].atom[0].charge = -0.16; //    N
      AA[TRP].atom[1].charge = 0.16; //    CA
      AA[TRP].atom[2].charge = 0.51; //    C
      AA[TRP].atom[3].charge = -0.51; //    O
      AA[TRP].atom[4].charge = -0.00; //    CB
      AA[TRP].atom[5].charge = -0.03; //    CG
      AA[TRP].atom[6].charge =  0.03; //    CD1
      AA[TRP].atom[7].charge =  0.10; //    CD2
      AA[TRP].atom[8].charge = -0.23; //    NE1
      AA[TRP].atom[9].charge =  0.13; //    CE2
      AA[TRP].atom[10].charge = 0.0; //    CE3
      AA[TRP].atom[11].charge = 0.0; //    CZ2
      AA[TRP].atom[12].charge = 0.0; //    CZ3
      AA[TRP].atom[13].charge = 0.0; //    CH2
      AA[TRP].atom[14].charge = 0.00; //    H
      AA[TRP].atom[15].charge = 0.00; //    HE1
      AA[TRP].atom[16].charge = 0.0; //    HD1
      AA[TRP].atom[17].charge = 0.0; //    HZ2
      AA[TRP].atom[18].charge = 0.0; //    HH2
      AA[TRP].atom[19].charge = 0.0; //    HZ3
      AA[TRP].atom[20].charge = 0.0; //    HE3
      AA[TRP].atom[21].charge = 0.00; //    HA
      AA[TRP].atom[22].charge = 0.00; //    2HB
      AA[TRP].atom[23].charge = 0.00; //    3HB

      AA[TRP].atom[24].charge = -0.67; //   0XT
      AA[TRP].atom[25].charge = 0.00; //  1H
      AA[TRP].atom[26].charge = 0.00; //  2H

// TYR

      AA[TYR].atom[0].charge = -0.16; //    N
      AA[TYR].atom[1].charge = 0.16; //    CA
      AA[TYR].atom[2].charge = 0.51; //    C
      AA[TYR].atom[3].charge = -0.51; //    O
      AA[TYR].atom[4].charge = 0.0; //    CB
      AA[TYR].atom[5].charge = 0.00; //    CG
      AA[TYR].atom[6].charge = 0.00; //    CD1
      AA[TYR].atom[7].charge = 0.00; //    CD2
      AA[TYR].atom[8].charge = 0.00; //    CE1
      AA[TYR].atom[9].charge = 0.00; //    CE2
      AA[TYR].atom[10].charge = 0.11; //    CZ
      AA[TYR].atom[11].charge = -0.11; //    OH
      AA[TYR].atom[12].charge = 0.00; //    H
      AA[TYR].atom[13].charge = 0.00; //    HH
      AA[TYR].atom[14].charge = 0.00; //    HD1
      AA[TYR].atom[15].charge = 0.00; //    HE1
      AA[TYR].atom[16].charge = 0.00; //    HE2
      AA[TYR].atom[17].charge = 0.00; //    HD2
      AA[TYR].atom[18].charge = 0.00; //    HA
      AA[TYR].atom[19].charge = 0.00; //   2HB
      AA[TYR].atom[20].charge = 0.00; //   3HB

      AA[TYR].atom[21].charge = -0.67; //   0XT
      AA[TYR].atom[22].charge = 0.00; //  1H
      AA[TYR].atom[23].charge = 0.00; //  2H

// ASH --- CHARMM ASPP

      AA[ASH].atom[0].charge = -0.16; //    N  <--- Nt(-0.30)
      AA[ASH].atom[1].charge = 0.16; //    CA  <--- Nt(0.21)
      AA[ASH].atom[2].charge = 0.51; //    C  <--- Ct(0.34)
      AA[ASH].atom[3].charge = -0.51; //    O  <--- Ct(-0.67)
      AA[ASH].atom[4].charge = 0.00; //    CB
      AA[ASH].atom[5].charge = 0.70; //    CG
      AA[ASH].atom[6].charge = -0.50; //    OD1
      AA[ASH].atom[7].charge = -0.20; //    OD2
      AA[ASH].atom[8].charge = 0.00; //     H
      AA[ASH].atom[9].charge = 0.00; //     HA
      AA[ASH].atom[10].charge = 0.00; //    2HB
      AA[ASH].atom[11].charge = 0.00; //    3HB
      AA[ASH].atom[12].charge = 0.00; //    2HD NEW ADDED - like CHARMM ASPP

      AA[ASH].atom[13].charge = -0.67; //   0XT
      AA[ASH].atom[14].charge = 0.00; //  1H
      AA[ASH].atom[15].charge = 0.00; //  2H

// CYX

    AA[CYX].atom[0].charge = -0.16; //   N  <--- Nt(-0.30)
    AA[CYX].atom[1].charge = 0.16; //   CA  <--- Nt(0.21)
    AA[CYX].atom[2].charge = 0.51; //   C  <--- Ct(0.34)
    AA[CYX].atom[3].charge = -0.51; //   O  <--- Ct(-0.67)
    AA[CYX].atom[4].charge =  0.08; //   CB
    AA[CYX].atom[5].charge = -0.08; //   SG
    AA[CYX].atom[6].charge = 0.00; //   H  <--- Nt(0.33)
    AA[CYX].atom[7].charge = 0.00; //   HA  <--- Nt(0.10)
    AA[CYX].atom[8].charge = 0.00; //   2HB
    AA[CYX].atom[9].charge = 0.00; //   3HB

    AA[CYX].atom[10].charge = -0.67; //   0XT
    AA[CYX].atom[11].charge = 0.00; //  1H
    AA[CYX].atom[12].charge = 0.00; //  2H

// GLH

      AA[GLH].atom[0].charge = -0.16; //    N
      AA[GLH].atom[1].charge = 0.16; //    CA
      AA[GLH].atom[2].charge = 0.51; //    C
      AA[GLH].atom[3].charge = -0.51; //    O
      AA[GLH].atom[4].charge = 0.00; //    CB
      AA[GLH].atom[5].charge = 0.0; //    CG
      AA[GLH].atom[6].charge = 0.70; //    CD
      AA[GLH].atom[7].charge = -0.50; //    OE1
      AA[GLH].atom[8].charge = -0.20; //    OE2
      AA[GLH].atom[9].charge = 0.00; //    H
      AA[GLH].atom[10].charge = 0.00; //    HA
      AA[GLH].atom[11].charge = 0.00; //    2HB
      AA[GLH].atom[12].charge = 0.00; //    3HB
      AA[GLH].atom[13].charge = 0.00; //    2HG
      AA[GLH].atom[14].charge = 0.00; //    3HG
      AA[GLH].atom[15].charge = 0.00; //    2HE NEW ADDED

      AA[GLH].atom[16].charge = -0.67; //   0XT
      AA[GLH].atom[17].charge = 0.00; //  1H
      AA[GLH].atom[18].charge = 0.00; //  2H

// HIP   HIP = HSC (Di-Protonated at ND1 & NE2, +)


      AA[HIP].atom[0].charge = -0.16; //    N
      AA[HIP].atom[1].charge = 0.16; //    CA
      AA[HIP].atom[2].charge = 0.51; //    C
      AA[HIP].atom[3].charge = -0.51; //    O
      AA[HIP].atom[4].charge = -0.00; //    CB
      AA[HIP].atom[5].charge =  0.00; //    CG
      AA[HIP].atom[6].charge =  0.40; //    ND1
      AA[HIP].atom[7].charge =  0.00; //    CD2
      AA[HIP].atom[8].charge =  0.00; //     CE1
      AA[HIP].atom[9].charge =  0.40; //    NE2
      AA[HIP].atom[10].charge = 0.00; //    H
      AA[HIP].atom[11].charge = 0.00; //    HE2
      AA[HIP].atom[12].charge = 0.00; //    HD1 //NEW ADDED
      AA[HIP].atom[13].charge = 0.00; //    HA
      AA[HIP].atom[14].charge = 0.00; //    2HB
      AA[HIP].atom[15].charge = 0.00; //    3HB
      AA[HIP].atom[16].charge = 0.00; //    HE1
      AA[HIP].atom[17].charge = 0.00; //    HD2

      AA[HIP].atom[18].charge = -0.67; //   0XT
      AA[HIP].atom[19].charge = 0.00; //  1H
      AA[HIP].atom[20].charge = 0.00; //  2H



// HID  (Protonated at ND1, neutral)

      AA[HIE].atom[0].charge = -0.16; //    N
      AA[HIE].atom[1].charge = 0.16; //    CA
      AA[HIE].atom[2].charge = 0.51; //    C
      AA[HIE].atom[3].charge = -0.51; //    O
      AA[HIE].atom[4].charge =  0.00; //    CB
      AA[HIE].atom[5].charge =  0.10; //    CG
      AA[HIE].atom[6].charge = -0.10; //    ND1 ---> HD1
      AA[HIE].atom[7].charge =  0.10; //    CD2
      AA[HIE].atom[8].charge =  0.30; //    CE1
      AA[HIE].atom[9].charge = -0.40; //    NE2
      AA[HIE].atom[10].charge = 0.00; //    H
      AA[HIE].atom[11].charge = 0.00; //    HD1 ---> aqui
      AA[HIE].atom[12].charge = 0.00; //    HA
      AA[HIE].atom[13].charge = 0.00; //    2HB
      AA[HIE].atom[14].charge = 0.00; //    3HB
      AA[HIE].atom[15].charge = 0.00; //    HE1
      AA[HIE].atom[16].charge = 0.00; //    HD2

      AA[HIE].atom[17].charge = -0.67; //   0XT
      AA[HIE].atom[18].charge = 0.00; //  1H
      AA[HIE].atom[19].charge = 0.00; //  2H



// HIE = HIS nuestra = HSD (Protonated at NE2, neutral)

      AA[HIE].atom[0].charge = -0.16; //    N
      AA[HIE].atom[1].charge = 0.16; //    CA
      AA[HIE].atom[2].charge = 0.51; //    C
      AA[HIE].atom[3].charge = -0.51; //    O
      AA[HIE].atom[4].charge =  0.00; //    CB
      AA[HIE].atom[5].charge =  0.10; //    CG
      AA[HIE].atom[6].charge = -0.40; //    ND1
      AA[HIE].atom[7].charge =  0.10; //    CD2
      AA[HIE].atom[8].charge =  0.30; //    CE1
      AA[HIE].atom[9].charge = -0.10; //    NE2 --> HE1
      AA[HIE].atom[10].charge = 0.00; //    H
      AA[HIE].atom[11].charge = 0.00; //    HE2
      AA[HIE].atom[12].charge = 0.00; //    HA
      AA[HIE].atom[13].charge = 0.00; //    2HB
      AA[HIE].atom[14].charge = 0.00; //    3HB
      AA[HIE].atom[15].charge = 0.00; //    HE1 --->AQUI
      AA[HIE].atom[16].charge = 0.00; //    HD2

      AA[HIE].atom[17].charge = -0.67; //   0XT
      AA[HIE].atom[18].charge = 0.00; //  1H
      AA[HIE].atom[19].charge = 0.00; //  2H

//
// TYM



      AA[TYM].atom[0].charge = -0.16; //    N
      AA[TYM].atom[1].charge = 0.16; //    CA
      AA[TYM].atom[2].charge = 0.51; //    C
      AA[TYM].atom[3].charge = -0.51; //    O
      AA[TYM].atom[4].charge = 0.0; //    CB
      AA[TYM].atom[5].charge = 0.00; //    CG
      AA[TYM].atom[6].charge = 0.00; //    CD1
      AA[TYM].atom[7].charge = 0.00; //    CD2
      AA[TYM].atom[8].charge = 0.00; //    CE1
      AA[TYM].atom[9].charge = 0.00; //    CE2
      AA[TYM].atom[10].charge = 0.11; //    CZ
      AA[TYM].atom[11].charge = -0.11; //    OH
      AA[TYM].atom[12].charge = 0.00; //    H
      AA[TYM].atom[13].charge = 0.00; //    HD1
      AA[TYM].atom[14].charge = 0.00; //    HE1
      AA[TYM].atom[15].charge = 0.00; //    HE2
      AA[TYM].atom[16].charge = 0.00; //    HD2
      AA[TYM].atom[17].charge = 0.00; //    HA
      AA[TYM].atom[18].charge = 0.00; //   2HB
      AA[TYM].atom[19].charge = 0.00; //   3HB

      AA[TYM].atom[20].charge = -0.67; //   0XT
      AA[TYM].atom[21].charge = 0.00; //  1H
      AA[TYM].atom[22].charge = 0.00; //  2H

// LYN

			//correcto
      AA[LYN].atom[0].charge = -0.16; //    N
      AA[LYN].atom[1].charge = 0.16; //    CA
      AA[LYN].atom[2].charge = 0.51; //    C
      AA[LYN].atom[3].charge = -0.51; //    O
      AA[LYN].atom[4].charge =  0.00; //    CB
      AA[LYN].atom[5].charge =  0.00; //    CG
      AA[LYN].atom[6].charge =  0.00; //    CD
      AA[LYN].atom[7].charge =  0.31; //    CE
      AA[LYN].atom[8].charge =  0.69; //    NZ
      AA[LYN].atom[9].charge =  0.00; //    H
      AA[LYN].atom[10].charge = 0.00; //    1HZ
      AA[LYN].atom[11].charge = 0.00; //    2HZ
      AA[LYN].atom[12].charge = 0.00; //    HA
      AA[LYN].atom[13].charge = 0.00; //    2HB
      AA[LYN].atom[14].charge = 0.00; //    3HB
      AA[LYN].atom[15].charge = 0.00; //    2HG
      AA[LYN].atom[16].charge = 0.00; //    3HG
      AA[LYN].atom[17].charge = 0.00; //    2HD
      AA[LYN].atom[18].charge = 0.00; //    3HD
      AA[LYN].atom[19].charge = 0.00; //    2HE
      AA[LYN].atom[20].charge = 0.00; //    3HE

//  MSE

      AA[MSE].atom[0].charge = -0.16; //    N
      AA[MSE].atom[1].charge = 0.16; //    CA
      AA[MSE].atom[2].charge = 0.51; //    C
      AA[MSE].atom[3].charge = -0.51; //    O
      AA[MSE].atom[4].charge =  0.00; //    CB
      AA[MSE].atom[5].charge =  0.04; //    CG
      AA[MSE].atom[6].charge = -0.09; //    SD
      AA[MSE].atom[7].charge =  0.05; //    CE
      AA[MSE].atom[8].charge = 0.00; //    H
      AA[MSE].atom[9].charge = 0.00; //    HA
      AA[MSE].atom[10].charge = 0.00; //   2HB
      AA[MSE].atom[11].charge = 0.00; //   3HB
      AA[MSE].atom[12].charge = 0.00; //   2HG
      AA[MSE].atom[13].charge = 0.00; //   3HG
      AA[MSE].atom[14].charge = 0.00; //   1HE
      AA[MSE].atom[15].charge = 0.00; //   2HE
      AA[MSE].atom[16].charge = 0.00; //   3HE

      AA[MSE].atom[17].charge = -0.67; //   0XT
      AA[MSE].atom[18].charge = 0.00; //  1H
      AA[MSE].atom[19].charge = 0.00; //  2H


}
