#define N_ATOM_MAX 54
#define N_AMINO 40
#define MAX_VDW 300 //number of steps in the vdw computation
#define VDW_STEP 0.1 //Stepsize in amstrongs
//#define EMAX 1.5 //Soft Limit of the vdw values
#define EMAX 1.2 //Soft Limit of the vdw values
#define ELEC_MAX 1.5
#define MAX_ELE 300
#define ELE_STEP 0.1

// Unique Residue (aminoacid and nucleotids) identifiers
short int const ALA =  0 ;
short int const CYS =  1 ;
short int const ASP =  2 ;
short int const GLU =  3 ;
short int const PHE =  4 ;
short int const GLY =  5 ;
short int const HIS =  6 ;
short int const ILE =  7 ;
short int const LYS =  8 ;
short int const LEU =  9 ;
short int const MET =  10 ;
short int const ASN =  11 ;
short int const PRO =  12 ;
short int const GLN =  13 ;
short int const ARG =  14 ;
short int const SER =  15 ;
short int const THR =  16 ;
short int const VAL =  17 ;
short int const TRP =  18 ;
short int const TYR =  19 ;
short int const ASH =  20 ;       //Neutral ASP
short int const CYX =  21 ;       //SS-bonded CYS
short int const CYM =  22 ;       //Negative CYS
short int const GLH =  23 ;       //Neutral GLU
short int const HIP =  24 ;       //Positive HIS
short int const HID =  25 ;       //Neutral HIS, proton HD1 present
short int const HIE =  26 ;       //Neutral HIS, proton HE2 present
short int const LYN =  27 ;       //Neutral LYS
short int const TYM =  28 ;       //Negative TYR
short int const MSE =  29 ; 		// Seleno-Methionine
//
short int const NtE =  30 ;
short int const CtE =  31 ;

// Mon modified (2/12/2009)
short int const DGUA =  32 ;
short int const DADE =  33 ;
short int const DCYT =  34 ;
short int const DTHY =  35 ;
short int const GUA =  36 ;
short int const ADE =  37 ;
short int const CYT =  38 ;
short int const URA =  39 ;

#ifndef Conventions
#define Conventions
/**
 * Convention of residues/atoms related values (vdw, charge, etc.)
 */
enum Convention {Rosseta, ICM, EEF1,Sybil};
/**
 * Name Convention of residues/atoms
 */
enum ConventionNames{iupac,pdb,sybil};
#endif

typedef int Tpair[2];
typedef char entero4[4];

//find the matirx (mat) and vector (vec) that rotate chi deg about an axis
//define by a1-->a2 rotation is right-handed. That is, it is a clockwise
//rotation when looking from a1 to a2 chi in radians
void getrotS_bk( float * a1, float * a2, float chi, float * * mat, float * vec );

void angle_bk(float *a1, float *a2, float *b1, float *b2, float & ang);
void angle_bk(double *a1, double *a2, double *b1, double *b2, double & ang);
/// Aminoacid prototype
typedef struct {
public:

  char atom_name[5];
  float charge;  // [natoms]
  int fullatom_type;   //  [natoms]
  int nbonded_neighbors; // number of bonded atoms per atom (intra residue bonding)
  int bonded_neighbor[4];


  //float vdw; //Incluido por Nacho a partir de las tablas de Julio


  bool hastemplate; //the atom is heavy
  // template atoms used for placing atoms
    int ta [3];
  // icoor_complete - xyz coordinates of template amino acids
  float icoor[3];

  int numHydrogens_atm;
  int *hydrogens_atm;

} t_aa_atom_p;


///List of Protityped aminoacids
typedef struct {
  public:
    char aa_name3[4];
    char aa_name1;

// Booleans for detection
     bool       is_protein;
     bool       is_RNA;
     bool       is_DNA;

// Booleans basic properties
     bool       aa_is_polar;
     bool       aa_is_nonpolar;
     bool       aa_is_aromatic;
     bool       aa_is_charged;

// lista de propiedades de cada aa

    int	        natoms; // number of total atoms
    int 	nheavyatoms; // number of heavy atoms
    int 	nchi; // number of chi angles
    float       mass; // total mass
    float       nelec; // total number of electrons

   t_aa_atom_p *atom;

// lista chi parameters
   entero4 *chi_atoms;
   bool **chi_required;
   int *chi_types;

  // atom number for backbone HN
  int HNpos;
  // atom number for backbone HA
  int HApos;

  // number of polar hydrogens
  int nH_polar_complete;
  // atom numbers for polar H
  int *Hpos_polar_complete;

  // number of aromatic hydrogens
  int nH_aromatic_complete;
  // atom number for apolar hydrogens
  int *Hpos_aromatic_complete;

  // number of apolar hydrogens
  int nH_apolar_complete;
  //atom number for apolar hydrogens
  int *Hpos_apolar_complete;

  // number of acceptors
  int nacceptors;
  //acceptor information
  Tpair *aBase;
  Tpair *accpt_pos;

  //Number of Hydrogens connection
  int nH_hydrogen_connexions;
  //atoms hydrogens are connected too
  Tpair *Hbase;

  // number of atoms for EEF1 potencial
  int natoms_EEF1;

} t_aa;

///Initialization of a list with all the descriptions of the different types of aminoacids
///@param opt: Rosseta or ICM
void init_aminoacids(Convention opt=Rosseta, ConventionNames opt2=iupac);

// number to convention name
Convention int2cff (int name_convention);

// number to convention name
ConventionNames int2cName (int name_convention);




///Initialization of a list with all the descriptions of the different types of nucleotids
void nucleotids(Convention opt=Rosseta); // Mon added (1/12/2009)

#ifndef RESINI_H
#define RESINI_H

int resnum_from_resname(char *);
// Returns the residue name (in 3/4 characteres format). Watch out! You should allocate 5 chars (last one is '\0')
void resname_from_resnum(int, char *);



extern float VDW_C[N_ATOM_MAX][MAX_VDW];
extern float VDW_H[N_ATOM_MAX][MAX_VDW];
extern float ELE[MAX_ELE];

extern t_aa AA[N_AMINO];
///Global variable to know if the aminoacid table is created
extern bool aminoacids_is_init;
///Inicialization of a table with diferent Carbon-large VDW factors for all the atoms at different distances
void init_vdw_c(float probeRad=1.8,float probeEmax=1.2,float r=1.0,float f1=1.0,float f2=1.0);
//void init_vdw_c(float radius);
///Inicialization of a table with diferent Hydrogen VDW factors for all the atoms at different distances
void init_vdw_h(void);
///Inicialization of a table with diferent Electrostatic factors for all the atoms at different distances
void init_ele(void);

void aa_iupac2pdb();
void aa_iupac2pqr();
void aa_pdb2iupac();

void change_charges_withoutH();
#endif
