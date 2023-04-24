
#ifndef ATOMO_H
#define ATOMO_H


#include <stdio.h>
#include <math.h>
#include <string.h>


#include "CommonTypes.h"
#include "container.h"
///         atomo types for energy calculation....

#define NUM_ELEMENTS 29 // Added by Mon (28/12/2012)

//car hybridization
int const H_HYBRID = { 0 };
int const SP2_HYBRID = { 1 };
int const SP3_HYBRID = { 2 };
int const RING_HYBRID = { 3 };

// polar...
int const APOLAR = { 0 };
int const POLAR = { 1 };




// atom types
/**
* *@author Jose Ignacio Garzon
*
* Store the information of a atom that could be found in an aminoacid
*/
typedef  struct  {

  /// type of atom
  int at;
  /// Nombre
  char name[5];
  /// Hybrid
  int hyb;
  /// Polar
  int polar;
  /// Charge
  float chrge;
  /// vdw full
  float vdw;
  /// soft
  float vdws;
  ///deep
  float vdwd;
  // acceptor
  bool acp;
  // donor
  bool dor;

	//NEW TYPES FOR SOLVATATION
	float Avol;
	float Asolpar;

}  Atom_type;

///Pointer to the selected list of different atom types in the aminoacids
extern Atom_type *atom_types;
///List of different atom types in Rosseta convention
extern Atom_type atom_types_Rosseta[54];
///List of different atom types in ICM convention
extern Atom_type atom_types_ICM[54];
///List of different atom types in EEF1 convention
extern Atom_type atom_types_EEF1[31];
///List of different atom types in Sybil (Rosseta Variation) convention
extern Atom_type atom_types_Sybil[45];
/// Number of atom types
extern int num_atom_type;

/**
* Representation of a basic element(Oxygen, Carbon, etc.)
*/
 class Element {
 public:
  /*Element attributes*/
  /// Name
  Tname name;
  /// Symbol
  Symbol symbol;
  char sym[2];
  /// Group
  int group;
  /// Period
  int period;
  /// Number
  float number;
  /// weight
  float weight;
  /// Numero of atoms
  float atomic;
  /// Covalence
  float cov;
  /// Van Der Walls
  float vdw;
  float en;


  ///Basic constructor. Do nothing
  Element()
  {
  };

/**Constructor
* @param n: name of the Element
* @param s: symbol
* @param s2: full symbol
* @param g: group
* @param w: weight
* @param p: period
* @param num: Number of atoms
* @param a: atomic weight???
* @param c: Covalence
* @param e: Energy
*/
 Element(Tname n,Symbol s,char s2[2], int g, int p, float num, float w, float a,float c, float v,float e)
  {
    strcpy(name,n);
    symbol=s;
    strcpy(sym,s2);
    group=g;
    period=p;
    number=num;
    weight=w;
    atomic=a;
    vdw=v;
    en=e;
    cov=c;

  };

 ///Constructor. Copy of another element
 ///@param ele: Reference object to be copied
  Element(Element *ele)
  {
    strcpy(name,ele->name);
    symbol=ele->symbol;
    strcpy(sym,ele->sym);
    group=ele->group;
    period=ele->period;
    number=ele->number;
    weight=ele->weight;
    atomic=ele->atomic;
    vdw=ele->vdw;
    en=ele->en;
  };
  ///Destructor
  ~Element()
  {
  };
};



/**
* Static class with the collection of different basic elements
*/
class Table_Elements{


 public:
   ///Access to an element of the table by position
   static Element* getElement (int index)
   {
     return &table_Elements[index];
   };

  ///Access to an element of the table by symbol
   static Element* getElement (char sym[3] )
   {
     int i;
     Element *e;
     for(i=0;i<NUM_ELEMENTS;i++)
     {
      e=getElement(i);
//      fprintf(stderr,"getElement> e->sym= %s   sym= %s\n",e->sym, sym);
	  if(strcmp(e->sym,sym)==0)
        return &table_Elements[i];
     }
     //fprintf(stderr,"getElement> ERROR:Element does not exist: %s \n",sym);
     return NULL;
	 /*if(strcmp(sym,"C ")==0)
       return table_Elements[0];
     if(strcmp(sym,"H ")==0)
       return table_Elements[1];
     if(strcmp(sym,"N ")==0)
       return table_Elements[2];
     if(strcmp(sym,"O ")==0)
       return table_Elements[3];
     */
   };

 private:
   /// Static table of elements
   static Element table_Elements[];

};

class Bond;

/**
* An Atom. It stores information about an Atom: Element that compounds it, spatial position, etc.
* The class also allows to modify this information
* For more information about functionability see classes Contained
*/
class Atom : public PDB_Contained{
public:

  ///Simple Constructor
  Atom();
  /// Constructor that sets all the information about the atom
  Atom(Element *e, Tcoor p,float ch,char name[5],int serial, float occ, float fact);
  /// Constructor by template
  Atom(Atom_type atom_type, char name[5], float *p, int serial);
  /// Copy Constructor. Makes a exact copy of another atom
  Atom(Atom *a);

  ///Destructor
  virtual ~Atom();

  ///Returns spatial position of the atom
  void getPosition(Tcoor coor);

  ///Returns charge of the atom
  float getCharge();
  ///Returns a pointer to the Element that the atom is made of
  Element* getElement();
  ///Returns the serial number of the atom in the system
  int getPdbSerial();
  ///returns name of the atom in the system
  char* getPdbName();
  float getPdbocc();
  float getPdbfact();


  /*Modify information*/
  ///Sets the position of the atom
  void setPosition(Tcoor pos);
  ///Sets the charge of the atom
  void setCharge(float ch);
  ///Sets the serial number of the atom in the system
  void setPdbSerial(int serial);
  ///Sets the serial number of the atom in the system
  void setPdbName(char n[5]);
  ///Sets value in occupancy column
  void setPdbocc(float occ);
  ///Sets value in bfact column
  void setPdbfact(float fact);
  ///Sets the Element of the atom
  void setEtype(int etype);
  ///Translates the atom an offset in the Space
  bool move(Tcoor offset);
  ///Negative Translation of the atom an offset in the Space
  bool moven(Tcoor offset);

  ///Returns the Name of the atom
  char *getName();
  /// Initialization of list of Elements. Atom has no elements, so it
  /// only returns  true
  bool initAll();
  /// Move all elements. Calls to funtion move
  bool moveAll(Tcoor offset);
  /// Move all elements. Calls to funtion moven
  bool moveAlln(Tcoor offset);
  ///Returns the identifier of the class. In this case: Atom
  TElement getClass();
  ///Returns the distance from the current atom to the referenced by argument
  inline float dist(Atom *a2)
  {
    Tcoor pos2;

    a2->getPosition(pos2);
    return sqrt( pow( position[0] - pos2[0], 2 ) +
                 pow( position[1] - pos2[1], 2 ) +
                 pow( position[2] - pos2[2], 2 ) );

  };

  // PABLO, ESTO NO FUNCIONA...
  ///Modify the position of the atom by a rotation and a translation
//  inline void rotate_trans(double **R, Tcoor t)
//  {
//  	position[0] = R[0][0]*position[0] + R[0][1]*position[1] + R[0][2]*position[2] + t[0];
//      position[1] = R[1][0]*position[0] + R[1][1]*position[1] + R[1][2]*position[2] + t[1];
//  	position[2] = R[2][0]*position[0] + R[2][1]*position[1] + R[2][2]*position[2] + t[2];
//  };

  /// Modify the position of the atom by a rotation and a translation
  inline void rotate_trans(double **R, Tcoor t, Tcoor p)
  {
	  p[0] = R[0][0]*position[0] + R[0][1]*position[1] + R[0][2]*position[2] + t[0];
	  p[1] = R[1][0]*position[0] + R[1][1]*position[1] + R[1][2]*position[2] + t[1];
	  p[2] = R[2][0]*position[0] + R[2][1]*position[1] + R[2][2]*position[2] + t[2];
	  position[0] = p[0];
	  position[1] = p[1];
	  position[2] = p[2];
  };


  /// Returns the Molecule type that contains the object
  TMOL getMolType();

  int get_numBonds();
  Bond* getBond(int i);
  void insertBond(Bond *b);
  bool removeBond(Bond *b);


private:

  /*Attributes of the Atom*/
  /// Element
  Bond **bonds;
  int num_bonds;
  Element *element;

  ///Charge
  float charge;

  /// Space position
  Tcoor position;

  /// Serial number in the system
  int pdbSerialNumber;
  /// Name
  char pdbName[5];
  float pdbocc;
  float pdbfact;

};


class Bond{
private:
	Atom *init;
	Atom *final;
	int link;
public:
	Bond(Atom* i, Atom *f, int l);
	~Bond();
	int getLink();
	Atom *getInit();
	Atom *getFinal();
};

#endif
