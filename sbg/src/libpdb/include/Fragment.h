
#ifndef FRAGMENT_H
#define FRAGMENT_H


/*#include <string.h>*/
#include <stdlib.h>
#include "Atomo.h"

typedef float( * COORDS ) [3];

/**
* A fragment. A set of atoms. Any minimum structure that contains atoms must be
* a sub-class of this (Residue,etc)
* For more information about functionability see classes Contained and Container
*/
 class Fragment : public PDB_Container {
 protected:
	int npos;
	char letter;
	char ss;

 public:

	char fragid; // Mon: Fragment (AminoAcid or Nucleotide) unique IDentifier

   ///Simple Constructor
   Fragment();
   ///Constructor
   /**
    *@param name: Name of the object
    *@param i_nid: Fragment identifier
    *@param i_npos: Fragment position
    *@param i_letter: Associated letter in PDB File
    *@param t: Molecule datatype
    */
   Fragment(Tname name,int i_nid=0,int i_npos=0,char i_letter=' ');
   ///Copy Constructor
   ///@param old: Reference Fragment. The new object will be a exact copy of this object
   ///@param with_elements: If true , it copies the array of elements in the buffer of the new object, otherwise the buffer is empty
   Fragment(Fragment *old,bool with_elements=true);
   /// Destructor
   virtual ~Fragment();
	///Returns the order position of the Fragment
	int get_pos();
	///Returns the associated letter of the fragment in the PDB file
	char get_letter();
	///Returns secondary structure character
	char get_ss();
	///Sets secondary structure character
	void set_ss(char ch);
	/// Returns the Molecule type that contains the object
	virtual TMOL getMolType()=0;


};
#endif

