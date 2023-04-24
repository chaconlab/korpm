
#ifndef MOLECULE_H
#define MOLECULE_H

#include "CommonTypes.h"
#include "container.h"
#include "Chain.h"




/**
* A Molecule. Any kind of structure formed by atoms (Proteins, ADN, ARN, etc)
*/
 class Molecule: public PDB_Container
{
 public:
  ///Simple Constructor
  Molecule();
  /**
  * Simple Constructor
  * @param name: Name of the Molecule
  * @param i_nid: Molecule identifier
  */
  Molecule(Tname name,int i_nid=0);
  /**
  * Copy Constructor
  * @param old: Reference Molecule to be copied
  * @param with_elements: If true , it copies the array of elements in the buffer of the new object, otherwise the buffer is empty
  */
  Molecule(Molecule *old,bool with_elements=true);
  /// Destructor
  virtual ~Molecule();
  /// Returns the Molecule type that contains the object
  virtual TMOL getMolType()=0;
};


#endif
