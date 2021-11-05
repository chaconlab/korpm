/*
 * NAcid.h
 *
 *  Created on: Feb 4, 2010
 *      Author: nacho
 */

#ifndef NACID_H_
#define NACID_H_

#include "Molecule.h"


/// A NucleAcid: A Molecule formed by Chains of Nucleotids. It contains Chains as
/// sub-elements
/// For more information about functionability see classes Contained and Container

class NAcid:public Molecule{

 public:
    ///Simple Constructor
	 NAcid();
   /**
   * Simple Constructor
   * @param name: Name of the Protein
   * @param i_nid: Nucleacid identifier
   */
	 NAcid(Tname name,int i_nid=0);
   /**
   * Copy Constructor
   * @param old: Reference Nucleacid to be copied
   * @param with_elements: If true , it copies the array of elements in the buffer of the new object, otherwise the buffer is empty
   */
	 NAcid(NAcid *old,bool with_elements=true);

  ///Return the class of the object
  TElement getClass();
  /// Returns the Molecule type that contains the object
  TMOL getMolType();
};
#endif /* NACID_H_ */
