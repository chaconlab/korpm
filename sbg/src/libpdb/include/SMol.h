/*
 * SMol.h
 *
 *  Created on: Feb 4, 2010
 *      Author: nacho
 */

#ifndef SMOL_H_
#define SMOL_H_

#include "Molecule.h"
/// A Small Molecule: A Molecule formed by acollection of heteroatoms
/// For more information about functionability see classes Contained and Container

class SMol:public Molecule{
 protected:
	int npos;
	Bond **bonds;
	int num_bonds;
 public:
    ///Simple Constructor
	SMol();
   /**
   * Simple Constructor
   * @param name: Name of the Protein
   * i_nid: Small Molecule identifier
   * i_pos: Order position
   */
	SMol(Tname name,int i_nid=0,int i_pos=0);
   /**
   * Copy Constructor
   * @param old: Reference small molecule to be copied
   * @param with_elements: If true , it copies the array of elements in the buffer of the new object, otherwise the buffer is empty
   */
	SMol(SMol *old,bool with_elements=true);

	~SMol();

	///Returns the order position of the Molecule
	int get_pos();
	///Return the class of the object
	TElement getClass();
	/// Returns the Molecule type that contains the object
	TMOL getMolType();

	//float maxLength();

  void insertBond(Bond *b);
  Bond *getBond(int i);
  bool deleteBond(Bond *b);
  int get_numBonds();
  void deleteHydrogenBonds();


    ///Maximun length between atoms of the system to the Origin of coordinates.
};

#endif /* SMOL_H_ */
