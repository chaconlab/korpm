
#ifndef PROTEIN_H
#define PROTEIN_H


#include "Molecule.h"


/// A Protein: A Molecule formed by Chains of Aminoacids. It contains Chains as
/// sub-elements
/// For more information about functionability see classes Contained and Container

 class Protein:public Molecule{

 public:
    ///Simple Constructor
   Protein();
   /**
   * Simple Constructor
   * @param name: Name of the Protein
   * @param i_nid: Protein identifier
   */
   Protein(Tname name,int i_nid=0);
   /**
      *  Copy Constructor
      *  @param old: Reference Protein to be copied
      *  @param with_elements: If true , it copies the array of elements in the buffer of the new object, otherwise the buffer is empty
      */
   Protein(Protein *old, bool with_elements=true);



  ///Return the class of the object
  TElement getClass();
  /// Returns the Molecule type that contains the object
  TMOL getMolType();
};

#endif
