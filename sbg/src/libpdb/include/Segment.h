//#ifndef SEGMENT_H
//#define SEGMENT_H

#include "Residue.h"
#include "Nucleotide.h"


 /**
 * A segment: It represents a a full sub-Chain (not broken) of a chain of
 * aminoiacids. It contains Residues (Aminoacids) as sub-elements
 * To see more information about functionability. See classes Contained and PDB_Container
 */
 class Segment: public PDB_Container {

 public:
   /// Simple Constructor
   Segment();
   /// Simple Constructor
   /// @param name: Name of the Segment
   /// @param in_nid: Segment identifier
   /// @param t: Molecule datatype
   Segment(Tname name, int in_nid=0);
   /// Copy Constructor
   ///@param old: reference fragment to be copied
   ///@param with_elements: If true , it copies the array of elements in the buffer of the new object, otherwise the buffer is empty
   Segment(Segment *old,bool with_elements=true);
  ///Destructor
  virtual ~Segment();
  ///Return the class of the object
  TElement getClass();
  /// Returns the Molecule type that contains the object
   TMOL getMolType();
};
//#endif
