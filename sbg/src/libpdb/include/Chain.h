
#ifndef CHAIN_H
#define CHAIN_H

#include "Segment.h"


 /**
 * A Chain: It represents a chain of aminoacids. It contains Segments as
 * sub-elements because it is possible that a chain is "broken" (a aminoacid is lost)
 * For more information about functionability see classes Contained and Container
 */
 class Chain: public PDB_Container {

 public:
   /// Simple Constructor
   Chain();
   /// Simple Constructor
   ///@param name: Name of the chain
   ///@param i_nid: Chain identifier
   ///@param t: Molecule datatype
   Chain(Tname name,int i_nid=0 );
   /**
    *  Copy Constructor
    *  @param old: Reference Chain to be copied
    *  @param with_elements: If true , it copies the array of elements in the buffer of the new object, otherwise the buffer is empty
    */
   Chain(Chain *old,bool with_elements=true);
  ///Destructor
  virtual ~Chain();

  ///Return the class of the object
  TElement getClass();
  /// Returns the Molecule type that contains the object
  TMOL getMolType();

};
#endif
