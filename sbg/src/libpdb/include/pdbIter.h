#ifndef PDBITER_H
#define PDBITER_H

#include "Macromolecule.h"

/**
*Iterator to fast access to the different elements (Atoms, Residues, chains, etc)
* of a system of Macromolecules compounded by proteins
*/
class pdbIter{
public:
  ///Simple constructor
  pdbIter();
  /**
   * Constructor over a Container Object
   *@param: c: Object over the Iterator is created. The iterator allows access to all elements and sub-elements of this object
   *@param: include_smol: Determines if small molecules (hereoatoms) must be included in the iterator
   *@param: include_nacid: Determines if nucleacids must be included in the iterator
   *@param: include_protein: Determines if proteins must be included in the iterator
   */
  pdbIter(PDB_Container *c, bool include_smol=true,bool include_nacid=true, bool include_protein=true, bool in_special_iter=false);
  ///Copy contructor
  ///@param i2: Reference iterator to be copied
  pdbIter(pdbIter *i2);
  ///Destructor
  ~pdbIter();
  //Delete virtual segments if created
  void clean_virtual();
  ///Increases the accessible elements width the elements of other iterator
  ///@param other: Iterator which the the elements are taken from
  void add(pdbIter *other);
  ///Number of accessible macromoelecules
  int num_macro();
  ///Number of accessible chains
  int num_chain();
  ///Number of accessible molecules
  int num_molecule();
  ///Number of accessible segments
  int num_segment();
  ///Number of accessible fragments
  int num_fragment();
  ///Number of accessible atoms
  int num_atom();

   ///returns current accessible Macromolecule
   Macromolecule* get_macro();
   ///returns current accessible chain
   Chain* get_chain();
   ///returns current accessible Molecule
   Molecule* get_molecule();
   ///returns current accessible Segment
   Segment* get_segment();
   ///returns current accessible Fragment
   Fragment* get_fragment();
   ///returns current accessible Atom
   Atom* get_atom();

  ///Current accessible Macromolecule
  int pos_macro;
 ///Current accessible Chain
  int  pos_chain;
  ///Current accessible Molecule
  int  pos_molecule;
  ///Current accessible Segment
  int  pos_segment;
  ///Current accessible Fragment
  int  pos_fragment;
  ///Current accessible Atom
  int  pos_atom;

   ///Points to the next accessible Macromolecule
   bool next_macro();
   ///Points to the next accessible Chain
   bool next_chain();
   ///Points to the next accessible Molecule
   bool next_molecule();
   ///Points to the next accessible Segment
   bool next_segment();
   ///Points to the next accessible Fragment
   bool next_fragment();
   ///Points to the next accessible Atom
   bool next_atom();

   /// Returns true if the pointer of the list of accessible Macromolecules is at the end
   bool gend_macro();
   /// Returns true if the pointer of the list of accessible Chains is at the end
   bool gend_chain();
      /// Returns true if the pointer of the list of accessible Molecules is at the end
   bool gend_molecule();
      /// Returns true if the pointer of the list of accessible Segments is at the end
   bool gend_segment();
      /// Returns true if the pointer of the list of accessible Fragments is at the end
   bool gend_fragment();
      /// Returns true if the pointer of the list of accessible Atoms is at the end
   bool gend_atom();

  ///Forward movement of the pointer of accessible atoms
  bool operator ++ (int);
  ///Backward movement of the pointer of accessible atoms
  bool operator -- (int);

  /// Returns true if the pointer of the list of accessible Macromolecules is at the beginning
   bool ginit_macro();
   /// Returns true if the pointer of the list of accessible Chains is at the beginning
   bool ginit_chain();
   /// Returns true if the pointer of the list of accessible Molecules is at the beginning
   bool ginit_molecule();
   /// Returns true if the pointer of the list of accessible Segments is at the beginning
   bool ginit_segment();
   /// Returns true if the pointer of the list of accessible Fragments is at the beginning
   bool ginit_fragment();
   /// Returns true if the pointer of the list of accessible Atoms is at the beginning
   bool ginit_atom();

   ///Points to the previous accessible Macromolecule
   bool back_macro();
   ///Points to the previous accessible Chain
   bool back_chain();
   ///Points to the previous accessible Molecule
   bool back_molecule();
   ///Points to the previous accessible Segment
   bool back_segment();
   ///Points to the previous accessible Fragment
   bool back_fragment();
   ///Points to the previous accessible Atom
   bool back_atom();


  private:
  /// List of Accessible Macromolecules
  Macromolecule **L_macro;
  /// List of Accessible Chains
  Chain **L_chain;
  /// List of Accessible Molecules
  Molecule **L_molecule;
  /// List of Accessible Segments
  Segment **L_segment;
  /// List of Accessible Fragments
  Fragment **L_fragment;
  /// List of Accessible Atoms
  Atom **L_atom;

  /// Number of accessible Macromolecules
  int end_macro;
  /// Number of accessible Chains
  int  end_chain;
  /// Number of accessible Molecules
  int  end_molecule;
  /// Number of accessible Segments
  int  end_segment;
  /// Number of accessible Fragments
  int  end_fragment;
  /// Number of accessible Atoms
   int end_atom;

   ///Introduce SMOL as Fragments inside virtual segments.
   bool special_iter;

};




#endif
