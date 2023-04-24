/***************************************************************************
                          container.h  -  description
                             -------------------
    begin                : Wed Jul 21 2004
    copyright            : (C) 2004 by Jose Ignacio Garzon
    email                :
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef CONTAINER_H
#define CONTAINER_H



#include <string.h>
#include <stdlib.h>
#include "CommonTypes.h"

enum TElement {pdb_null,pdb_world,pdb_macro,pdb_protein,pdb_nacid,pdb_smol, pdb_chain,pdb_segment,pdb_residue,pdb_nucleotide,pdb_atom};
enum TMOL {tmol_null,tmol_protein,tmol_nacid,tmol_smol};
/**
  *@author Jose Ignacio Garzon
  *
  *Virtual Class. All the Class to store information about the different parts
  * of a macromsystem (Protein,Chain, Residue, Atom, etc) must be subClass of this.
  * This Class declares all the functions neccesary by a Contained Element to be accesible
  * For more information about functionability see classes Contained and Container
 */
class PDB_Contained {
protected:
  ///Pointer to Father
  void *father;
public:

  virtual ~PDB_Contained();

  virtual char *getName()=0;
  ///Virtual Function. Must move the element (and the sub-elements)
  ///@param offset: Movement to be applied
  virtual bool moveAll(Tcoor offset)=0;
  ///Virtual Function. Must move the element (and the sub-elements)
  ///@param offset: Inverse Movement to be applied
  virtual bool moveAlln(Tcoor offset)=0;
  ///Virtual Function. Must return the name of the class
  virtual TElement getClass()=0;
  /// Virtual Function. Sets the pointer in the first sub-element of the list. The sub-elements
  /// , if they have sub-sub-elements, points to the first sub-sub-element also
  virtual bool initAll()=0;
  /// Returns the Molecule type that contains the object
  virtual TMOL getMolType()=0;
  /// Returns the pointer to the Container object that contains this object
  void *getFather();
  /// Sets the pointer to the Container object that contains this object
  bool setFather(void *f);
  };

#include "Atomo.h"


/**
*@author Jose Ignacio Garzon
*
* Implements all the funtionability of the classes that store sub-elements
* (Proteins contains Chains, Chains contains Segments, etc.)
*/

class PDB_Container: public PDB_Contained {
protected:
  ///Identifier of the obhject
  Tname id;
  ///Identifier number
  int nid;
  ///List of sub-lements
  PDB_Contained **elements;

  int limit;
  ///Current Element that is pointed
  int currentE;
  /// True if the sub-elments are not a sub-Class of Container. False in another case
  bool leafs;
  /// Returns the position of the first element with the name indicated
  ///param name: name to be search
  int getPosition(Tname name);

public:

  //~PDB_Container();
  ///Introduces a new sub-element at the end of the list
  ///@param newE: Pointer to the new sub-element
  bool add(PDB_Contained *newE);
  ///Returns the sub-element currently pointed
  void *getCurrent();
  ///Points to the last sub-element
  int end();
  /// Erases the sub-element in the position. The sub-element is completly lost
  /// @param pos: position of the sub-element to be deleted
  bool erase(int pos);
  /// Erases the sub-element with the name indicated. The sub-element is completly lost
  /// @param pos: name of the sub-element to be deleted
  bool erase(Tname name);
  /// Erases all the sub-elements of the object
  bool eraseAll();
  /// Returns the number of sub-elements
  int getLimit();
  /// Returns a pointer to the sub-element in the position
  ///@param pos: Position of the sub-element returned
  PDB_Contained *getE(int pos);
  /// Sets the pointer in the first sub-element of the list
  bool init();
  /// Sets the pointer in the first sub-element of the list. The sub-elements
  /// , if they have sub-sub-elements, points to the first sub-sub-element also
  bool initAll();
  /// applies a movement to all the sub-elements
  /// @param Movement to be applied
  bool moveAll(Tcoor offset);
  /// applies a movement to all the sub-elements
  /// @param Inverse Movement to be applied
  bool moveAlln(Tcoor offset);
  /// Points to the next sub-element
  bool next();
  /// Points to the previous sub-element
  bool previous();
  /// removes from the list a sub-element. The sub-element is not erased
  /// @param pos: position of the sub_element removed
  bool remove(int pos);
  /// removes from the list a sub-element. The sub-element is not erased
  /// @param pos: name of the sub_element removed
  bool remove(Tname name);
  /// Removes the total list of sub-elements. The sub-elements are not deleted
  bool removeAll();
  /// Erases all the structure over a level. Type level cannot be the same as the level of the object
  ///@param level: Class type that must not be deleted (subelements will be kept also)
  bool erase_level(TElement level);
  /// Returns the name of the object
  char *getName();
  /// Changes the name of the object
  /// @param name: new name
  void modifyName(char *name);
  ///Changes positions of two sub-elemens
  ///@param: i,j: positions of the subelements to be changed
  void exchange(int i, int j);
  /// Returns the id number
  int getIdNumber();
  /// Checks whether the types of the elements fits with the Molecules where they are included
  bool check_consistency();
  ///Returns the number of atoms included in the container
    int num_atoms();
  ///Erases all the hydrogen included in the container
  void delete_hydrogens();
  ///Erases all hetero
  void delete_hetero();
  ///Erases all waters
  void delete_water();
  // Delete contained duplicate atoms preserving the first occurrence (Mon 21/02/2018)
  void delete_duplicates();

  /// Returns the Molecule type that contains the object
  virtual TMOL getMolType()=0;

};

#endif
