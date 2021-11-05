/***************************************************************************
                          Condition.h  -  description
                             -------------------
    begin                : Fri Aug 6 2004
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
#ifndef CONDITION_H
#define CONDITION_H

/**
* Representation of a condition to select Elements in a system. Selects atoms
* in function of the name and the position in the lists of Elements
*/
class Condition{
 private:
  ///List of valid names for the selected atoms
  char **atoms;
  ///Number of names in the list "atoms"
  int limit;
  ///Searches limits in the collections of Proteins, Chains, Residues, and Atoms
  int lowP,highP,lowS,highS,lowCh,highCh,lowR,highR,lowA,highA;
  /**
  *Checks if a position is valid in the list of Proteins
  *@param pos: Position to check
  */
  bool isValidP(int pos);
  /**
  *Checks if a position is valid in the list of Segments
  *@param pos: Position to check
  */
  bool isValidS(int pos);
  /**
  *Checks if a position is valid in the list of Chains
  *@param pos: Position to check
  */
  bool isValidCh(int pos);
  /**
  *Checks if a position is valid in the list of Residues
  *@param pos: Position to check
  */
  bool isValidR(int pos);
  /**
  *Checks if a position is valid in the list of Atoms
  *@param pos: Position to check
  */
  bool isValidA(int pos);

 public:
 /**
 * Constructor
 * @param lP,hP,lS,hS,lCh,hCh,lR,hR,lA,hA: Search limits
 * WARNING: Segment (lS,hS) goes before Chain (lCh,hCh), and this should be reversed!!!
 */

  Condition(int lP,int hP,int lS,int hS,int lCh,int hCh,int lR,int hR,int lA,int hA);
  ~Condition();
  /**
  * Introduces names of atoms to search
  * @param cad: list of atoms
  */
  void add(char *cad);

  /**
  * Checks if a atom name and a position are valid for the selection
  * @param posP,posCh,posR,posA: Positions in the list of elements
  * @param cad: Name of the atom
  */
  bool fitCondition(int posP,int posS,int posCh, int posR, int posA, char *cad);

 } ;

/**
* Collection of elements of class Condition
*/
 class Conditions{
  private:
    ///List of elements of class Condition
    Condition **list;
  public:
    ///Number of conditions stored
    int limit;

    ///Constructor
    Conditions();
    ///Destructor
    ~Conditions();
    ///Adds a Condition to the collection
    ///@param c: Pointer to Condition
    void add(Condition *c);
    ///Returns pointer to the ith Condition in the collection
    ///@param i: Position of the returned Condition
    Condition *getCondition(int i);
 };


#endif
