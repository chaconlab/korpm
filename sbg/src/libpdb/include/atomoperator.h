/***************************************************************************
                          atomoperator.h  -  description
                             -------------------
    begin                : Thu May 6 2004
    copyright            : (C) 2004 by
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

#ifndef ATOMOPERATOR_H
#define ATOMOPERATOR_H

#include "Atomo.h"


/**
* Virtual class. All the operations over Atom must inherit from this one
*/
class AtomOperator {
public:
        ///Constructor
	AtomOperator();
        ///Destructor
virtual ~AtomOperator()=0;

///Application of the operator over the Atom passed by argument
virtual bool apply(Atom *a)=0;

/**
* Application of the operator on a referenced atom b. The result of the operation
* is applied over the atom a
*
* @param a: Resulted atom
* @param b: Initial atom
*/
virtual bool apply(Atom *a,Atom *b)=0;
};


/**
* Rotation Operator. The atom is rotated
*/
class EulerRot:public AtomOperator{
public:
  /**
  *Constructor using Euler angles
  *
  *@param psi,theta,phi: Euler angles
  *@param opcion: Convention of Euler angles. 0 (ZXZ) 1 (XYZ) 2 (ZYZ).
  */
  EulerRot(float psi, float theta, float phi,int opcion=0);
  ///Destructor(ZXZ)
  ~EulerRot();
  /**
   * Application of the operator
   *@param a: atom to be rotated
   */
  bool apply(Atom *a);
  ///Application of the operator from a reference atom to another atom
  ///@param a: atom to be rotated
  ///@param b: atom that gives the initial position
  bool apply(Atom *a,Atom *b);

//private:
  /// Rotation matrix (Mon: it has to be public to easily output rotation matrix if needed)
  float matrix[3][3];
};

/**
* Rotation Operator through one axis
*/
class AxisRot:public AtomOperator{
public:
  /**
  * Constructor
  *
  * @param _axis: axe through the atom is moved
  * @param _angle: Grades of the rotation
  */
  AxisRot(char _axis, float _angle);
  ///Destructor
  ~AxisRot();
  /**
   * Application of the operator
   *@param a: atom to be rotated
   */
  bool apply(Atom*a);
  /**
   * Application of the operator from a reference atom to another atom
   *@param a: atom to be rotated
   *@param b: atom that gives the initial position
   */
  bool apply(Atom *a,Atom *b);
private:
  /// axe through the atom is moved
  char axis;
  /// Grades of the rotation
  float angle;
};



/**
* Rotation and translation operator using a Rotation Matrix
*/
class M4Rot: public AtomOperator{
public:
   ///Constructor
   ///@param in: Rotation/translation Matrix to apply
   M4Rot(float in[4][4]);
  /**
   * Constructor
   *@param rot: Euler angles that define the rotation
   *@param option: Euler convention.0 (ZXZ) 1 (XYZ) 2 (ZYZ).
   */
  M4Rot(float *rot, int option,float *offset);
  ~M4Rot();
  ///Application of the operator
  ///@param a: atom to be rotated
  bool apply(Atom *a);
  ///Application of the operator from a reference atom to another atom
  ///@param a: atom to be rotated
  ///@param b: atom that gives the initial position
  bool apply(Atom *a,Atom *b);

private:
  ///Rotation Matrix
  float matrix[4][4];

};





#endif
