/***************************************************************************
                          atomoperator.cpp  -  description
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

#include "atomoperator.h"
#include <math.h>
#include <stdio.h>


using namespace std;

AtomOperator::AtomOperator(){
}
AtomOperator::~AtomOperator(){
}

EulerRot::EulerRot(float psi, float theta, float phi,int opcion){
  float sin_psi = sin(psi);
  float cos_psi = cos(psi);
  float sin_theta= sin(theta);
  float cos_theta= cos(theta);
  float sin_phi=sin(phi);
  float cos_phi=cos(phi);

  if(opcion==0)
  {
    matrix[0][0]= cos_psi*cos_phi - cos_theta*sin_phi*sin_psi;
    matrix[0][1]= cos_psi*sin_phi + cos_theta*cos_phi*sin_psi;
    matrix[0][2]= sin_psi*sin_theta;
    matrix[1][0]= -sin_psi*cos_phi - cos_theta*sin_phi*cos_psi;
    matrix[1][1]= -sin_psi*sin_phi + cos_theta*cos_phi*cos_psi;
    matrix[1][2]= cos_psi*sin_theta;
    matrix[2][0]= sin_theta*sin_phi;
    matrix[2][1]= -sin_theta*cos_phi;
    matrix[2][2]= cos_theta;
    return;
  }
  if(opcion==1)
  {
    matrix[0][0]= cos_theta*cos_phi;
    matrix[0][1]= cos_theta*sin_phi;
    matrix[0][2]= -sin_theta;
    matrix[1][0]= sin_psi*sin_theta*cos_phi - cos_psi*sin_phi;
    matrix[1][1]= sin_psi*sin_theta*sin_phi + cos_phi*cos_psi;
    matrix[1][2]= cos_theta*sin_psi;
    matrix[2][0]= cos_psi*sin_theta*cos_phi + sin_psi*sin_phi;
    matrix[2][1]= cos_psi*sin_theta*sin_phi - sin_psi*cos_phi;
    matrix[2][2]= cos_theta*cos_psi;
    return;
  }
  if(opcion==2)
  {
    matrix[0][0]= -sin_psi*sin_phi+cos_theta*cos_phi*cos_psi;
    matrix[0][1]= sin_psi*cos_phi+cos_theta*sin_phi*cos_psi;
    matrix[0][2]= -cos_psi*sin_theta;
    matrix[1][0]= -cos_psi*sin_phi-cos_theta*cos_phi*sin_psi;
    matrix[1][1]= cos_psi*cos_phi-cos_theta*sin_phi*sin_psi;
    matrix[1][2]= sin_psi*sin_theta;
    matrix[2][0]= sin_theta*cos_phi;
    matrix[2][1]= sin_theta*sin_phi;
    matrix[2][2]= cos_theta;
    return;
  }
  printf("mrot> Error: Bad option creating matrix\n");


}

EulerRot::~EulerRot(){
}

bool EulerRot::apply(Atom *a){

  Tcoor newPos;
  Tcoor pos;
  a->getPosition(pos);
  float currx=pos[0];
  float curry=pos[1];
  float currz=pos[2];

  newPos[0]= currx *matrix[0][0] + curry*matrix[0][1] + currz*matrix[0][2];
  newPos[1]= currx *matrix[1][0] + curry*matrix[1][1] + currz*matrix[1][2];
  newPos[2]= currx *matrix[2][0] + curry*matrix[2][1] + currz*matrix[2][2];

  a->setPosition(newPos);
  return true;
}

bool EulerRot::apply(Atom *a,Atom *b){

  Tcoor newPos;
  Tcoor pos;
  a->getPosition(pos);
  float currx=pos[0];
  float curry=pos[1];
  float currz=pos[2];

  newPos[0]= currx *matrix[0][0] + curry*matrix[0][1] + currz*matrix[0][2];
  newPos[1]= currx *matrix[1][0] + curry*matrix[1][1] + currz*matrix[1][2];
  newPos[2]= currx *matrix[2][0] + curry*matrix[2][1] + currz*matrix[2][2];

  b->setPosition(newPos);
  return true;
}


AxisRot::AxisRot(char _axis, float _angle)
{
  axis=_axis;
  angle=_angle;

  if(axis!='X' && axis!='Y' && axis!='Z')
      fprintf(stdout,"Error_AxisRot: Incorrect coordenate identifier (X,Y,Z)\n");

}
AxisRot::~AxisRot()
{
}

bool AxisRot::apply(Atom *a)
{
  float sint= sin(angle);
  float cost= cos(angle);
  Tcoor newPos;
  Tcoor pos;
  a->getPosition(pos);
  float currx=pos[0];
  float curry=pos[1];
  float currz=pos[2];

  switch(axis)
  {
    case('X'):
        newPos[0]=currx;
        newPos[1]=cost*curry + sint*currz;
        newPos[2]=cost*currz - sint*curry;
      break;
    case('Y'):
        newPos[0]=cost*currx + sint*currz;
        newPos[1]=curry;
        newPos[2]=cost*currz - sint*currx;
      break;
    case('Z'):
        newPos[0]=cost*currx - sint*curry;
        newPos[1]=cost*curry + sint*currx;
        newPos[2]=currz;
      break;
     default:
       return false;
       break;
  }
  a->setPosition(newPos);
  return(true);


}

bool AxisRot::apply(Atom *a,Atom *b)
{
  float sint= sin(angle);
  float cost= cos(angle);
  Tcoor newPos;
  Tcoor pos;
  a->getPosition(pos);
  float currx=pos[0];
  float curry=pos[1];
  float currz=pos[2];

  switch(axis)
  {
    case('X'):
        newPos[0]=currx;
        newPos[1]=cost*curry + sint*currz;
        newPos[2]=cost*currz - sint*curry;
      break;
    case('Y'):
        newPos[0]=cost*currx + sint*currz;
        newPos[1]=curry;
        newPos[2]=cost*currz - sint*currx;
      break;
    case('Z'):
        newPos[0]=cost*currx - sint*curry;
        newPos[1]=cost*curry + sint*currx;
        newPos[2]=currz;
      break;
     default:
       return false;
       break;
  }
  b->setPosition(newPos);
  return(true);

}


M4Rot::M4Rot(float in[4][4]){

  for(int a=0;a<4;a++)
      for(int b=0;b<4;b++)
        matrix[a][b]=in[a][b];
}

M4Rot::M4Rot(float *rot, int option,float *offset){

	int i;
	float sin_psi = sin(rot[0]);
  float cos_psi = cos(rot[0]);
  float sin_theta= sin(rot[1]);
  float cos_theta= cos(rot[1]);
  float sin_phi=sin(rot[2]);
  float cos_phi=cos(rot[2]);

  if(option==0)
  {
    matrix[0][0]= cos_psi*cos_phi - cos_theta*sin_phi*sin_psi;
    matrix[0][1]= cos_psi*sin_phi + cos_theta*cos_phi*sin_psi;
    matrix[0][2]= sin_psi*sin_theta;
    matrix[1][0]= -sin_psi*cos_phi - cos_theta*sin_phi*cos_psi;
    matrix[1][1]= -sin_psi*sin_phi + cos_theta*cos_phi*cos_psi;
    matrix[1][2]= cos_psi*sin_theta;
    matrix[2][0]= sin_theta*sin_phi;
    matrix[2][1]= -sin_theta*cos_phi;
    matrix[2][2]= cos_theta;

  }
  if(option==1)
  {
    matrix[0][0]= cos_theta*cos_phi;
    matrix[0][1]= cos_theta*sin_phi;
    matrix[0][2]= -sin_theta;
    matrix[1][0]= sin_psi*sin_theta*cos_phi - cos_psi*sin_phi;
    matrix[1][1]= sin_psi*sin_theta*sin_phi + cos_phi*cos_psi;
    matrix[1][2]= cos_theta*sin_psi;
    matrix[2][0]= cos_psi*sin_theta*cos_phi + sin_psi*sin_phi;
    matrix[2][1]= cos_psi*sin_theta*sin_phi - sin_psi*cos_phi;
    matrix[2][2]= cos_theta*cos_psi;

  }
  if(option==2)
  {
    matrix[0][0]= -sin_psi*sin_phi+cos_theta*cos_phi*cos_psi;
    matrix[0][1]= sin_psi*cos_phi+cos_theta*sin_phi*cos_psi;
    matrix[0][2]= -cos_psi*sin_theta;
    matrix[1][0]= -cos_psi*sin_phi-cos_theta*cos_phi*sin_psi;
    matrix[1][1]= cos_psi*cos_phi-cos_theta*sin_phi*sin_psi;
    matrix[1][2]= sin_psi*sin_theta;
    matrix[2][0]= sin_theta*cos_phi;
    matrix[2][1]= sin_theta*sin_phi;
    matrix[2][2]= cos_theta;

  }

	for(i=0;i<3;i++)
	{
		matrix[3][i]=0;
		matrix[i][3]=offset[i];
	}

	matrix[3][3]=1.0;

}

M4Rot::~M4Rot(){
}

bool M4Rot::apply(Atom *a){

  Tcoor newPos;
  Tcoor pos;
  a->getPosition(pos);
  float currx=pos[0];
  float curry=pos[1];
  float currz=pos[2];

  newPos[0]= currx *matrix[0][0] + curry*matrix[0][1] + currz*matrix[0][2]+ matrix[0][3];
  newPos[1]= currx *matrix[1][0] + curry*matrix[1][1] + currz*matrix[1][2]+ matrix[1][3];
  newPos[2]= currx *matrix[2][0] + curry*matrix[2][1] + currz*matrix[2][2]+ matrix[2][3];

  a->setPosition(newPos);
  return true;
}

bool M4Rot::apply(Atom *a, Atom *b){

  Tcoor newPos;
  Tcoor pos;
  a->getPosition(pos);
  float currx=pos[0];
  float curry=pos[1];
  float currz=pos[2];

  newPos[0]= currx *matrix[0][0] + curry*matrix[0][1] + currz*matrix[0][2]+ matrix[0][3];
  newPos[1]= currx *matrix[1][0] + curry*matrix[1][1] + currz*matrix[1][2]+ matrix[1][3];
  newPos[2]= currx *matrix[2][0] + curry*matrix[2][1] + currz*matrix[2][2]+ matrix[2][3];

  b->setPosition(newPos);
  return true;
}
