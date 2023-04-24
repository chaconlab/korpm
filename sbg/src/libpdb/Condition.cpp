/***************************************************************************
                          Condition.cpp  -  description
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
#include <stdio.h>

#include <stdlib.h>
#include <string.h>

#include "Condition.h"



Condition::Condition(int lP,int hP,int lS,int hS,int lCh,int hCh,int lR,int hR,int lA,int hA)
{
  lowP=lP;
  highP=hP;
  lowS=lS;
  highS=hS;
  lowCh=lCh;
  highCh=hCh;
  lowR=lR;
  highR=hR;
  lowA=lA;
  highA=hA;
  atoms=(char**)malloc(sizeof(char*)*0);
  limit=0;
}

Condition::~Condition()
{
  for(int i=0;i<limit;i++)
    free(atoms[i]);
  free(atoms);
}

bool Condition::isValidP(int pos)
{
  if(lowP==-1)
  {
    if(highP==-1)
      return true;
    else
      if(highP>=pos)
        return true;
      else
        return false;
  }
  else
  {
    if(highP==-1)
      if(lowP<=pos)
        return true;
      else
        return false;
    else
      if(highP>=pos && lowP<=pos)
        return true;
      else
        return false;
  }
}

bool Condition::isValidS(int pos)
{
  if(lowS==-1)
  {
    if(highS==-1)
      return true;
    else
      if(highS>=pos)
        return true;
      else
        return false;
  }
  else
  {
    if(highS==-1)
      if(lowS<=pos)
        return true;
      else
        return false;
    else
      if(highS>=pos && lowS<=pos)
        return true;
      else
        return false;
  }
}

bool Condition::isValidCh(int pos)
{
  if(lowCh==-1)
  {
    if(highCh==-1)
      return true;
    else
      if(highCh>=pos)
        return true;
      else
        return false;
  }
  else
  {
    if(highCh==-1)
      if(lowCh<=pos)
        return true;
      else
        return false;
    else
      if(highCh>=pos && lowCh<=pos)
        return true;
      else
        return false;
  }
}

bool Condition::isValidR(int pos)
{
  if(lowR==-1)
  {
    if(highR==-1)
      return true;
    else
      if(highR>=pos)
        return true;
      else
        return false;
  }
  else
  {
    if(highR==-1)
      if(lowR<=pos)
        return true;
      else
        return false;
    else
      if(highR>=pos && lowR<=pos)
        return true;
      else
        return false;
  }
}


bool Condition::isValidA(int pos)
{
  if(lowA==-1)
  {
    if(highA==-1)
      return true;
    else
      if(highA>=pos)
        return true;
      else
        return false;
  }
  else
  {
    if(highA==-1)
      if(lowA<=pos)
        return true;
      else
        return false;
    else
      if(highA>=pos && lowA<=pos)
        return true;
      else
        return false;
  }
}

void Condition::add(char *cad)
{
  atoms=(char**)realloc(atoms,sizeof(char*)*limit + sizeof(char*));
  atoms[limit]=(char*)malloc(sizeof(char)*6);
  strcpy(atoms[limit],cad);
  limit++;
}

bool Condition::fitCondition(int posP,int posS,int posCh, int posR, int posA, char *cad)
{
  bool out=false;
  if( isValidP(posP) && isValidS(posS) && isValidCh(posCh) && isValidR(posR) && isValidA(posA))
  {
	 if(limit<=0)
		 return true;

    for(int i=0;i<limit && out==false ;i++)
    {
      //fprintf(stderr,"comparando [%s] con [%s]\n",cad,atoms[i]);
      if(strcmp(cad,atoms[i])==0)
      {
          out=true;
        //cout <<"Bingo"<<endl;
      }
    }
  }

  return out;
}


Conditions::Conditions()
{
  list=(Condition**)malloc(sizeof(Condition*)*0);
  limit=0;
}
Conditions::~Conditions()
{
  for(int i=0;i<limit;i++)
  {
  		//delete(list[i]);
  	  free(list[i]);
  }
  free(list);
}
void Conditions::add(Condition *c)
{
 list=(Condition**)realloc(list,sizeof(Condition*)*limit + sizeof(Condition*));
 list[limit]=c;
  limit++;
}

Condition * Conditions::getCondition(int i)
{
  return list[i];
}
