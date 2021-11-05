/*
 * NAcid.cpp
 *
 *  Created on: Feb 4, 2010
 *      Author: nacho
 */

#include <stdio.h>

#include "NAcid.h"

NAcid::NAcid(Tname name,int i_nid):Molecule(name,i_nid)
{
}
NAcid::NAcid(NAcid *old, bool with_elements):Molecule(old,with_elements)
{
	}

NAcid::NAcid():Molecule()
{
 strcpy(id,"NucleAcid");
}




TElement NAcid::getClass()
{
  return pdb_nacid;
}
TMOL NAcid::getMolType()
{
	return tmol_nacid;
}

