#include <stdio.h>

#include "Protein.h"

Protein::Protein(Tname name,int i_nid):Molecule(name,i_nid)
{
}
Protein::Protein(Protein *old, bool with_elements):Molecule(old,with_elements)
{
}

Protein:: Protein():Molecule()
{
 strcpy(id,"Protein");
}




TElement Protein::getClass()
{
  return pdb_protein;
}

TMOL Protein::getMolType()
{
	return tmol_protein;
}

