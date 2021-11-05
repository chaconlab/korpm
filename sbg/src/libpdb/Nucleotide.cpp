/*
 * Nucleotide.cpp
 *
 *  Created on: Feb 5, 2010
 *      Author: nacho
 */

#include <Nucleotide.h>

Nucleotide::Nucleotide( Tname name, int in_nid, int in_npos, char i_letter ):Fragment( name, in_nid, in_npos, i_letter )
	  {
	  };

Nucleotide::Nucleotide( Nucleotide * old ,bool with_elements):Fragment( old ,with_elements)
	  {
	  };

TMOL Nucleotide::getMolType()
{
	return tmol_nacid;
}

TElement Nucleotide::getClass()
{
  return pdb_nucleotide;
}
