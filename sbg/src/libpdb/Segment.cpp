
#include <stdio.h>

#include "Segment.h"

using namespace std;

Segment:: Segment()
{
  strcpy(id,"Segment");
  elements=(PDB_Contained**)malloc(0);
  limit=0;
  currentE=0;
  nid=0;
  father=NULL;
}

Segment:: Segment(Tname name, int in_nid )
{
  // strcpy(id,name);
  memcpy(id, name, strlen(name) + 1);
  elements=(PDB_Contained**)malloc(0);
  limit=0;
  currentE=0;
  nid=in_nid;
  father=NULL;

}

Segment::Segment(Segment *old,bool with_elements)
{
  //strcpy(id,old->getName());
  memcpy(id, old->getName(), strlen(old->getName()) + 1);

  father=NULL;
  nid=old->getIdNumber();

  if(with_elements)
  {
	  limit=old->getLimit();
	  elements=(PDB_Contained**)malloc(sizeof(PDB_Contained*)*limit);
	  if(elements==NULL)
		  fprintf(stdout,"Memory allocation error\n");
	  for(int i=0;i<limit;i++)
	  {
		  if(old->getMolType()==tmol_protein)
			  elements[i]=new  Residue((Residue*)old->getE(i));
		  else
			  elements[i]=new  Nucleotide((Nucleotide*)old->getE(i));

		  ((PDB_Contained*)(elements[i]))->setFather(this);
	  }
  }
  else
  {
	  limit=0;
	  elements=(PDB_Contained**)malloc(0);
  }

  if(limit!=0)
    currentE=0;
  else
    currentE=-1;
}

Segment:: ~Segment()
{
  int i;

  for(i=0;i<limit;i++)
  {
	  delete (Fragment*)elements[i];
  }
  free(elements);

}

TElement Segment::getClass()
{
  return pdb_segment;
}

TMOL Segment::getMolType()
{
	if(limit>0)
	{
		return elements[0]->getMolType();
	}
	else
	{
		PDB_Container *f;
		f=(PDB_Container*)getFather();
		if(f!=NULL)
			return f->getMolType();
		else
			return tmol_null;
	}
}
