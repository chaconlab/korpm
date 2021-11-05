
#include <stdio.h>

#include "Chain.h"

using namespace std;


Chain:: Chain()
{
  strcpy(id,"Chain");
  elements=(PDB_Contained**)malloc(0);
  limit=0;
  currentE=0;
  nid=0;
  father=NULL;

}

Chain:: Chain(Tname name,int i_nid )
{
  strcpy(id,name);
  elements=(PDB_Contained**)malloc(0);
  limit=0;
  currentE=0;
  nid=i_nid;
  father=NULL;


}

Chain::Chain(Chain *old, bool with_elements)
{
  int i;
  strcpy(id,old->getName());
  nid=old->getIdNumber();
  father=NULL;

  if(with_elements)
  {
	  limit=old->getLimit();
	  elements=(PDB_Contained**)malloc(sizeof(PDB_Contained*)*limit);
	  if(elements==NULL)
		  fprintf(stdout,"Error in memory allocation\n");

	  for(i=0;i<limit;i++)
	  {
		  elements[i]=new  Segment((Segment*)old->getE(i));
		  ((PDB_Contained*)(elements[i]))->setFather(this);
	  }
  }
  else
  {
	  elements=(PDB_Contained**)malloc(0);
	  limit=0;
  }

  if(limit!=0)
    currentE=0;
  else
    currentE=-1;
}

Chain:: ~Chain()
{
  int i;

  for(i=0;i<limit;i++)
  {
	delete (Segment*)elements[i];
  }
  free(elements);
}

TElement Chain::getClass()
{
  return pdb_chain;
}

TMOL Chain::getMolType()
{
	PDB_Container *f;
	f=(PDB_Container*)getFather();
	if(f!=NULL)
		return f->getMolType();
	else
		if(limit>0)
		{
			return elements[0]->getMolType();
		}
		else
		{
			return tmol_null;
		}
}

