#include <stdio.h>

#include "Molecule.h"

Molecule::Molecule(Molecule *old,bool with_elements)
{

  strcpy(id,old->getName());
  nid=old->getIdNumber();
  father=NULL;
  if(with_elements)
    {
	  limit=old->getLimit();
	  elements=(PDB_Contained**)malloc(sizeof(PDB_Contained*)*limit);
	  if(elements==NULL)
		  fprintf(stderr,"Error de asignacion de memoria\n");
      for(int i=0;i<limit;i++)
	  {
    	if(old->getClass()==pdb_protein || old->getClass()==pdb_nacid)
    	  elements[i]=(PDB_Contained*)new  Chain((Chain*)old->getE(i));
    	else
    	  elements[i]=(PDB_Contained*)new  Atom((Atom*)old->getE(i));

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

Molecule:: Molecule()
{
  elements=(PDB_Contained**)malloc(0);
  limit=0;
  currentE=0;
  nid=0;
  father=NULL;
}

Molecule::Molecule(Tname name, int in_nid)
{
  // strcpy(id,name);
  memcpy(id, name, strlen(name) + 1);

  elements=(PDB_Contained**)malloc(0);
  limit=0;
  currentE=0;
  nid=in_nid;
  father=NULL;
}

Molecule::~Molecule()
{

  for(int i=0;i<limit;i++)
  {
	  switch(elements[i]->getClass())
	  {
		  case pdb_atom:
			  delete (Atom*)elements[i];
			  break;
		  case pdb_chain:
			  delete (Chain*)elements[i];
			  break;
		  default:
		    break;
	  }

  }
  free(elements);
}




