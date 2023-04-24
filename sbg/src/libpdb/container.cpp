/***************************************************************************
                          container.cpp  -  description
                             -------------------
    begin                : Wed Jul 21 2004
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

#include "container.h"
#include <stdio.h>


PDB_Contained::~PDB_Contained()
{

}



/*
PDB_Container::~PDB_Container()
{
  int i;
  for(i=0;i<limit;i++)
  {
	delete elements[i];
  }
  free(elements);
  limit=0;
}
*/


void *PDB_Contained::getFather()
{
	return father;
}

bool PDB_Contained::setFather(void *f)
{
    father=f;
	return true;
}

int PDB_Container::getPosition(Tname name)
{
  int i;

  for(i=0;i<(limit);i++)
  {
    if(strcmp(name,elements[i]->getName())==0)
      return i;
  }
  return -1;
}

bool PDB_Container::add(PDB_Contained *newE)
{
 elements=(PDB_Contained**)realloc(elements,(sizeof(PDB_Contained*)*limit)+(sizeof(PDB_Contained*)));
  if(elements==NULL)
    return false;
  elements[limit]=newE;
  limit++;
  if(newE->getFather()==NULL)
  {
	  newE->setFather(this);
  }
  return true;
}

void* PDB_Container::getCurrent()
{
  if(elements==NULL)
     return NULL;
  if(currentE!=-1)
    return elements[currentE];
  else
    return NULL;
}

int PDB_Container::end()
{
  if(limit!=0)
    currentE=limit-1;
  else
    currentE=-1;
  return currentE;
}

bool PDB_Container::erase(int pos)
{
  if((pos>=0)&&(pos<limit))
  {
    delete elements[pos];
//	  free(elements[pos]);
    return remove(pos);
  }
  return false;
}

bool PDB_Container::erase(Tname name)
{
  int index;
  index= getPosition(name);

  if (index!=-1)
    return erase(index);
  else
    return false;
}

bool PDB_Container::eraseAll()
{
  int i;
  for(i=0;i<limit;i++)
    delete elements[i];
//	  free(elements[i]);
   elements=(PDB_Contained**)realloc(elements,0);
  limit=0;
  currentE=-1;
  return true;
}




int PDB_Container::getLimit()
{
  return limit;
}

PDB_Contained* PDB_Container::getE(int pos)
{
  if((pos>=0)&&(pos<limit))
    return elements[pos];
  else
    return NULL;
}

bool PDB_Container::init()
{
  if (limit!=0)
  {
    currentE=0;
    return true;
  }
  else
  {
    currentE=-1;
    return false;
  }
}

bool PDB_Container::initAll()
{
  int i;

  for(i=0;i<limit;i++)
  {
      elements[i]->initAll();
  }
  return this->init();

}

bool PDB_Container::moveAll(Tcoor offset)
{
  int i;
  for(i=0;i<limit;i++)
      if(elements[i]->moveAll(offset)==false)
        return false;

  return true;
}

bool PDB_Container::moveAlln(Tcoor offset)
{
  int i;
  for(i=0;i<limit;i++)

      if(elements[i]->moveAlln(offset)==false)
        return false;

  return true;
}

bool PDB_Container::next()
{
 if(currentE<(limit-1))
  {
    currentE++;
    return true;
  }
  return false;
}

bool PDB_Container::previous()
{
 if(currentE>0)
  {
    currentE--;
    return true;
  }
  return false;
}

bool PDB_Container::remove(int pos)
{
 int i;

  if((pos>=0)&&(pos<limit))
  {
    for(i=pos;i<(limit-1);i++)
     elements[i]=elements[i+1];


   elements=(PDB_Contained**)realloc(elements,(sizeof(PDB_Contained*)*limit)-(sizeof(PDB_Contained*)));


   if(currentE>pos)
     currentE--;

   limit--;

   if(currentE>=limit)
     currentE=-1;

   if(elements==NULL)
     return false;
   return true;
  }
  return false;
}

bool PDB_Container::remove(Tname name)
{
  int index;
  index= getPosition(name);

  if (index!=-1)
    return remove(index);
  else
    return false;
}

bool PDB_Container::removeAll()
{
  elements=(PDB_Contained**)realloc(elements,0);
  limit=0;
  currentE=-1;
  return true;
}

// MON: This function leads to memory leakage... It may be obsolete??? Detected in Macromolecule::secondary_structure
// MON: We've checked that "erase_level" is not used across the whole library... Remove!
bool PDB_Container::erase_level(TElement level)
{
	int i;
	PDB_Contained *son;
	if(this->getClass()==level)
	{
		fprintf(stderr,"Error: erase_level function. Level specified is the same as the object.\n ");
		return false;
	}
	else
	{
		for(i=0;i<limit;i++)
		{
			son=elements[i];
			if(son->getClass()!=level)
			{
				if(son->getClass()==pdb_atom)
				{
					fprintf(stderr,"Alert: erase_level function. Atom level reached and erased, are you sure you introduced a correct level?.\n ");
				}
				else
				{
					((PDB_Container*)son)->erase_level(level);
				}
				delete son;
			}
		}
		elements=(PDB_Contained**)realloc(elements,0);
		limit=0;
		currentE=-1;
	}
	return true;
}

void PDB_Container::modifyName(char *name)
{
	strcpy(id,name);
}

char* PDB_Container::getName()
{
  return id;
}

void PDB_Container::exchange(int i, int j)
{
  PDB_Contained *tmp;
  tmp=elements[i];
  elements[i]=elements[j];
  elements[j]=tmp;
}

int PDB_Container::getIdNumber()
{
	return nid;
}


bool PDB_Container::check_consistency()
{
	int i;
	PDB_Container *father, *grandfather,*grandgrandfather;

	printf("OBJECT TYPE NUMBERS:\n");
	printf("pdb_world:\t0\n");
	printf("pdb_macro:\t1\n");
	printf("pdb_protein:\t2\n");
	printf("pdb_nacid:\t3\n");
	printf("pdb_smol:\t4\n");
	printf("pdb_chain:\t5\n");
	printf("pdb_segment:\t6\n");
	printf("pdb_residue:\t7\n");
	printf("pdb_nucleotide\t8\n");
	printf("pdb_atom:\t9\n");
	printf("pdb_hetero:\t10\n");
	printf("\n");

	printf("MOLECULE TYPE NUMBERS:\n");
	printf("tmol_null:\t0\n");
	printf("tmol_protein:\t1\n");
	printf("tmol_nacid:\t2\n");
	printf("tmol_smol:\t3\n");
	printf("\n");

	//ALL THIS BLOCK IS REDUNDANT CHECKING
	if(getClass()==pdb_residue || getClass()==pdb_nacid )
	{
		father=(PDB_Container*)getFather();
		grandfather=(PDB_Container*)father->getFather();
		grandgrandfather=(PDB_Container*)grandfather->getFather();


		if(getClass()==pdb_residue)
		{

			if(grandgrandfather->getClass()!=pdb_protein || grandgrandfather->getMolType()!=tmol_protein)
			{
				fprintf(stderr,"Consistency Error: Residue %s with id %d is not in a Protein but in a Molecule of type %d\n",
					getName(),getIdNumber(),grandgrandfather->getMolType());
				return false;
			}
		}
		if(getClass()==pdb_nacid)
		{
			if(grandgrandfather->getClass()!=pdb_nacid || grandgrandfather->getMolType()!=tmol_nacid)
			{
				fprintf(stderr,"Consistency Error: Nucleotide %s with id %d is not in a Nucleacid but in a Molecule of type %d\n",
				getName(),getIdNumber(),grandgrandfather->getMolType());
				return false;
			}
		}
	}

	 for(i=0;i<limit;i++)
		{
			 if(elements[i]->getMolType()!=getMolType() && this->getClass()!=pdb_macro)
			 {
				 fprintf(stderr,"Consistency Error: Element %s with identifier %d of object type %d is in Molecule of type %d while its son number %d of object type %d is in a Molecule of type %d\n",
							           getName(),getIdNumber(),getClass(),getMolType(),i,elements[i]->getClass(),elements[i]->getMolType());
				 return false;
			 }
			 else
			 {
				if(getClass()!=pdb_residue && getClass()!=pdb_nacid )
					if( !( ((PDB_Container*)elements[i])->check_consistency() ) )
						return false;

			 }
		 }

	 return true;
}

int PDB_Container::num_atoms()
{
	int i, cont=0;

	if(getClass()==pdb_residue || getClass()==pdb_nucleotide || getClass()==pdb_smol)
		return limit;
	else
	{
		for(i=0;i<limit;i++)
		{
			cont+=((PDB_Container*)elements[i])->num_atoms();
		}
		return cont;
	}
}

void PDB_Container::delete_hydrogens()
{
	int i;
	bool *list_erase;
	Element *elem;

	if(getClass()==pdb_residue || getClass()==pdb_nucleotide || getClass()==pdb_smol)
	{
		list_erase=(bool*)malloc(sizeof(bool)*limit);
		for(i=0;i<limit;i++)
		{
			//fprintf(stderr,"atom: %s\n",((Atom*)elements[i])->getName());
//			if ( (((Atom*)elements[i])->getElement())->symbol==H )
			elem = ((Atom*)elements[i])->getElement();
			if ( elem->symbol==H || elem->symbol==D )
			{
				list_erase[i]=true;
			}
			else
			{
				list_erase[i]=false;
			}
		}
		for(i=limit-1;i>=0;i--)
		{
			if(list_erase[i])
			{
				erase(i);
			}
		}

		free(list_erase);
	}
	else // MON: check this (is recursively senseless...) Consider deletion below lines...
	{
		for(i=0;i<limit;i++)
		{
			((PDB_Container*)elements[i])->delete_hydrogens();
		}
	}
}

void PDB_Container::delete_hetero()
{
	int i;
	bool *list_erase;

	list_erase=(bool*)malloc(sizeof(bool)*limit);

	//fprintf(stderr,"atom-->: %s\n", getName());
	for(i=0;i<limit;i++)  {

		if ( ((PDB_Container*)elements[i])->getClass()==pdb_smol )
		{
			//fprintf(stderr,"atom--> %d %s %d\n",i, ((PDB_Container*)elements[i])->getName(), ((PDB_Container*)elements[i])->getClass());
			list_erase[i]=true;
		} else
		{
			list_erase[i]=false;
		}
	}

	for(i=limit-1;i>=0;i--)
		if(list_erase[i])
			erase(i);

	free(list_erase);
}

void PDB_Container::delete_water()
{
	int i;
	bool *list_erase;

	list_erase=(bool*)malloc(sizeof(bool)*limit);

	//fprintf(stderr,"atom-->: %s\n", getName());
	for(i=0;i<limit;i++)  {

		if (( strcmp(((PDB_Container*)elements[i])->getName(), "HOH")==0 ) ||
				( strcmp(((PDB_Container*)elements[i])->getName(), "DOD")==0 ) ||
				( strcmp(((PDB_Container*)elements[i])->getName(), "WAT")==0 ) ||
				( strcmp(((PDB_Container*)elements[i])->getName(), "TIP")==0 ))


		{
			//fprintf(stderr,"atom--> %d %s %d\n",i, ((PDB_Container*)elements[i])->getName(), ((PDB_Container*)elements[i])->getClass());
			list_erase[i]=true;
		} else
		{
			list_erase[i]=false;
		}
	}

	for(i=limit-1;i>=0;i--)
		if(list_erase[i])
			erase(i);

	free(list_erase);
}

// Delete contained duplicate atoms preserving the first occurrence (Mon 21/02/2018).
void PDB_Container::delete_duplicates()
{
	int i,j;
	bool *list_erase;

	list_erase = (bool*) malloc(sizeof(bool)*limit);
	list_erase[0] = false;

	// fprintf(stderr,"\n%s %d\n",this->getName(),this->getIdNumber());

	for(i=0;i<limit;i++)
		for(j=i+1;j<limit;j++)
		{
			// fprintf(stderr,"i is %s and j is %s\n", ((PDB_Container*)elements[i])->getName(), ((PDB_Container*)elements[j])->getName() );

			// If contained elements have the same name one of them should be deleted... (deleting the j-th occurrence)
			if( strcmp( ((PDB_Container*)elements[i])->getName(), ((PDB_Container*)elements[j])->getName() ) == 0 )
			{
				list_erase[j]=true;
				// fprintf(stderr,"%s is the same as %s\n", ((PDB_Container*)elements[i])->getName(), ((PDB_Container*)elements[j])->getName() );
				fprintf(stderr,"Warning duplicate atoms within the same residue found! Residue: %s %d\n",this->getName(),this->getIdNumber());
				// exit(1);
			}
			else
				list_erase[j]=false;
		}

	for(i=limit-1;i>=0;i--)
		if(list_erase[i])
			erase(i);

	free(list_erase);
}



