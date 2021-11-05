/*
 * SMol.cpp
 *
 *  Created on: Feb 4, 2010
 *      Author: nacho
 */
#include <stdio.h>

#include "SMol.h"

///Auxiliar function to find the position of an atom in a Smol
int return_pos_aux_smol(SMol *smol,Atom *at);

SMol::SMol(Tname name,int i_nid,int i_pos):Molecule(name,i_nid)
{
	strcpy(id,name);
	npos=i_pos;
	bonds=NULL;
	num_bonds=0;
}
SMol::SMol(SMol *old, bool with_elements):Molecule(old,with_elements)
{
	bonds=NULL;
	num_bonds=0;

	//INTRODUCE BONDS
	int i,pos_init,pos_final;
	Bond *bond,*new_bond;
	Atom *at_init,*at_final;

	//Find all bonds in old Smol and introduce in the new Smol
	for(i=0;i<old->get_numBonds();i++)
	{
		//Bond in old Smol
		//fprintf(stderr,"\n%d\n",i);
		//fprintf(stderr,"Bond in old Smol\n");
		bond=old->getBond(i);

		//Positions of the bonded atoms in old SMol
		//fprintf(stderr,"Positions of the bonded atoms\n");
		pos_init=return_pos_aux_smol(old,bond->getInit());
		pos_final=return_pos_aux_smol(old,bond->getFinal());

		//Corresponding Atoms bonded in new SMol
		//fprintf(stderr,"Atoms bonded (%d %d)\n",pos_init,pos_final);
		at_init=(Atom*)(this->getE(pos_init));
		at_final=(Atom*)(this->getE(pos_final));

		//Creation of new Bond
		//fprintf(stderr,"Creation of new Bond\n");
		new_bond=new Bond(at_init,at_final,bond->getLink());

		//Link the atoms to the bond
		//fprintf(stderr,"Link the atoms to the bond (%s %s)\n",at_init->getPdbName(),at_final->getPdbName());
		at_init->insertBond(new_bond);
		at_final->insertBond(new_bond);

		//Introduction of the bond in Smol
		//fprintf(stderr,"Introduction of the bond in Smol\n");
		this->insertBond(new_bond);
	}

}

SMol::SMol():Molecule()
{
 strcpy(id,"SMolecule");
 bonds=NULL;
 num_bonds=0;
}

SMol::~SMol()
{
	int i;
	if(num_bonds>0)
	{
		for(i=0;i<num_bonds;i++)
			delete bonds[i];
		free(bonds);
	}
}

int SMol::get_pos()
{
	return npos;
}


TElement SMol::getClass()
{
  return pdb_smol;
}
TMOL SMol::getMolType()
{
	return tmol_smol;
}

void SMol::insertBond(Bond *b)
{
	num_bonds++;

	bonds=(Bond**)realloc(bonds,sizeof(Bond*)*num_bonds);
	bonds[num_bonds-1]=b;
}

bool SMol::deleteBond(Bond *b)
{
	int i,j;

	for(i=0;i<num_bonds;i++)
	{
		if(bonds[i]==b)
		{
			(b->getInit())->removeBond(b);
			(b->getFinal())->removeBond(b);

			for(j=i;j<num_bonds-1;j++)
			{
				bonds[j]=bonds[j+1];
			}
			num_bonds--;
			bonds=(Bond**)realloc(bonds,sizeof(Bond*)*num_bonds);
			delete b;
			return true;
		}
	}

	return false;
}

void SMol::deleteHydrogenBonds()
{
	int i;
	Bond *b;

	for(i=num_bonds-1;i>=0;i--)
	{
		b=getBond(i);

//		if( ((b->getFinal())->getElement())->symbol==H || ((b->getInit())->getElement())->symbol==H )
		if( ((b->getFinal())->getElement())->symbol==H || ((b->getInit())->getElement())->symbol==H || ((b->getFinal())->getElement())->symbol==D || ((b->getInit())->getElement())->symbol==D )
			deleteBond(b);
	}
}

Bond *SMol::getBond(int i)
{
	if( i>num_bonds )
			return NULL;
	else
	{
		return bonds[i];
	}
}

int SMol::get_numBonds()
{
	return num_bonds;
}

int return_pos_aux_smol(SMol *smol,Atom *at)
{

	int i;
	Atom *at2;
	for(i=0;smol->getLimit();i++)
	{
		at2=(Atom*)smol->getE(i);
		if (at2==at)
			return i;
	}
	return -1;
}



