#include "Macromolecule.h"
#include "pdbIter.h"
#include <stdio.h>


pdbIter::pdbIter()
{
  end_macro=0,
  end_chain=0,
  end_molecule=0,
  end_segment=0,
  end_fragment=0,
  end_atom=0;

  pos_macro=-1,
  pos_chain=-1,
  pos_molecule=-1,
  pos_segment=-1,
  pos_fragment=-1,
  pos_atom=-1;

  L_macro=(Macromolecule**)malloc(0),
  L_chain=(Chain**)malloc(0),
  L_molecule=(Molecule**)malloc(0),
  L_segment=(Segment**)malloc(0),
  L_fragment=(Fragment**)malloc(0),
  L_atom=(Atom**)malloc(0);

}

pdbIter::pdbIter(pdbIter *i2)
{
  /*L_macro=(Macromolecule**)malloc(sizeof(Macromolecule*)*i2->end_macro);
  L_chain=(Chain**)malloc(sizeof(Chain*)*i2->end_chain);
  L_molecule=(Molecule**)malloc(sizeof(Molecules*)*i2->end_molecule);
  L_segment=(Segment**)malloc(sizeof(Segment*)*i2->end_segment);
  L_fragment=(Fragment**)malloc(sizeof(Fragment*)*i2->end_fragment);
  L_atom=(Atom**)malloc(sizeof(Atom*)*i2->end_atom);

  memcpy(L_macro,i2->L_macro,sizeof(Macromolecule*)*i2->end_macro);
  memcpy(L_chain,i2->L_chain,sizeof(Chain*)*i2->end_chain);
  memcpy(L_molecule,i2->L_molecule,sizeof(Molecule*)*i2->end_molecule);
  memcpy(L_segment,i2->L_segment,sizeof(Segment*)*i2->end_segment);
  memcpy(L_fragment,i2->L_fragment,sizeof(Fragment*)*i2->end_fragment);
  memcpy(L_atom,i2->L_atom,sizeof(Atom*)*i2->end_atom);
  */

  L_macro=i2->L_macro;
  L_molecule=i2->L_molecule;
  L_chain=i2->L_chain;
  L_segment=i2->L_segment;
  L_fragment=i2->L_fragment;
  L_atom=i2->L_atom;

  end_macro=i2->end_macro,
  end_chain=i2->end_chain,
  end_molecule=i2->end_molecule,
  end_segment=i2->end_segment,
  end_fragment=i2->end_fragment,
  end_atom=i2->end_atom;

  if(i2->pos_macro>-1)
    pos_macro=0;
  if(i2->pos_molecule>-1)
    pos_molecule=0;
  if(i2->pos_chain>-1)
    pos_chain=0;
  if(i2->pos_segment>-1)
    pos_segment=0;
  if(i2->pos_fragment>-1)
    pos_fragment=0;
  if(i2->pos_atom>-1)
    pos_atom=0;
}

pdbIter::pdbIter(PDB_Container *c, bool include_smol,bool include_nacid, bool include_protein,bool in_special_iter)
{
   int i;
   int cont=0;

  end_macro=0,
  end_chain=0,
  end_molecule=0,
  end_segment=0,
  end_fragment=0,
  end_atom=0;

  pos_macro=-1,
  pos_chain=-1,
  pos_molecule=-1,
  pos_segment=-1,
  pos_fragment=-1,
  pos_atom=-1;
  special_iter=in_special_iter;

  L_macro=(Macromolecule**)malloc(0),
  L_chain=(Chain**)malloc(0),
  L_molecule=(Molecule**)malloc(0),
  L_segment=(Segment**)malloc(0),
  L_fragment=(Fragment**)malloc(0),
  L_atom=(Atom**)malloc(0);

  TElement pdb_class= c->getClass();
  TMOL moltype;
  Segment *virtual_segment;
  int end=c->getLimit();

  if(end<=0)
    return;

  c->initAll();
  pdbIter *son;

  do
  {
	  moltype= ((PDB_Contained*)(c->getCurrent()))->getMolType();
	  if( (moltype==tmol_protein && include_protein) ||  (moltype==tmol_nacid && include_nacid) || (moltype==tmol_smol && include_smol) )
		  cont++;
  } while(c->next()!=false);

  switch(pdb_class)
  {
    case pdb_world:
      end_macro=cont;
      L_macro=(Macromolecule**)realloc(L_macro,sizeof(Macromolecule*)*cont);
      pos_macro=0;
      break;

    case pdb_macro:
      end_molecule=cont;
      L_molecule=(Molecule**)realloc(L_molecule,sizeof(Molecule*)*cont);
      pos_molecule=0;
      break;

    case pdb_protein:
    case pdb_nacid:
      end_chain=cont;
      L_chain=(Chain**)realloc(L_chain,sizeof(Chain*)*cont);
      pos_chain=0;
      break;

    case pdb_chain:
      end_segment=cont;
      L_segment=(Segment**)realloc(L_segment,sizeof(Segment*)*cont);
      pos_segment=0;
      break;

    case pdb_segment:
       	if(c->getName()[0] != 'V')
    	{
    		end_fragment=cont;
    		L_fragment=(Fragment**)realloc(L_fragment,sizeof(Fragment*)*cont);
    	}
    	pos_fragment=0;
      break;

    case pdb_smol:
    	if(special_iter)
    	{
        	L_segment=(Segment**)realloc(L_segment,sizeof(Segment*)*(end_segment+1));
        	virtual_segment=new Segment((char *)"V");
    		virtual_segment->add(c);
    		L_segment[end_segment]=virtual_segment;
    		end_segment++;
    	}
    	if(!(L_fragment=(Fragment**)realloc(L_fragment,sizeof(Fragment*)*(end_fragment+1))))
    		fprintf(stderr,"L_fragment NOT allocated by realloc (end_fragment) %d)\n",end_fragment);
    	L_fragment[end_fragment]=(Fragment*)c;
    	end_fragment++;
    case pdb_residue:
    case pdb_nucleotide:
      end_atom=cont;
      L_atom=(Atom**)realloc(L_atom,sizeof(Atom*)*cont);
      pos_atom=0;
      break;
	default:
	  break;
  }

  c->initAll();
  i=0;
  do{
	  moltype= ((PDB_Contained*)(c->getCurrent()))->getMolType();
	 if( (moltype==tmol_protein && include_protein) ||  (moltype==tmol_nacid && include_nacid) || (moltype==tmol_smol && include_smol) )
	 {
			switch(pdb_class)
			{
			case pdb_world:
			  L_macro[i++]=(Macromolecule*)c->getCurrent();
			  break;

			case pdb_macro:
			  L_molecule[i++]=(Molecule*)c->getCurrent();
			  break;

			case pdb_protein:
			case pdb_nacid:
			  L_chain[i++]=(Chain*)c->getCurrent();
			  break;

			case pdb_chain:
				L_segment[i++]=(Segment*)c->getCurrent();
				break;

			case pdb_segment:
				if(c->getName()[0] != 'V') // Pablo fixed bug: L_fragment should only be assigned if it is a "non-virtual" segment
										   // Segment name is always "Segment" unless it be a Virtual_segment ('V')
					L_fragment[i++]=(Fragment*)c->getCurrent();
				break;

			case pdb_smol:
			case pdb_residue:
			case pdb_nucleotide:
				L_atom[i++]=(Atom*)c->getCurrent();
				break;
			default:
			  break;
			}
			if(pdb_class!=pdb_residue && pdb_class!=pdb_nucleotide && pdb_class!=pdb_smol)
			{
				son= new pdbIter((PDB_Container*)c->getCurrent(),include_smol,include_nacid, include_protein,in_special_iter);
				add(son);
				delete son;
			}
	 }

  }while(c->next()!=false);
  c->initAll();
}

void pdbIter::clean_virtual()
{
	int i;
	if(special_iter)
		for(i=0; i< end_segment;i++)
		{
			Segment *s=L_segment[i];
			if(s->getName()[0] == 'V')
			{
				s->removeAll();
				delete s;
			}
		}
}

pdbIter::~pdbIter()
{


  free(L_macro);
  free(L_chain);
  free(L_molecule);
  free(L_segment);
  free(L_fragment);
  free(L_atom);
}

void pdbIter::add(pdbIter *other)
{

 if(other->num_macro()>0)
 {
   L_macro=(Macromolecule**)realloc(L_macro,sizeof(Macromolecule*)*(end_macro+other->num_macro()));
   for(other->pos_macro=0;!other->gend_macro();other->next_macro())
    L_macro[end_macro++]=other->get_macro();
 }

 if(other->num_molecule()>0)
 {
   L_molecule=(Molecule**)realloc(L_molecule,sizeof(Molecule*)*(end_molecule+other->num_molecule()));
   for(other->pos_molecule=0;!other->gend_molecule();other->next_molecule())
    L_molecule[end_molecule++]=other->get_molecule();
 }

 if(other->num_chain()>0)
 {
  L_chain=(Chain**)realloc(L_chain,sizeof(Chain*)*(end_chain+other->num_chain()));
  for(other->pos_chain=0;!other->gend_chain();other->next_chain())
    L_chain[end_chain++]=other->get_chain();
 }

 if(other->num_segment()>0)
 {
  L_segment=(Segment**)realloc(L_segment,sizeof(Segment*)*(end_segment+other->num_segment()));
  for(other->pos_segment=0;!other->gend_segment();other->next_segment())
    L_segment[end_segment++]=other->get_segment();
 }

 if(other->num_fragment()>0)
 {
   L_fragment=(Fragment**)realloc(L_fragment,sizeof(Fragment*)*(end_fragment+other->num_fragment()));
   for(other->pos_fragment=0;!other->gend_fragment();other->next_fragment())
    L_fragment[end_fragment++]=other->get_fragment();
 }

 if(other->num_atom()>0)
 {
    L_atom=(Atom**)realloc(L_atom,sizeof(Atom*)*(end_atom+other->num_atom()));
    for(other->pos_atom=0;!other->gend_atom();other->next_atom())
      L_atom[end_atom++]=other->get_atom();
 }
}

bool pdbIter::operator ++ (int)
{
  return next_atom();

};

bool pdbIter::operator -- (int)
{
   return back_atom();
};

int pdbIter::num_macro()
{
 return end_macro;
}

int pdbIter::num_chain()
{
 return end_chain;
}

int pdbIter::num_molecule()
{
 return end_molecule;
}

int pdbIter::num_segment()
{
 return end_segment;
}

int pdbIter::num_fragment()
{
 return end_fragment;
}

int pdbIter::num_atom()
{
 return end_atom;
}

Macromolecule* pdbIter::get_macro()
{
 if(pos_macro>=0 && pos_macro<end_macro)
   return ((Macromolecule*)L_macro[pos_macro]);
 else
  return NULL;
}

Chain* pdbIter::get_chain()
{
 if(pos_chain>=0 && pos_chain<end_chain)
   return ((Chain*)L_chain[pos_chain]);
 else
   return NULL;
}

Molecule* pdbIter::get_molecule()
{
 if(pos_molecule>=0 && pos_molecule<end_molecule)
   return ((Molecule*)L_molecule[pos_molecule]);
 else
    return NULL;
}

Segment* pdbIter::get_segment()
{
 if(pos_segment>=0 && pos_segment<end_segment)
   return ((Segment*)L_segment[pos_segment]);
 else
    return NULL;
}

Fragment* pdbIter::get_fragment()
{
 if(pos_fragment>=0 && pos_fragment<end_fragment)
   return ((Fragment*)L_fragment[pos_fragment]);
 else
    return NULL;
}

Atom* pdbIter::get_atom()
{
 if(pos_atom>=0 && pos_atom<end_atom)
   return ((Atom*)L_atom[pos_atom]);
 else
    return NULL;

}

bool pdbIter::next_macro()
{
 if(pos_macro<end_macro)
 {
   pos_macro++;
   return true;
 }
 return false;
}

bool pdbIter::next_chain()
{
 if(pos_chain<end_chain)
 {
   pos_chain++;
   return true;
 }
 return false;
}

bool pdbIter::next_molecule()
{
 if(pos_molecule<end_molecule)
 {
   pos_molecule++;
   return true;
 }
 return false;
}

bool pdbIter::next_segment()
{
 if(pos_segment<end_segment)
 {
   pos_segment++;
   return true;
 }
 return false;
}

bool pdbIter::next_fragment()
{
 if(pos_fragment<end_fragment)
 {
   pos_fragment++;
   return true;
 }
 return false;
}

bool pdbIter::next_atom()
{
 if(pos_atom<end_atom)
 {
   pos_atom++;
   return true;
 }
 return false;
}

bool pdbIter::gend_macro()
{
 if(pos_macro<end_macro)
   return false;

 return true;
}

bool pdbIter::gend_chain()
{
 if(pos_chain<end_chain)
   return false;

 return true;
}

bool pdbIter::gend_molecule()
{
 if(pos_molecule<end_molecule)
   return false;

 return true;
}

bool pdbIter::gend_segment()
{
 if(pos_segment<end_segment)
   return false;

 return true;
}

bool pdbIter::gend_fragment()
{
 if(pos_fragment<end_fragment)
   return false;

 return true;
}

bool pdbIter::gend_atom()
{
 if(pos_atom<end_atom)
   return false;

 return true;
}



bool pdbIter::ginit_macro()
{
 if(pos_macro<0)
   return true;

 return false;
}
bool pdbIter::ginit_chain()
{
 if(pos_chain<0)
   return true;

 return false;
}
bool pdbIter::ginit_molecule()
{
 if(pos_molecule<0)
 return true;

return false;
}
bool pdbIter::ginit_segment()
{
 if(pos_segment<0)
   return true;

 return false;
}
bool pdbIter::ginit_fragment()
{
 if(pos_fragment<0)
   return true;

 return false;
}
bool pdbIter::ginit_atom()
{
 if(pos_atom<0)
   return true;

 return false;
}

bool pdbIter::back_macro()
{
 if(pos_macro>=0)
 {
   pos_macro--;
   return true;
 }
 return false;
}

bool pdbIter::back_chain()
{
 if(pos_chain>=0)
 {
   pos_chain--;
   return true;
 }
 return false;
}

bool pdbIter::back_molecule()
{
 if(pos_molecule>=0)
 {
   pos_molecule--;
   return true;
 }
 return false;
}

bool pdbIter::back_segment()
{
 if(pos_segment>=0)
 {
   pos_segment--;
   return true;
 }
 return false;
}

bool pdbIter::back_fragment()
{
 if(pos_fragment>=0)
 {
   pos_fragment--;
   return true;
 }
 return false;
}

bool pdbIter::back_atom()
{
 if(pos_atom>=0)
 {
   pos_atom--;
   return true;
 }
 return false;
}




