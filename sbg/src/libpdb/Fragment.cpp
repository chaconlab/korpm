
#include <stdio.h>

#include "Fragment.h"
#include "ResIni.h"

using namespace std;
Fragment:: Fragment()
{
  strcpy(id,"Fragment");

  elements=(PDB_Contained**)malloc(0);
  limit=0;
  currentE=0;
  nid=0;
  npos=0;
  father=NULL;
}

Fragment:: Fragment(Tname name, int i_nid, int i_npos,char i_letter)
{
  // strcpy(id,name);
  memcpy(id, name, strlen(name) + 1);

  fragid = resnum_from_resname(id); // Mon added...
  if(fragid >= 24 && fragid <= 26)
	  fragid = 6; // HIP, HID, or HIE --> HIS.
  else if(fragid >= 21 && fragid <= 22)
	  fragid = 1; // CYX or CYM --> CYS
  else if(fragid == 29)
	  fragid = 10; // MSE --> MET
  else if(fragid == 20)
	  fragid = 2; // ASH --> ASP
  else if(fragid == 23)
	  fragid = 3; // GLH --> GLU
  else if(fragid == 28)
	  fragid = 19; // TYM --> TYR
  elements=(PDB_Contained**)malloc(0);
  limit=0;
  currentE=0;
  nid=i_nid;
  npos=i_npos;
  letter=i_letter;
  father=NULL;

}
Fragment:: Fragment(Fragment *old, bool with_elements)
{
  // strcpy(id,old->getName());
  memcpy(id, old->getName(), strlen(old->getName()) + 1);

  fragid = resnum_from_resname(id); // Mon added...
  if(fragid >= 24 && fragid <= 26)
	  fragid = 6; // HIP, HID, or HIE --> HIS.
  else if(fragid >= 21 && fragid <= 22)
	  fragid = 1; // CYX or CYM --> CYS
  else if(fragid == 29)
	  fragid = 10; // MSE --> MET
  else if(fragid == 20)
	  fragid = 2; // ASH --> ASP
  else if(fragid == 23)
	  fragid = 3; // GLH --> GLU
  else if(fragid == 28)
	  fragid = 19; // TYM --> TYR
  nid=old->getIdNumber();
  npos=old->get_pos();
  letter=old->get_letter();
  if(with_elements)
  {
	  limit=old->getLimit();
	  elements=(PDB_Contained**)malloc(sizeof(PDB_Contained*)*limit);
	  if(elements==NULL)
		  fprintf(stdout,"Error in memory allocation\n");

	  for(int i=0;i<limit;i++)
	  {
		  elements[i]=new Atom((Atom*)old->getE(i));
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


Fragment:: ~Fragment()
{

  for(int i=0;i<limit;i++)
  {
	  delete (Atom*)elements[i];
  }

  free(elements);

}


int Fragment::get_pos()
{
	return npos;
}

char Fragment::get_letter()
{
	return letter;
}


char Fragment::get_ss()
{
	return ss;
}

void Fragment::set_ss(char ch)
{
	ss=ch;
}
