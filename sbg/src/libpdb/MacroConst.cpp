#include "Macromolecule.h"
#include <stdio.h>

using namespace std;

Macromolecule::Macromolecule()
{
  strcpy( id, (const char *)"Macromolecule" );
  elements = ( PDB_Contained * * ) malloc( 0 );
  limit = 0;
  currentE = 0;

  nid=0;
  father=NULL;


}

Macromolecule::Macromolecule( Tname name, int in_nid)
{
  //strcpy( id, (const char *)name );
  memcpy(id, name, strlen(name) + 1);

  elements = ( PDB_Contained * * ) malloc( 0 );
  limit = 0;
  currentE = 0;
  nid = in_nid;
  father=NULL;
}

Macromolecule::Macromolecule( Macromolecule * old )
{
  int i;
  Molecule * m;

  // strcpy( id, (const char *)old->getName() );
  memcpy(id, old->getName(), strlen(old->getName()) + 1);

  limit = old->getLimit();
  nid=old->getIdNumber();
  father=NULL;
  elements = ( PDB_Contained * * ) malloc( sizeof( Molecule * ) * limit );
  if ( elements == NULL )
    fprintf(stdout,"Error in memory allocation\n");

  for ( i = 0; i < limit; i++ )
  {
    m = ( Molecule * ) old->getE( i );
    switch (m->getClass())
	{
    case pdb_protein:
    	//printf("proteina\n");
    	elements[i] = new Protein( ( Protein * ) m );
    	break;
    case pdb_nacid:
    	//printf("nacid\n");
    	elements[i] = new NAcid( ( NAcid * ) m );
    	break;
    case pdb_smol:
    	//printf("smol\n");
    	elements[i] = new SMol( ( SMol * ) m );
    	break;
    default:
      	fprintf(stdout,"Error molecule type not found\n");
	}

    ((PDB_Contained*)(elements[i]))->setFather(this);
  }

  if ( limit != 0 )
    currentE = 0;
  else
    currentE = -1;




  macromoleculeDimension.geometricCenter[0] = old->getDimension()->geometricCenter[0];
  macromoleculeDimension.geometricCenter[1] = old->getDimension()->geometricCenter[1];
  macromoleculeDimension.geometricCenter[2] = old->getDimension()->geometricCenter[2];
  //geoBox();
}

Macromolecule::Macromolecule( Macromolecule **list, int num_m )
{
  int i,j,k,k2;
  Molecule * m;

  strcpy( id, (const char *)list[0]->getName() );
  nid=list[0]->getIdNumber();
  father=NULL;

  limit=0;
  for(i=0;i<num_m;i++)
    limit += list[i]->getLimit();
  elements = ( PDB_Contained * * ) malloc( sizeof( Molecule * ) * limit );
  if ( elements == NULL )
    fprintf(stdout,"Error in memory allocation\n");

  k2=0;
  for(j=0;j<num_m;j++)
  {
    k=list[j]->getLimit();
    for ( i = 0; i < k; i++ )
    {
      m = ( Molecule * ) list[j]->getE( i );

      switch (m->getClass())
      	{
          case pdb_protein:
          	elements[k2] = new Protein( ( Protein * ) m );
          	break;
          case pdb_nacid:
          	elements[k2] = new NAcid( ( NAcid * ) m );
          	break;
          case pdb_smol:
          	elements[k2] = new SMol( ( SMol * ) m );
          	break;
          default:
          	fprintf(stdout,"Error molecule type not found\n");
      	}

      ((PDB_Contained*)(elements[k2++]))->setFather(this);
    }
  }

  if ( limit != 0 )
    currentE = 0;
  else
    currentE = -1;



  geoBox();
}


Macromolecule::~Macromolecule()
{
  int i;

  for ( i = 0; i < limit; i++ )
  {

    //( ( Protein * ) elements[i] )->destructor();
    delete (Protein*)elements[i];
  }
  free( elements );




  macromoleculeDimension.~Dimensions();
}

TElement Macromolecule::getClass()
{
  return pdb_macro;
}

TMOL Macromolecule::getMolType()
{
	return tmol_null;
}
