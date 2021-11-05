
#include "Macromolecule.h"

//Funciones para HeteroAtomos
bool Macromolecule::addH( Contained * newH )
{
  heteros = ( HeteroMolecule * * ) realloc( heteros, ( sizeof( Contained * ) * limitH ) + ( sizeof( Contained * ) ) );
  if ( heteros == NULL )
    return false;
  heteros[limitH] = ( HeteroMolecule * ) newH;
  limitH++;
  return true;
}

void * Macromolecule::getCurrentH()
{
  if ( heteros == NULL )
    return NULL;
  if ( currentH != -1 )
    return heteros[currentH];
  else
    return NULL;
}

int Macromolecule::endH()
{
  if ( limitH != 0 )
    currentH = limitH - 1;
  else
    currentH = -1;
  return currentH;
}


bool Macromolecule::eraseH( int pos )
{
  if ( ( pos >= 0 ) && ( pos < limitH ) )
  {
    delete heteros[pos];
    return removeH( pos );
  }
  return false;
}


bool Macromolecule::eraseAllH()
{
  int i;
  for ( i = 0; i < limitH; i++ )
    delete heteros[i];
  heteros = ( HeteroMolecule * * ) realloc( heteros, 0 );
  limitH = 0;
  currentH = -1;
  return true;
}

int Macromolecule::getLimitH()
{
  return limitH;
}

Contained * Macromolecule::getH( int pos )
{
  if ( ( pos >= 0 ) && ( pos < limitH ) )
    return heteros[pos];
  else
    return NULL;
}

bool Macromolecule::initH()
{
  if ( limitH != 0 )
  {
    int i;
    currentH = 0;

    for ( i = 0; i < limitH; i++ )
    {
      heteros[i]->initAll();
    }

    return true;
  }
  else
  {
    currentH = -1;
    return false;
  }
}


bool Macromolecule::moveH( Tcoor offset )
{
  int i;
  for ( i = 0; i < limitH; i++ )
    if ( heteros[i]->moveAll( offset ) == false )
      return false;

  return true;
}

bool Macromolecule::nextH()
{
  if ( currentH < ( limitH - 1 ) )
  {
    currentH++;
    return true;
  }
  return false;
}

bool Macromolecule::previousH()
{
  if ( currentH > 0 )
  {
    currentH--;
    return true;
  }
  return false;
}

bool Macromolecule::removeH( int pos )
{
  int i;

  if ( ( pos >= 0 ) && ( pos < limitH ) )
  {
    for ( i = pos; i < ( limitH - 1 ); i++ )
      heteros[i] = heteros[i + 1];


    heteros = ( HeteroMolecule * * ) realloc( heteros, ( sizeof( Contained * ) * limitH ) - ( sizeof( Contained * ) ) );


    if ( currentH > pos )
      currentH--;

    limitH--;

    if ( currentH >= limitH )
      currentH = -1;

    if ( heteros == NULL )
      return false;
    return true;
  }
  return false;
}


bool Macromolecule::removeAllH()
{
  heteros = ( HeteroMolecule * * ) realloc( heteros, 0 );
  limitH = 0;
  currentH = -1;
  return true;
}
