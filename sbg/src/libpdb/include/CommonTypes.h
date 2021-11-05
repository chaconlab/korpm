#ifndef COMMONTYPES_H
#define COMMONTYPES_H

#define PDB_PI 3.1415927
///Declacion de tipos auxiliares
typedef char Tname[20];
typedef float Tcoor[3];


/**
* *@author Jose Ignacio Garzon
*
* Different types of elements
*/
enum Symbol {C,H,N,O,S,P,CA,FE,MG,MN,NA,ZN,NI,CU,K,CO,AL,BR,CL,CR,SI,CD,AU,AG,PT,HG,I,F,D};

/**
* *@author Jose Ignacio Garzon
*
* Object to store the atributes of a region in the Space. This class allows to determine
* the spacital region that is ocupied by a Molecular system
*/
class Dimensions
{
public:

  ///Constructor
  Dimensions()
  {
    Xmax=-1000000;
    Ymax=-1000000;
    Zmax=-1000000;

    Xmin=1000000;
    Ymin=1000000;
    Zmin=1000000;
    geometricCenter[0]=0;
    geometricCenter[1]=0;
    geometricCenter[2]=0;

  };

  ///Copy Constructor
  ///@param d: Reference object to be copied
  Dimensions(Dimensions *d)
  {
    Xmax=d->Xmax;
    Ymax=d->Ymax;
    Zmax=d->Zmax;

    Xmin=d->Xmin;
    Ymin=d->Ymin;
    Zmin=d->Zmin;
    geometricCenter[0]=d->geometricCenter[0];
    geometricCenter[1]=d->geometricCenter[1];
    geometricCenter[2]=d->geometricCenter[2];

  }
  /// X Maximum limit in the space of the region
  float Xmax;
  /// Y Maximum limit in the space of the region
  float Ymax;
  /// Z Maximum limit in the space of the region
  float Zmax;
  /// X Minimum limit in the space of the region
  float Xmin;
  /// Y Minimum limit in the space of the region
  float Ymin;
  /// X Minimum limit in the space of the region
  float Zmin;
  /// Geometric center of the region
  float geometricCenter[3];
};
#endif
