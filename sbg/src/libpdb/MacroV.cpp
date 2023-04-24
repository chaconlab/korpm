

#if defined  _WIN32 || defined _WIN64
#define roundf floor
#endif

#include <libpdb/include/Macromolecule.h>
#include <libpdb/include/pdbIter.h>
// PDB atoms to prevent char warnings
char *Nstr2   = (char *)" N  ";
char *CAstr2  = (char *)" CA ";
char *Cstr2   = (char *)" C  ";
char *Ostr2   = (char *)" O  ";
char *CBstr2  = (char *)" CB ";
char *CGstr2  = (char *)" CG ";
char *CDstr2  = (char *)" CD ";
char *CEstr2  = (char *)" CE ";
char *CZstr2  = (char *)" CZ ";
char *HAstr2  = (char *)" HA ";
char *Hstr2   = (char *)" H  ";

#ifdef VOLUME_INCLUDED

#include <libvolume/include/vlvoliter_linear.h>


namespace FOPS
{
  extern  vlVolume* crop(vlVolume *vol, float threshold, bool self=false);
  extern  vlVolume* createMask(vlVolume *vol, float cutoff);
  extern  vlVolume* projectSurface( vlVolume *vol,float PR, bool inner=false );
  extern   void thresholdUp(vlVolume *vol, float limit);
  extern   void thresholdDown(vlVolume *vol, float limit);
  extern   void threshold(vlVolume *vol, float limit);
  extern   void mul(vlVolume *vol, float value);
  extern  void center_vol(vlVolume *vol, float *center);
  extern  void GaussFilter(vlVolume *vlm,float res);

	extern float calc_total_posneg(vlVolume *vol);
}
int get_element_type(Element *e);

float elect_limited(float point,float limit, float perc, float in300);

float original_ele(float dist, float ron, float roff, float radius, float charge, float charge2);

float autodock_ele(float dist, float ron, float roff, float radius, float charge, int type);




/* Creacion de un objeto volumen que tiene las dimensiones necesarias para almacenar la informacion interpolada del PDB */
vlVolume * Macromolecule::createVolume( float unit )
{
  float widthX, widthY, widthZ;
  int numX, numY, numZ;
  vlVolume * vlm;
	int aux[3];
	Tcoor original_center;
  /* Calculo del centro de geometrico del PDB asi como de las coordenadas maximas y minimas del PDB */
  geoBox();
  geoCenter( original_center );
  // center();

  /*Establecimiento de las dimensiones del PDB que se van a pasar al volumen
  Las dimensiones se redondean por arriba para hacerlas divisibles por el tama�o de voxel
  de esta forma se obtiene un numero exacto de voxels para cubrir todo el espacio*/
  widthX = ceil( ( macromoleculeDimension.Xmax - macromoleculeDimension.Xmin ) /  unit  ) * unit;
  /*Calculo del numero de voxels en un direccion.Se divide la dimension redondeada
  por arriba por el tama�o de voxel y se a�ade un voxel mas (se quiere un voxel al principio y al final de  la dimension)*/
  numX = ( int )( widthX / unit ) + 3;

  widthY = ceil( ( macromoleculeDimension.Ymax - macromoleculeDimension.Ymin ) / unit ) * unit;
  numY = ( int )( widthY / unit ) + 3;

  widthZ = ceil( ( macromoleculeDimension.Zmax - macromoleculeDimension.Zmin ) / unit ) * unit;
  numZ = ( int )( widthZ / unit ) + 3;


	/*aumento para introducir centro*/
	aux[0]=ceil((macromoleculeDimension.geometricCenter[0]-original_center[0]));
	if(aux[0]<0)
		aux[0]*=-1;
	aux[1]=ceil((macromoleculeDimension.geometricCenter[1]-original_center[1]));
	if(aux[1]<0)
		aux[1]*=-1;
	aux[2]=ceil((macromoleculeDimension.geometricCenter[2]-original_center[2]));
	if(aux[2]<0)
		aux[2]*=-1;

	numX+=2*aux[0];
	numY+=2*aux[1];
	numZ+=2*aux[2];


	if(numX%2==0)
		numX++;
	if(numY%2==0)
		numY++;
	if(numZ%2==0)
		numZ++;

  /* Se crea el volumen con las dimensiones calculadas */
  vlm = new vlVolume( vlDim( numX, numY, numZ ), Float, vlUnit( unit, unit, unit ) );
  vlm->clear();
  /* Se devuelve el puntero al volumen */
  return vlm;

}

// create volume plus padding ceros
vlVolume * Macromolecule::createVolume_pad( float unit, int pad )
{
  float widthX, widthY, widthZ;
  int numX, numY, numZ;
  vlVolume * vlm;
	float aux[3];
	Tcoor original_center;
  /* Calculo del centro de geometrico del PDB asi como de las coordenadas maximas y minimas del PDB */
  geoBox();
  geoCenter( original_center );
  // center();

  /*Establecimiento de las dimensiones del PDB que se van a pasar al volumen
  Las dimensiones se redondean por arriba para hacerlas divisibles por el tama�o de voxel
  de esta forma se obtiene un numero exacto de voxels para cubrir todo el espacio*/
  widthX = ceil( ( macromoleculeDimension.Xmax - macromoleculeDimension.Xmin ) / unit ) * unit;
  /*Calculo del numero de voxels en un direccion.Se divide la dimension redondeada
  por arriba por el tama�o de voxel y se a�ade un voxel mas (se quiere un voxel al principio y al final de  la dimension)*/
  numX = ( int )( widthX / unit ) + 1 + pad;

  widthY = ceil( ( macromoleculeDimension.Ymax - macromoleculeDimension.Ymin ) / unit ) * unit;
  numY = ( int )( widthY / unit ) + 1 + pad;

  widthZ = ceil( ( macromoleculeDimension.Zmax - macromoleculeDimension.Zmin ) / unit ) * unit;
  numZ = ( int )( widthZ / unit ) + 1 + pad;

	/*aumento para introducir centro*/
	aux[0]=ceil((macromoleculeDimension.geometricCenter[0]-original_center[0]));
	if(aux[0]<0)
		aux[0]*=-1;
	aux[1]=ceil((macromoleculeDimension.geometricCenter[1]-original_center[1]));
	if(aux[1]<0)
		aux[1]*=-1;
	aux[2]=ceil((macromoleculeDimension.geometricCenter[2]-original_center[2]));
	if(aux[2]<0)
		aux[2]*=-1;

	numX+=2*aux[0];
	numY+=2*aux[1];
	numZ+=2*aux[2];


	if(numX%2==0)
		numX++;
	if(numY%2==0)
		numY++;
	if(numZ%2==0)
		numZ++;

  /* Se crea el volumen con las dimensiones calculadas */
  vlm = new vlVolume( vlDim( numX, numY, numZ ), Float, vlUnit( unit, unit, unit ) );
  vlm->clear();
  /* Se devuelve el puntero al volumen */
  return vlm;

}




/* Creacion de un objeto volumen que tiene las dimensiones necesarias para almacenar la informacion interpolada del PDB
 * Coordenadas de origen establecidas*/
vlVolume * Macromolecule::createVolumeCentered( float unit , bool Hyd)
{
  float widthX, widthY, widthZ;
  int numX, numY, numZ;
  vlVolume * vlm;

	int aux[3];
	Tcoor original_center;
	Tcoor volcenter;
  /* Calculo del centro de geometrico del PDB asi como de las coordenadas maximas y minimas del PDB */
  geoBox();
  if(Hyd)
  {
  	geoCenter( original_center );
  }
  else
  {
  	Macromolecule *pdb2= new Macromolecule(this);
  	pdb2->deleteHYDS();
  	pdb2->geoCenter( original_center );
  	free(pdb2);
  }
  // center();
  //fprintf(stderr,"original center:%f %f %f\n",original_center[0],original_center[1],original_center[2]);


  /*Establecimiento de las dimensiones del PDB que se van a pasar al volumen
  Las dimensiones se redondean por arriba para hacerlas divisibles por el tamanio de voxel
  de esta forma se obtiene un numero exacto de voxels para cubrir todo el espacio*/
  widthX = ceil( ( macromoleculeDimension.Xmax - macromoleculeDimension.Xmin ) / unit ) * unit;
  /*Calculo del numero de voxels en un direccion.Se divide la dimension redondeada
  por arriba por el tama�o de voxel y se a�ade un voxel mas (se quiere un voxel al principio y al final de  la dimension)*/
  numX = ( int )( widthX / unit ) + 3;

  widthY = ceil( ( macromoleculeDimension.Ymax - macromoleculeDimension.Ymin ) / unit ) * unit;
  numY = ( int )( widthY / unit ) + 3;

  widthZ = ceil( ( macromoleculeDimension.Zmax - macromoleculeDimension.Zmin ) / unit ) * unit;
  numZ = ( int )( widthZ / unit ) + 3;

	/*aumento para introdiucir centro*/
	aux[0]=ceil(fabs(macromoleculeDimension.geometricCenter[0]-original_center[0]));
	aux[1]=ceil(fabs(macromoleculeDimension.geometricCenter[1]-original_center[1]));
	aux[2]=ceil(fabs(macromoleculeDimension.geometricCenter[2]-original_center[2]));

	numX+=2*aux[0];
	numY+=2*aux[1];
	numZ+=2*aux[2];


	if(numX%2==0)
		numX++;
	if(numY%2==0)
		numY++;
	if(numZ%2==0)
		numZ++;
  /* Se crea el volumen con las dimensiones calculadas */
  vlm = new vlVolume( vlDim( numX, numY, numZ ), Float, vlUnit( unit, unit, unit ) );
  vlm->clear();
	FOPS::center_vol(vlm,volcenter);
	//fprintf(stderr,"volcenter:%f %f %f\n",volcenter[0],volcenter[1],volcenter[2]);

   /*Calculo del desplazamiento en cada dimension para determinar las nuevas posiciones
  de los atomos en el nuevo origen de coordenadas. Originalmente los atomos estan
  posicionados en funcion del origen de coordenadas real del PDB en amstrong.
  Es necesario posicionar los atomos en funcion del origen de coordenadas del volumen
  en numero de voxels. Para ello hay que tener en cuenta que el PDB esta centrado en
  el volumen en funcion del centro geometrico. Por tanto a cada coordenada de los atomos
  hay que restarle la posicion del centro geometrico (se le posiciona en funcion del centro
  geometrico) y sumarle la distancia que separa el centro geometrico del origen de coordenadas
  del volumen (se posiciona en funcion del origen del volumen. De esta forma se conocera la
  posicion de cada atomo en el volumen). Se calcula el desplazamiento en numero de Voxels */
  //shiftX = ( float )( vlm->dim().x() - 1 ) / 2.0 - macromoleculeDimension.geometricCenter[0] / unit;
  //shiftY = ( float )( vlm->dim().y() - 1 ) / 2.0 - macromoleculeDimension.geometricCenter[1] / unit;
  //shiftZ = ( float )( vlm->dim().z() - 1 ) / 2.0 - macromoleculeDimension.geometricCenter[2] / unit;

  float centerX, centerY, centerZ;

  /*Calculo de la posicion del origen del volumen en el PDB (en funcion del origen real) en
  amstrongs. Para ello se resta el desplamiento con respecto al centro geometrico y luego
  se le suma las coordenadas de dicho centro con respecto al origen real */
  //centerX = ( float ) - unit * ( ( ( vlm->dim().x() ) - 1 ) / 2.0 ) + macromoleculeDimension.geometricCenter[0];
  //centerY = ( float ) - unit * ( ( ( vlm->dim().y() ) - 1 ) / 2.0 ) + macromoleculeDimension.geometricCenter[1];
  //centerZ = ( float ) - unit * ( ( ( vlm->dim().z() ) - 1 ) / 2.0 ) + macromoleculeDimension.geometricCenter[2];
  centerX = ( float ) - unit * volcenter[0] + original_center[0];
  centerY = ( float ) - unit * volcenter[1] + original_center[1];
  centerZ = ( float ) - unit * volcenter[2] + original_center[2];
  vlm->setPosition( vlPoint3f( centerX, centerY, centerZ ) );

  /* Se devuelve el puntero al volumen */
  return vlm;

}


/*Funcion de interpolacion del PDB a un volumen. La funcion devuelve un objeto volumen que contiene una interpolacion del PDB.
El volumen tiene como origen de coordenadas la esquina inferior (en todas las dimensiones)
del volumen, por lo que es necesario realizar una conversion de coordenadas del PDB al volumen
El volumen esta centrado en el volumen con respecto a su CENTRO GEOMETRICO*/
vlVolume * Macromolecule::fillVolumeNE( float unit )
{
  bool end = true;
  Tcoor atomCoor;
  float gx, gy, gz;

  Atom * at;
  vlVolume * vlm;
  Tcoor geometricCenter;
  Tcoor volcenter;
  float weight;
  vlm = createVolume( unit );
  float shiftX, shiftY, shiftZ;
	geoCenter( geometricCenter );
	FOPS::center_vol(vlm,volcenter);
  /*Calculo del desplazamiento en cada dimension para determinar las nuevas posiciones
  de los atomos en el nuevo origen de coordenadas. Originalmente los atomos estan
  posicionados en funcion del origen de coordenadas real del PDB en amstrong.
  Es necesario posicionar los atomos en funcion del origen de coordenadas del volumen
  en numero de voxels. Para ello hay que tener en cuenta que el PDB esta centrado en
  el volumen en funcion del centro geometrico. Por tanto a cada coordenada de los atomos
  hay que restarle la posicion del centro geometrico (se le posiciona en funcion del centro
  geometrico) y sumarle la distancia que separa el centro geometrico del origen de coordenadas
  del volumen (se posiciona en funcion del origen del volumen. De esta forma se conocera la
  posicion de cada atomo en el volumen). Se calcula el desplazamiento en numero de Voxels */
  shiftX = volcenter[0] - geometricCenter[0] / unit;
  shiftY = volcenter[1] - geometricCenter[1] / unit;
  shiftZ = volcenter[2] - geometricCenter[2] / unit;
  float centerX, centerY, centerZ;

  /*Calculo de la posicion del origen del volumen en el PDB (en funcion del origen real) en
  amstrongs. Para ello se resta el desplamiento con respecto al centro geometrico y luego
  se le suma las coordenadas de dicho centro con respecto al origen real */
  centerX = - unit * volcenter[0] + geometricCenter[0];
  centerY = - unit * volcenter[1] + geometricCenter[1];
  centerZ = - unit * volcenter[2] + geometricCenter[2];
  vlm->setPosition( vlPoint3f( centerX, centerY, centerZ ) );


 //  cout << "shift: " << centerX << " " << centerY << " " << centerZ << endl;
 //   cout << "geometric: " << macromoleculeDimension.geometricCenter[0] << " " << macromoleculeDimension.geometricCenter[1]
 //  << " " << macromoleculeDimension.geometricCenter[2] << endl;

  /* Realizacion de la interpolacion trilineal en el volumen para todos los atomos del PDB */
  initAll();

  while ( end != false )
  {
    /* Obtengo un atomo */
    at = getCurrentAtom();
    /* Obtencion de la posicion real del atomo en el PDb en amstrongs */
    at->getPosition(atomCoor);

    /* Calculo de la posicion relativa del atomo en el volumen en numero de voxels */
    gx = ( float )( atomCoor[0] / unit + shiftX );
    gy = ( float )( atomCoor[1] / unit + shiftY );
    gz = ( float )( atomCoor[2] / unit + shiftZ );

    /* Peso atomico del atomo */
    weight = ( at->getElement() )->number;//densidad electronica
		//weight = ( at->getElement() )->weight;//peso molecular

		//if((at->getElement())->symbol==C || (at->getElement())->symbol==O || (at->getElement())->symbol==N )
    /* Para cada vecino calculo trilineal de la cantidad de peso atomico del atomo que le afecta */
	    vlm->interpolateVoxel(vlPoint3f(gx,gy,gz), weight);
			//vlm->setVoxel(vlPoint3ui(gx,gy,gz), weight);
    end = nextAtom();
  }

  return vlm;
}

vlVolume * Macromolecule::fillVolumeBFS( float unit, float *bfs )
{
  bool end = true;
  Tcoor atomCoor;
  float gx, gy, gz;

  Atom * at;
  vlVolume * vlm;
  Tcoor geometricCenter;
  Tcoor volcenter;
  vlm = createVolume( unit );
  float shiftX, shiftY, shiftZ;
	geoCenter( geometricCenter );
	FOPS::center_vol(vlm,volcenter);
  /*Calculo del desplazamiento en cada dimension para determinar las nuevas posiciones
  de los atomos en el nuevo origen de coordenadas. Originalmente los atomos estan
  posicionados en funcion del origen de coordenadas real del PDB en amstrong.
  Es necesario posicionar los atomos en funcion del origen de coordenadas del volumen
  en numero de voxels. Para ello hay que tener en cuenta que el PDB esta centrado en
  el volumen en funcion del centro geometrico. Por tanto a cada coordenada de los atomos
  hay que restarle la posicion del centro geometrico (se le posiciona en funcion del centro
  geometrico) y sumarle la distancia que separa el centro geometrico del origen de coordenadas
  del volumen (se posiciona en funcion del origen del volumen. De esta forma se conocera la
  posicion de cada atomo en el volumen). Se calcula el desplazamiento en numero de Voxels */
  shiftX = volcenter[0] - geometricCenter[0] / unit;
  shiftY = volcenter[1] - geometricCenter[1] / unit;
  shiftZ = volcenter[2] - geometricCenter[2] / unit;
  float centerX, centerY, centerZ;

  /*Calculo de la posicion del origen del volumen en el PDB (en funcion del origen real) en
  amstrongs. Para ello se resta el desplamiento con respecto al centro geometrico y luego
  se le suma las coordenadas de dicho centro con respecto al origen real */
  centerX = - unit * volcenter[0] + geometricCenter[0];
  centerY = - unit * volcenter[1] + geometricCenter[1];
  centerZ = - unit * volcenter[2] + geometricCenter[2];
  vlm->setPosition( vlPoint3f( centerX, centerY, centerZ ) );


 //  cout << "shift: " << centerX << " " << centerY << " " << centerZ << endl;
 //   cout << "geometric: " << macromoleculeDimension.geometricCenter[0] << " " << macromoleculeDimension.geometricCenter[1]
 //  << " " << macromoleculeDimension.geometricCenter[2] << endl;

  /* Realizacion de la interpolacion trilineal en el volumen para todos los atomos del PDB */
  initAll();
  int cont=0;
  while ( end != false )
  {
    /* Obtengo un atomo */
    at = getCurrentAtom();
    /* Obtencion de la posicion real del atomo en el PDb en amstrongs */
    at->getPosition(atomCoor);

    /* Calculo de la posicion relativa del atomo en el volumen en numero de voxels */
    gx = ( float )( atomCoor[0] / unit + shiftX );
    gy = ( float )( atomCoor[1] / unit + shiftY );
    gz = ( float )( atomCoor[2] / unit + shiftZ );

    /* Para cada vecino calculo trilineal de la cantidad de peso atomico del atomo que le afecta */
	vlm->interpolateVoxel(vlPoint3f(gx,gy,gz), bfs[cont]);
	cont++;
			//vlm->setVoxel(vlPoint3ui(gx,gy,gz), weight);
    end = nextAtom();
  }

  return vlm;
}



/*Funcion de interpolacion del PDB a un volumen. La funcion devuelve un objeto volumen que contiene una interpolacion del PDB.
El volumen tiene como origen de coordenadas la esquina inferior (en todas las dimensiones)
del volumen, por lo que es necesario realizar una conversion de coordenadas del PDB al volumen
El volumen esta centrado en el volumen con respecto a su centro geometrico*/
vlVolume * Macromolecule::fillVolumeNE_3BB2R( float unit )
{
  bool end = true;
  Tcoor atomCoor;
  float gx, gy, gz;

  Atom * at;
  vlVolume * vlm;
  float weight;
  vlm = createVolume( unit );
  Tcoor volcenter;
  FOPS::center_vol(vlm,volcenter);
  float shiftX, shiftY, shiftZ;
  float centerX, centerY, centerZ;

  /*Calculo del desplazamiento en cada dimension para determinar las nuevas posiciones
  de los atomos en el nuevo origen de coordenadas. Originalmente los atomos estan
  posicionados en funcion del origen de coordenadas real del PDB en amstrong.
  Es necesario posicionar los atomos en funcion del origen de coordenadas del volumen
  en numero de voxels. Para ello hay que tener en cuenta que el PDB esta centrado en
  el volumen en funcion del centro geometrico. Por tanto a cada coordenada de los atomos
  hay que restarle la posicion del centro geometrico (se le posiciona en funcion del centro
  geometrico) y sumarle la distancia que separa el centro geometrico del origen de coordenadas
  del volumen (se posiciona en funcion del origen del volumen. De esta forma se conocera la
  posicion de cada atomo en el volumen). Se calcula el desplazamiento en numero de Voxels */
  shiftX = volcenter[0] - macromoleculeDimension.geometricCenter[0] / unit;
  shiftY = volcenter[1] - macromoleculeDimension.geometricCenter[1] / unit;
  shiftZ = volcenter[2] - macromoleculeDimension.geometricCenter[2] / unit;

  /*Calculo de la posicion del origen del volumen en el PDB (en funcion del origen real) en
  amstrongs. Para ello se resta el desplamiento con respecto al centro geometrico y luego
  se le suma las coordenadas de dicho centro con respecto al origen real */
  centerX = ( float ) - unit * volcenter[0] + macromoleculeDimension.geometricCenter[0];
  centerY = ( float ) - unit * volcenter[1] + macromoleculeDimension.geometricCenter[1];
  centerZ = ( float ) - unit * volcenter[2] + macromoleculeDimension.geometricCenter[2];
  vlm->setPosition( vlPoint3f( centerX, centerY, centerZ ) );

  //  cout << "shift: " << centerX << " " << centerY << " " << centerZ << endl;
  //  cout << "geometric: " << macromoleculeDimension.geometricCenter[0] << " " << macromoleculeDimension.geometricCenter[1]
  // << " " << macromoleculeDimension.geometricCenter[2] << endl;

  /* Realizacion de la interpolacion trilineal en el volumen para todos los atomos del PDB */
  initAll();
  while ( end != false )
  {
    /* Obtengo un atomo */
    at = getCurrentAtom();
    /* Obtencion de la posicion real del atomo en el PDb en amstrongs */
    at->getPosition(atomCoor);

    /* Calculo de la posicion relativa del atomo en el volumen en numero de voxels */
    gx = ( float )( atomCoor[0] / unit + shiftX );
    gy = ( float )( atomCoor[1] / unit + shiftY );
    gz = ( float )( atomCoor[2] / unit + shiftZ );

    /* Peso atomico del atomo */
	weight = at->getPdbfact(); // Number of Electrons from current atomic model
//	weight = at->getPdbocc(); // Mass from current atomic model
//	printf("Reading atomic charge: %s --> %f\n",at->getName(),weight);

//	weight = at->getPdbocc(); // Mass from the 3BB2R reduced model
//		printf("Reading atom weight: %s --> %f\n",at->getName(),weight);

	// weight = ( at->getElement() )->number;//densidad electronica
	// weight = ( at->getElement() )->weight;//peso molecular

		//if((at->getElement())->symbol==C || (at->getElement())->symbol==O || (at->getElement())->symbol==N )
    /* Para cada vecino calculo trilineal de la cantidad de peso atomico del atomo que le afecta */
	    vlm->interpolateVoxel(vlPoint3f(gx,gy,gz), weight);
			//vlm->setVoxel(vlPoint3ui(gx,gy,gz), weight);
    end = nextAtom();
  }

  return vlm;
}

/*Funcion de interpolacion del PDB a un volumen. La funcion devuelve un objeto volumen que contiene una interpolacion del PDB.
El volumen tiene como origen de coordenadas la esquina inferior (en todas las dimensiones)
del volumen, por lo que es necesario realizar una conversion de coordenadas del PDB al volumen
El volumen esta centrado en el volumen con respecto a su centro geometrico*/
vlVolume * Macromolecule::fillVolumeMW( float unit )
{
  bool end = true;
  Tcoor atomCoor;
  float gx, gy, gz;

  Atom * at;
  vlVolume * vlm;
  float weight;
  vlm = createVolume( unit );
  Tcoor volcenter;
  FOPS::center_vol(vlm,volcenter);
  float shiftX, shiftY, shiftZ;

  /*Calculo del desplazamiento en cada dimension para determinar las nuevas posiciones
  de los atomos en el nuevo origen de coordenadas. Originalmente los atomos estan
  posicionados en funcion del origen de coordenadas real del PDB en amstrong.
  Es necesario posicionar los atomos en funcion del origen de coordenadas del volumen
  en numero de voxels. Para ello hay que tener en cuenta que el PDB esta centrado en
  el volumen en funcion del centro geometrico. Por tanto a cada coordenada de los atomos
  hay que restarle la posicion del centro geometrico (se le posiciona en funcion del centro
  geometrico) y sumarle la distancia que separa el centro geometrico del origen de coordenadas
  del volumen (se posiciona en funcion del origen del volumen. De esta forma se conocera la
  posicion de cada atomo en el volumen). Se calcula el desplazamiento en numero de Voxels */
  shiftX = volcenter[0] - macromoleculeDimension.geometricCenter[0] / unit;
  shiftY = volcenter[1] - macromoleculeDimension.geometricCenter[1] / unit;
  shiftZ = volcenter[2] - macromoleculeDimension.geometricCenter[2] / unit;
  float centerX, centerY, centerZ;

  /*Calculo de la posicion del origen del volumen en el PDB (en funcion del origen real) en
  amstrongs. Para ello se resta el desplamiento con respecto al centro geometrico y luego
  se le suma las coordenadas de dicho centro con respecto al origen real */
  centerX = ( float ) - unit * volcenter[0] + macromoleculeDimension.geometricCenter[0];
  centerY = ( float ) - unit * volcenter[1] + macromoleculeDimension.geometricCenter[1];
  centerZ = ( float ) - unit * volcenter[2] + macromoleculeDimension.geometricCenter[2];
  vlm->setPosition( vlPoint3f( centerX, centerY, centerZ ) );


  //  cout << "shift: " << centerX << " " << centerY << " " << centerZ << endl;
  //  cout << "geometric: " << macromoleculeDimension.geometricCenter[0] << " " << macromoleculeDimension.geometricCenter[1]
  // << " " << macromoleculeDimension.geometricCenter[2] << endl;

  /* Realizacion de la interpolacion trilineal en el volumen para todos los atomos del PDB */
  initAll();
  while ( end != false )
  {
    /* Obtengo un atomo */
    at = getCurrentAtom();
    /* Obtencion de la posicion real del atomo en el PDb en amstrongs */
    at->getPosition(atomCoor);

    /* Calculo de la posicion relativa del atomo en el volumen en numero de voxels */
    gx = ( float )( atomCoor[0] / unit + shiftX );
    gy = ( float )( atomCoor[1] / unit + shiftY );
    gz = ( float )( atomCoor[2] / unit + shiftZ );

    /* Peso atomico del atomo */
    //weight = ( at->getElement() )->number;//densidad electronica
		weight = ( at->getElement() )->weight;//peso molecular

		//if((at->getElement())->symbol==C || (at->getElement())->symbol==O || (at->getElement())->symbol==N )
    /* Para cada vecino calculo trilineal de la cantidad de peso atomico del atomo que le afecta */
	    vlm->interpolateVoxel(vlPoint3f(gx,gy,gz), weight);
			//vlm->setVoxel(vlPoint3ui(gx,gy,gz), weight);
    end = nextAtom();
  }

  return vlm;
}


/*Funcion de interpolacion del PDB a un volumen. La funcion devuelve un objeto volumen que contiene una interpolacion del PDB.
El volumen tiene como origen de coordenadas la esquina inferior (en todas las dimensiones)
del volumen, por lo que es necesario realizar una conversion de coordenadas del PDB al volumen
El volumen esta centrado en el volumen con respecto a su centro geometrico*/
vlVolume * Macromolecule::fillVolume( float maxLenght, float unit )
{
  bool end = true;
  Tcoor atomCoor;
  float gx, gy, gz;
  int aux;
  float shiftX, shiftY, shiftZ;
  Atom * at;
  vlVolume * vlm;
  float weight;


  aux=ceil(maxLenght/unit);
  vlm=new vlVolume( vlDim( aux, aux, aux ), Float, vlUnit( unit, unit, unit ) );
  Tcoor volcenter;
  FOPS::center_vol(vlm,volcenter);


  /*Calculo del desplazamiento en cada dimension para determinar las nuevas posiciones
  de los atomos en el nuevo origen de coordenadas. Originalmente los atomos estan
  posicionados en funcion del origen de coordenadas real del PDB en amstrong.
  Es necesario posicionar los atomos en funcion del origen de coordenadas del volumen
  en numero de voxels. Para ello hay que tener en cuenta que el PDB esta centrado en
  el volumen en funcion del centro geometrico. Por tanto a cada coordenada de los atomos
  hay que restarle la posicion del centro geometrico (se le posiciona en funcion del centro
  geometrico) y sumarle la distancia que separa el centro geometrico del origen de coordenadas
  del volumen (se posiciona en funcion del origen del volumen. De esta forma se conocera la
  posicion de cada atomo en el volumen). Se calcula el desplazamiento en numero de Voxels */
  shiftX = volcenter[0] - macromoleculeDimension.geometricCenter[0] / unit;
  shiftY = volcenter[1] - macromoleculeDimension.geometricCenter[1] / unit;
  shiftZ = volcenter[2] - macromoleculeDimension.geometricCenter[2] / unit;
  float centerX, centerY, centerZ;

  /*Calculo de la posicion del origen del volumen en el PDB (en funcion del origen real) en
  amstrongs. Para ello se resta el desplamiento con respecto al centro geometrico y luego
  se le suma las coordenadas de dicho centro con respecto al origen real */
  centerX = ( float ) - unit * volcenter[0] + macromoleculeDimension.geometricCenter[0];
  centerY = ( float ) - unit * volcenter[1] + macromoleculeDimension.geometricCenter[1];
  centerZ = ( float ) - unit * volcenter[2] + macromoleculeDimension.geometricCenter[2];
  vlm->setPosition( vlPoint3f( centerX, centerY, centerZ ) );


  //  cout << "shift: " << centerX << " " << centerY << " " << centerZ << endl;
  //  cout << "geometric: " << macromoleculeDimension.geometricCenter[0] << " " << macromoleculeDimension.geometricCenter[1]
  // << " " << macromoleculeDimension.geometricCenter[2] << endl;

  /* Realizacion de la interpolacion trilineal en el volumen para todos los atomos del PDB */
  initAll();
  while ( end != false )
  {
    /* Obtengo un atomo */
    at = getCurrentAtom();
    /* Obtencion de la posicion real del atomo en el PDb en amstrongs */
    at->getPosition(atomCoor);

    /* Calculo de la posicion relativa del atomo en el volumen en numero de voxels */
    gx = atomCoor[0] / unit + shiftX;
    gy = atomCoor[1] / unit + shiftY;
    gz = atomCoor[2] / unit + shiftZ;

    /* Peso atomico del atomo */
    weight = ( float )( at->getElement() )->number;

    vlm->interpolateVoxel(vlPoint3f(gx,gy,gz), weight);

    end = nextAtom();
  }
  return vlm;
}

/*Funcion de interpolacion de valores Beta del PDB a un volumen. La funcion devuelve un objeto volumen que contiene una interpolacion de los valores
 de Beta (fact) del PDB.
El volumen tiene como origen de coordenadas la esquina inferior (en todas las dimensiones)
del volumen, por lo que es necesario realizar una conversion de coordenadas del PDB al volumen
El volumen esta centrado en el volumen con respecto a su centro geometrico*/
vlVolume * Macromolecule::fillVolumeBeta( float unit )
{
  bool end = true;
  Tcoor atomCoor;
  float gx, gy, gz;

  Atom * at;
  vlVolume * vlm;
  float weight;
  vlm = createVolume( unit );
  float shiftX, shiftY, shiftZ;
  Tcoor volcenter;
  FOPS::center_vol(vlm,volcenter);

  /*Calculo del desplazamiento en cada dimension para determinar las nuevas posiciones
  de los atomos en el nuevo origen de coordenadas. Originalmente los atomos estan
  posicionados en funcion del origen de coordenadas real del PDB en amstrong.
  Es necesario posicionar los atomos en funcion del origen de coordenadas del volumen
  en numero de voxels. Para ello hay que tener en cuenta que el PDB esta centrado en
  el volumen en funcion del centro geometrico. Por tanto a cada coordenada de los atomos
  hay que restarle la posicion del centro geometrico (se le posiciona en funcion del centro
  geometrico) y sumarle la distancia que separa el centro geometrico del origen de coordenadas
  del volumen (se posiciona en funcion del origen del volumen. De esta forma se conocera la
  posicion de cada atomo en el volumen). Se calcula el desplazamiento en numero de Voxels */
  shiftX = volcenter[0] - macromoleculeDimension.geometricCenter[0] / unit;
  shiftY = volcenter[1] - macromoleculeDimension.geometricCenter[1] / unit;
  shiftZ = volcenter[2] - macromoleculeDimension.geometricCenter[2] / unit;
  float centerX, centerY, centerZ;

  /*Calculo de la posicion del origen del volumen en el PDB (en funcion del origen real) en
  amstrongs. Para ello se resta el desplamiento con respecto al centro geometrico y luego
  se le suma las coordenadas de dicho centro con respecto al origen real */
  centerX = ( float ) - unit * volcenter[0] + macromoleculeDimension.geometricCenter[0];
  centerY = ( float ) - unit * volcenter[1] + macromoleculeDimension.geometricCenter[1];
  centerZ = ( float ) - unit * volcenter[2] + macromoleculeDimension.geometricCenter[2];
  vlm->setPosition( vlPoint3f( centerX, centerY, centerZ ) );


  //  cout << "shift: " << centerX << " " << centerY << " " << centerZ << endl;
  //  cout << "geometric: " << macromoleculeDimension.geometricCenter[0] << " " << macromoleculeDimension.geometricCenter[1]
  // << " " << macromoleculeDimension.geometricCenter[2] << endl;

  /* Realizacion de la interpolacion trilineal en el volumen para todos los atomos del PDB */
  initAll();
  while ( end != false )
  {
    /* Obtengo un atomo */
    at = getCurrentAtom();
    /* Obtencion de la posicion real del atomo en el PDb en amstrongs */
    at->getPosition(atomCoor);

    /* Calculo de la posicion relativa del atomo en el volumen en numero de voxels */
    gx = ( float )( atomCoor[0] / unit + shiftX );
    gy = ( float )( atomCoor[1] / unit + shiftY );
    gz = ( float )( atomCoor[2] / unit + shiftZ );

    /* Peso atomico del atomo */
   // weight = ( at->getElement() )->number;
     weight = at->getPdbfact();

    /* Para cada vecino calculo trilineal de la cantidad de peso atomico del atomo que le afecta */
    vlm->interpolateVoxel(vlPoint3f(gx,gy,gz), weight);
    end = nextAtom();
  }

  return vlm;
}



void Macromolecule::project( vlVolume * vlm )
{
  bool end = true;
  Tcoor atomCoor;
  float gx, gy, gz;
  Atom * at;

  float weight;

  float shiftX, shiftY, shiftZ;
  Tcoor volcenter;
  FOPS::center_vol(vlm,volcenter);

  /*Calculo del desplazamiento en cada dimension para determinar las nuevas posiciones
  de los atomos en el nuevo origen de coordenadas. Originalmente los atomos estan
  posicionados en funcion del origen de coordenadas real del PDB en amstrong.
  Es necesario posicionar los atomos en funcion del origen de coordenadas del volumen
  en numero de voxels. Para ello hay que tener en cuenta que el PDB esta centrado en
  el volumen en funcion del centro geometrico. Por tanto a cada coordenada de los atomos
  hay que restarle la posicion del centro geometrico (se le posiciona en funcion del centro
  geometrico) y sumarle la distancia que separa el centro geometrico del origen de coordenadas
  del volumen (se posiciona en funcion del origen del volumen. De esta forma se conocera la
  posicion de cada atomo en el volumen). Se calcula el desplazamiento en numero de Voxels */
  shiftX = volcenter[0] - macromoleculeDimension.geometricCenter[0] / vlm->units().x();
  shiftY = volcenter[1] - macromoleculeDimension.geometricCenter[1] / vlm->units().y();
  shiftZ = volcenter[2] - macromoleculeDimension.geometricCenter[2] / vlm->units().z();
  float centerX, centerY, centerZ;

  /*Calculo de la posicion del origen del volumen en el PDB (en funcion del origen real) en
  amstrongs. Para ello se resta el desplamiento con respecto al centro geometrico y luego
  se le suma las coordenadas de dicho centro con respecto al origen real */
  centerX = ( float ) - vlm->units().x() * volcenter[0] + macromoleculeDimension.geometricCenter[0];
  centerY = ( float ) - vlm->units().y() * volcenter[1] + macromoleculeDimension.geometricCenter[1];
  centerZ = ( float ) - vlm->units().z() * volcenter[2] + macromoleculeDimension.geometricCenter[2];
  vlm->setPosition( vlPoint3f( centerX, centerY, centerZ ) );


  //  cout << "shift: " << centerX << " " << centerY << " " << centerZ << endl;
  //  cout << "geometric: " << macromoleculeDimension.geometricCenter[0] << " " << macromoleculeDimension.geometricCenter[1]
  // << " " << macromoleculeDimension.geometricCenter[2] << endl;

  /* Realizacion de la interpolacion trilineal en el volumen para todos los atomos del PDB */
  initAll();
  while ( end != false )
  {
    /* Obtengo un atomo */
    at = getCurrentAtom();
    /* Obtencion de la posicion real del atomo en el PDb en amstrongs */
    at->getPosition(atomCoor);

    /* Calculo de la posicion relativa del atomo en el volumen en numero de voxels */
    gx = atomCoor[0] / vlm->units().x() + shiftX;
    gy = atomCoor[1] / vlm->units().y() + shiftY;
    gz = atomCoor[2] / vlm->units().z() + shiftZ;

    /* Peso atomico del atomo */
    weight = ( float )( at->getElement() )->number;
    vlm->interpolateVoxel(vlPoint3f(gx,gy,gz), weight);
    end = nextAtom();
  }
}

void Macromolecule::project_Real( vlVolume * vlm )
{
  bool end = true;
  Tcoor atomCoor;
  float gx, gy, gz;
  Atom * at;

  float weight;

  vlPoint3f pos;
  vlm->getPosition(&pos);
  /* Realizacion de la interpolacion trilineal en el volumen para todos los atomos del PDB */
  initAll();
  while ( end != false )
  {
    /* Obtengo un atomo */
    at = getCurrentAtom();
    /* Obtencion de la posicion real del atomo en el PDb en amstrongs */
    at->getPosition(atomCoor);

    /* Calculo de la posicion relativa del atomo en el volumen en numero de voxels */
    gx = ((atomCoor[0]- pos.x()) / vlm->units().x()) ;
    gy = ((atomCoor[1]- pos.y()) / vlm->units().y()) ;
    gz = ((atomCoor[2]- pos.z()) / vlm->units().z()) ;

    /* Peso atomico del atomo */
    weight = ( float )( at->getElement() )->number;
    vlm->interpolateVoxel(vlPoint3f(gx,gy,gz), weight);
    end = nextAtom();
  }
}

void Macromolecule::project_unit( vlVolume * vlm )
{
  bool end = true;
  Tcoor atomCoor;
  float gx, gy, gz;
  Atom * at;

  float weight=1.0;

  vlPoint3f pos;
  vlm->getPosition(&pos);
  /* Realizacion de la interpolacion trilineal en el volumen para todos los atomos del PDB */
  initAll();
  while ( end != false )
  {
    /* Obtengo un atomo */
    at = getCurrentAtom();
    /* Obtencion de la posicion real del atomo en el PDb en amstrongs */
    at->getPosition(atomCoor);

    /* Calculo de la posicion relativa del atomo en el volumen en numero de voxels */
    gx = ((atomCoor[0]- pos.x()) / vlm->units().x()) ;
    gy = ((atomCoor[1]- pos.y()) / vlm->units().y()) ;
    gz = ((atomCoor[2]- pos.z()) / vlm->units().z()) ;

    vlm->interpolateVoxel(vlPoint3f(gx,gy,gz), weight);
    end = nextAtom();
  }
}

//// Mon made (29/5/2008)
//// Projects over an pre-existing vlVolume a Macromolecule with NE and without moving
//// the reference system.
//void Macromolecule::project_RealNE_3BB2R( vlVolume * vlm)
//{
//  bool end = true;
//  Tcoor atomCoor;
//  float gx, gy, gz;
//  Atom * at;
//
//  float weight;
//
//  vlPoint3f pos;
//  vlm->getPosition(&pos);
//  /* Realizacion de la interpolacion trilineal en el volumen para todos los atomos del PDB */
//  initAll();
//  while ( end != false )
//  {
//    /* Obtengo un atomo */
//    at = getCurrentAtom();
//    /* Obtencion de la posicion real del atomo en el PDb en amstrongs */
//    at->getPosition(atomCoor);
//
//    /* Calculo de la posicion relativa del atomo en el volumen en numero de voxels */
//    gx = ((atomCoor[0]- pos.x()) / vlm->units().x()) ;
//    gy = ((atomCoor[1]- pos.y()) / vlm->units().y()) ;
//    gz = ((atomCoor[2]- pos.z()) / vlm->units().z()) ;
//
//    /* Peso atomico del atomo */
////    weight = ( float )( at->getElement() )->number;
//	weight = at->getPdbfact(); // Number of Electrons from the 3BB2R reduced model
//
//    vlm->interpolateVoxel(vlPoint3f(gx,gy,gz), weight);
//    end = nextAtom();
//  }
//}


// Mon modified (2/11/2010)
// Projects over an pre-existing vlVolume a Macromolecule with NE and without moving
// the reference system.
// If fast=true, a fast method will be used instead of trilinear-interpolation (default=false)
void Macromolecule::project_RealNE_3BB2R( vlVolume * vlm, bool fast )
{
	bool end = true;
	Tcoor atomCoor;
	float gx, gy, gz;
	Atom *at;
	vlPoint3f pos;
	vlm->getPosition(&pos);
	vlStep m_step = vlm->stepping();
	vlDim m_dimLimit = vlm->dim();
	float *m_pData = (float*)vlm->getVoxelVoidPtr(vlPoint3ui(0,0,0));
	float voxel;
	int stepy = m_step.y();
	int stepz = m_step.z();
	int limitx = m_dimLimit.x();
	int limity = m_dimLimit.y();
	int limitz = m_dimLimit.z();
	int offset;
	int x0,x1,y0,y1,z0,z1;
	float a,b,c,ab,ab1,a1b,a1b1;
	float unitsx = vlm->units().x();
	float unitsy = vlm->units().y();
	float unitsz = vlm->units().z();

	initAll();

	if(fast) // Fast method, i.e. no-trilinear interpolation
	{
		while ( end != false )
		{
			at = getCurrentAtom(); // Get atom
			at->getPosition(atomCoor); // Get atom cartesian coordinates

			// Computing the relative position of current atom inside the volume (in voxel units)
			gx = ((atomCoor[0]- pos.x()) / unitsx);
			gy = ((atomCoor[1]- pos.y()) / unitsy);
			gz = ((atomCoor[2]- pos.z()) / unitsz);

			// Por si acaso...
			if( !(gx > limitx || gy > limity || gz > limitz
					|| gx < 0 || gy < 0 || gz < 0) )
			{
				/* Peso atomico del atomo */
				//    weight = ( float )( at->getElement() )->number;
				voxel = at->getPdbfact(); // Number of Electrons from current atomic model
//				voxel = at->getPdbocc(); // Mass from current atomic model
//			    voxel = ( at->getElement() )->number; // electron density (taken from atomic elements table)

				x0=(int)roundf(gx);
				y0=(int)roundf(gy);
				z0=(int)roundf(gz);
				*(m_pData + (z0*stepz+y0*stepy+x0)) += voxel;
			}
			end = nextAtom();
		}
	}
	else
	{
		// Trilinear interpolation for every atom in PDB (not so fast)
		while ( end != false )
		{
			at = getCurrentAtom(); // Get atom
			at->getPosition(atomCoor); // Get atom cartesian coordinates
			// Computing the relative position of current atom inside the volume (in voxel units)
			gx = ((atomCoor[0]- pos.x()) / unitsx);
			gy = ((atomCoor[1]- pos.y()) / unitsy);
			gz = ((atomCoor[2]- pos.z()) / unitsz);

			// Por si acaso...
			if( !(gx > limitx || gy > limity || gz > limitz
					|| gx < 0 || gy < 0 || gz < 0) )
			{
				/* Peso atomico del atomo */
				//    weight = ( float )( at->getElement() )->number;
				voxel = at->getPdbfact(); // Number of Electrons from current atomic model
//				voxel = at->getPdbocc(); // Mass from current atomic model
//			    voxel = ( at->getElement() )->number; // electron density (taken from atomic elements table)
				x0=(int)floorf(gx);
				y0=(int)floorf(gy);
				z0=(int)floorf(gz);
				x1=x0+1;
				y1=y0+1;
				z1=z0+1;
				a = x1 - gx; // = x0 + 1 - gx  <-- distance from gx to i+1 voxel center
				b = y1 - gy; // = y0 + 1 - gy  <-- distance from gy to j+1 voxel center
				c = z1 - gz; // = z0 + 1 - gz  <-- distance from gz to k+1 voxel center
				// (1 - a) = 1-x0-1+gx = gx-x0 --> distance from current voxel center to gx
				// (1 - b) = 1-y0-1+gy = gy-y0 --> distance from current voxel center to gy
				ab = a * b;
				ab1 = a * ( 1 - b );
				a1b = ( 1 - a ) * b;
				a1b1 = ( 1 - a ) * ( 1 - b );

				// now, dummy variables
				a = ( 1 - c ) * voxel;
				c = (z1 - gz) * voxel;
				offset = z0*stepz + y0*stepy + x0; // initial offset

				// set the voxel data
//				*(m_pData + (z0*stepz+y0*stepy+x0)) += ab*c*voxel;
				*(m_pData + offset) += ab*c;
				if(x1<limitx)
//					*(m_pData + (z0*stepz+y0*stepy+x1)) += a1b*c*voxel;
					*(m_pData + (offset + 1)) += a1b*c;

				if(y1<limity)
				{
//					*(m_pData + (z0*stepz+y1*stepy+x0)) += ab1*c*voxel;
					*(m_pData + (offset + stepy)) += ab1*c;
					if(x1<limitx)
//						*(m_pData + (z0*stepz+y1*stepy+x1)) += a1b1*c*voxel;
						*(m_pData + (offset + stepy + 1)) += a1b1*c;
				}

				if(z1<limitz)
				{
//					*(m_pData + (z1*stepz+y0*stepy+x0)) += ab*a*voxel;
					*(m_pData + (offset + stepz)) += ab*a;
					if(x1<limitx)
//						*(m_pData + (z1*stepz+y0*stepy+x1)) += a1b*a*voxel;
						*(m_pData + (offset + stepz + 1)) += a1b*a;
					if(y1<limity)
					{
//						*(m_pData + (z1*stepz+y1*stepy+x0)) += ab1*a*voxel;
						*(m_pData + (offset + stepz + stepy)) += ab1*a;
						if(x1<limitx)
//							*(m_pData + (z1*stepz+y1*stepy+x1)) += a1b1*a*voxel;
							*(m_pData + (offset + stepz + stepy + 1)) += a1b1*a;
					}
				}
			}
			end = nextAtom();
		}
	}
}

#ifdef USE_PTHREAD // Enables PThread parallel routines
// Mon made (29/05/2013)
// Projects over an pre-existing vlVolume a Macromolecule with NE and without moving
// the reference system in parallel. (Thread routine)
// If fast=true, a fast method will be used instead of trilinear-interpolation (default=false)
void *project_RealNE_3BB2R_thread( void *threadarg )
{
	Tcoor atomCoor;
	float gx, gy, gz;
	Atom *at;
	vlPoint3f pos;
	float voxel;
	int offset;
	int x0,x1,y0,y1,z0,z1;
	float a,b,c,ab,ab1,a1b,a1b1;

//	fprintf(stderr,"Hi, I've just born! Thread %u\n", (unsigned int) pthread_self() );

	//	initAll();

	while(true) // loops for ever...
	{
//		fprintf(stderr,"Hi, I'm waiting for work... Thread %u\n", (unsigned int) pthread_self() );
		pthread_mutex_lock(((project_RealNE_data *) threadarg)->p_mutex_begin);
		while( ((project_RealNE_data *) threadarg)->begin ) // Checking whether the thread would begin ("re-signal" for "wait" safety)
			pthread_cond_wait(((project_RealNE_data *) threadarg)->p_cond_begin, ((project_RealNE_data *) threadarg)->p_mutex_begin);
		((project_RealNE_data *) threadarg)->begin = true; // for next iteration... (enables "beginability")
		pthread_mutex_unlock(((project_RealNE_data *) threadarg)->p_mutex_begin);
//		fprintf(stderr,"Hi, I'm thread %u awaking\n", (unsigned int) pthread_self() );

		// Loading data from "thread argument"
		vlVolume *vlm = ((project_RealNE_data *) threadarg)->vlm;
		int stepy = ((project_RealNE_data *) threadarg)->stepy;
		int stepz = ((project_RealNE_data *) threadarg)->stepz;
		int limitx = ((project_RealNE_data *) threadarg)->limitx;
		int limity = ((project_RealNE_data *) threadarg)->limity;
		int limitz = ((project_RealNE_data *) threadarg)->limitz;
		float *m_pData = (float*)vlm->getVoxelVoidPtr(vlPoint3ui(0,0,0));
		float unitsx = ((project_RealNE_data *) threadarg)->unitsx;
		float unitsy = ((project_RealNE_data *) threadarg)->unitsy;
		float unitsz = ((project_RealNE_data *) threadarg)->unitsz;
		int first = ((project_RealNE_data *) threadarg)->first;
		int last = ((project_RealNE_data *) threadarg)->last;
		// WARNING: Macromolecule casting is mandatory because some "dichotomy" occurs in "project_RealNE_data" declaration...)
		pdbIter *iter = (pdbIter *)(((project_RealNE_data *) threadarg)->iter); // load Macromolecule iterator
		vlm->getPosition(&pos);

		// Projecting pdb into map by Trilinear interpolation
		for(iter->pos_atom = first; iter->pos_atom < last; iter->pos_atom++ )
		{
			//	while ( end != false )
			//	{
			at = iter->get_atom(); // Get atom
			at->getPosition(atomCoor); // Get atom cartesian coordinates
			// Computing the relative position of current atom inside the volume (in voxel units)
			gx = ((atomCoor[0]- pos.x()) / unitsx);
			gy = ((atomCoor[1]- pos.y()) / unitsy);
			gz = ((atomCoor[2]- pos.z()) / unitsz);

			// Por si acaso...
			if( !(gx > limitx || gy > limity || gz > limitz
					|| gx < 0 || gy < 0 || gz < 0) )
			{
				/* Peso atomico del atomo */
				//    weight = ( float )( at->getElement() )->number;
				voxel = at->getPdbfact(); // Number of Electrons from current atomic model
				//				voxel = at->getPdbocc(); // Mass from current atomic model
				//			    voxel = ( at->getElement() )->number; // electron density (taken from atomic elements table)
				x0=(int)floorf(gx);
				y0=(int)floorf(gy);
				z0=(int)floorf(gz);
				x1=x0+1;
				y1=y0+1;
				z1=z0+1;
				a = x1 - gx; // = x0 + 1 - gx  <-- distance from gx to i+1 voxel center
				b = y1 - gy; // = y0 + 1 - gy  <-- distance from gy to j+1 voxel center
				c = z1 - gz; // = z0 + 1 - gz  <-- distance from gz to k+1 voxel center
				// (1 - a) = 1-x0-1+gx = gx-x0 --> distance from current voxel center to gx
				// (1 - b) = 1-y0-1+gy = gy-y0 --> distance from current voxel center to gy
				ab = a * b;
				ab1 = a * ( 1 - b );
				a1b = ( 1 - a ) * b;
				a1b1 = ( 1 - a ) * ( 1 - b );

				// now, dummy variables
				a = ( 1 - c ) * voxel;
				c = (z1 - gz) * voxel;
				offset = z0*stepz + y0*stepy + x0; // initial offset

				// set the voxel data
				*(m_pData + offset) += ab*c;
				if(x1<limitx)
					*(m_pData + (offset + 1)) += a1b*c;

				if(y1<limity)
				{
					*(m_pData + (offset + stepy)) += ab1*c;
					if(x1<limitx)
						*(m_pData + (offset + stepy + 1)) += a1b1*c;
				}

				if(z1<limitz)
				{
					*(m_pData + (offset + stepz)) += ab*a;
					if(x1<limitx)
						*(m_pData + (offset + stepz + 1)) += a1b*a;
					if(y1<limity)
					{
						*(m_pData + (offset + stepz + stepy)) += ab1*a;
						if(x1<limitx)
							*(m_pData + (offset + stepz + stepy + 1)) += a1b1*a;
					}
				}
			}
		}
    	// Here, sending the end signal to main thread...
    	pthread_mutex_lock(((project_RealNE_data *) threadarg)->p_mutex_end);
    	(*(((project_RealNE_data *) threadarg)->p_nended))++; // increasing "nended" counter... for "wait" safety...
    	pthread_cond_signal(((project_RealNE_data *) threadarg)->p_cond_end); // Sending end signal
        pthread_mutex_unlock(((project_RealNE_data *) threadarg)->p_mutex_end);
    }
    pthread_exit(NULL);
}

// Mon made (29/05/2013)
// Projects over an pre-existing vlVolume a Macromolecule with NE and without moving
// the reference system in parallel. (Initialization routine)
// If fast=true, a fast method will be used instead of trilinear-interpolation (default=false)
void Macromolecule::project_RealNE_3BB2R_par_init(int nthreads, project_RealNE_data **p_threads_data, pthread_t **p_threads)
{
	bool debug = false;
	int i,j,rc;
	project_RealNE_data *threads_data;
	pthread_t *threads;
	pthread_mutex_t *p_mutex_map; // mutex pointer (to protect result map access)
	pthread_mutex_t *p_mutex_begin; // begin mutex pointer
	pthread_cond_t *p_cond_begin; // begin condition pointer
	pthread_mutex_t *p_mutex_end; // end mutex pointer
	pthread_cond_t *p_cond_end; // end condition pointer
	int *p_nended; // nended pointer

	// Mutexes and condition creation
	if( !(p_mutex_map = (pthread_mutex_t *) malloc(sizeof(pthread_mutex_t)) ) ||
			!(p_mutex_begin = (pthread_mutex_t *) malloc(sizeof(pthread_mutex_t)) ) ||
			!(p_mutex_end = (pthread_mutex_t *) malloc(sizeof(pthread_mutex_t)) ) ||
			!(p_cond_begin = (pthread_cond_t *) malloc(sizeof(pthread_cond_t)) ) ||
			!(p_cond_end = (pthread_cond_t *) malloc(sizeof(pthread_cond_t)) ) ||
			!(p_nended = (int *) malloc(sizeof(int)) ) )
	{
		fprintf(stderr,"Msg(project_RealNE_3BB2R_par_init): I'm sorry, mutex/condition/nended memory allocation failed!\nForcing exit!\n");
		exit(1);
	}

	// Initializing "nended"
	*p_nended = 0;

	// Allocating threads data
	if( !(threads_data = (project_RealNE_data *) malloc(sizeof(project_RealNE_data) * nthreads)) )
	{
		fprintf(stderr,"Msg(convoluteK_nopad_par_init): I'm sorry, thread memory allocation failed!\nForcing exit!\n");
		exit(1);
	}
	if( !(threads = (pthread_t *) malloc(sizeof(pthread_t) * nthreads)) )
	{
		fprintf(stderr,"Msg(convoluteK_nopad_par_init): I'm sorry, thread allocation failed!\nForcing exit!\n");
		exit(1);
	}
	*p_threads_data = threads_data;
	*p_threads = threads;

	// Initializing mutexes and conditions (only the first time)
	pthread_mutex_init(p_mutex_map, NULL); // mutex initialization
	pthread_mutex_init(p_mutex_begin, NULL); // mutex initialization
	pthread_cond_init(p_cond_begin, NULL); // condition initialization
	pthread_mutex_init(p_mutex_end, NULL); // mutex initialization
	pthread_cond_init(p_cond_end, NULL); // condition initialization


	for(i = 0; i < nthreads; i++)
	{
		threads_data[i].p_mutex_map = p_mutex_map; // passing by reference the map mutex
		threads_data[i].p_mutex_begin = p_mutex_begin; // passing by reference the begin mutex
		threads_data[i].p_cond_begin = p_cond_begin; // passing by reference the begin condition
		threads_data[i].p_mutex_end = p_mutex_end; // passing by reference the end mutex
		threads_data[i].p_cond_end = p_cond_end; // passing by reference the end condition
		threads_data[i].p_nended = p_nended; // passing by reference the nended variable
		threads_data[i].begin = true; // by default, thread is able to begin...
	}

	// Creating threads
	for(i = 0; i < nthreads; i++)
	{
		if(debug)
			fprintf(stderr,"Creating thread (convoluteK_nopad_par_init): %d\n", i);
		rc = pthread_create(&threads[i], NULL, project_RealNE_3BB2R_thread, (void *) &threads_data[i]);
		if (rc)
		{
			fprintf(stderr,"ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
	}
}

// Mon made (29/05/2013)
// Projects over an pre-existing vlVolume a Macromolecule with NE and without moving
// the reference system in parallel. (Controller routine)
// If fast=true, a fast method will be used instead of trilinear-interpolation (default=false)
void Macromolecule::project_RealNE_3BB2R_par( vlVolume *vlm, int nthreads, project_RealNE_data *threads_data)
{
	bool debug = false;
	int rc,i;
	void *status;
	vlPoint3f pos;
	vlm->getPosition(&pos);
	vlStep m_step = vlm->stepping();
	vlDim m_dimLimit = vlm->dim();
	int stepy = m_step.y();
	int stepz = m_step.z();
	int limitx = m_dimLimit.x();
	int limity = m_dimLimit.y();
	int limitz = m_dimLimit.z();

	int num_atoms = this->get_num_atoms();
	int natoms = num_atoms / nthreads; // Number of atoms scheduled per thread (to balance computational burden)

//	project_RealNE_data *threads_data = *p_threads_data;

	// Loading threads input data
	if(debug)
		fprintf(stderr,"Loading thread data:  num_atoms= %d   natoms= %d\n",num_atoms,natoms);
	for(i = 0; i < nthreads; i++)
	{
//		threads_data[i].iter = this; // it will be casted into Macromolecule when it be used by the threads...
		threads_data[i].iter = (pdbIter *) new pdbIter( this ); // create and cast iterator from current Macromolecule
		threads_data[i].vlm = vlm;
		threads_data[i].stepy = m_step.y();
		threads_data[i].stepz = m_step.z();
		threads_data[i].limitx = m_dimLimit.x();
		threads_data[i].limity = m_dimLimit.y();
		threads_data[i].limitz = m_dimLimit.z();
		threads_data[i].unitsx = vlm->units().x();
		threads_data[i].unitsy = vlm->units().y();
		threads_data[i].unitsz = vlm->units().z();
		threads_data[i].first = i*natoms; // First atom index
		threads_data[i].last = (i+1)*natoms; // Last atom index (not included)
		if(debug && i != nthreads-1)
			fprintf(stderr,"Loaded data for thread %d:  first= %d  last= %d\n",i,threads_data[i].first,threads_data[i].last);
	}
	threads_data[i-1].last = num_atoms; // Last atom index (not included)
	if(debug)
		fprintf(stderr,"Loaded data for thread %d:  i_initial= %d  i_final= %d\n",i-1,threads_data[nthreads-1].first,threads_data[nthreads-1].last);

	// Awaking all threads!
	if(debug)
		fprintf(stderr, "Awaking all threads! (project_RealNE_3BB2R_par)\n");
    pthread_mutex_lock(threads_data[0].p_mutex_end); // This prevents any end signal miss by controller thread.
    pthread_mutex_lock(threads_data[0].p_mutex_begin);
    for(int n=0; n<nthreads; n++)
    	threads_data[n].begin = false; // "re-signal" beginning...
	pthread_cond_broadcast(threads_data[0].p_cond_begin); // Send beginning  signal...
    pthread_mutex_unlock(threads_data[0].p_mutex_begin);
    // Waiting till all working threads finish.
	while( *(threads_data[0].p_nended) < nthreads )
	{
    	pthread_cond_wait(threads_data[0].p_cond_end, threads_data[0].p_mutex_end); // gets END signal from threads
    	if(debug)
    		fprintf(stderr,"Some thread has finished:  nended= %d\n",*(threads_data[0].p_nended));
	}

    pthread_mutex_unlock(threads_data[0].p_mutex_end);
    (*(threads_data[0].p_nended)) = 0; // reset "ended treads" counter...

    // Deleting created iterators...
	for(i = 0; i < nthreads; i++)
		delete( (pdbIter *)threads_data[i].iter );
}
#endif

void Macromolecule::writeVDW(char *filename, float voxel_size, int opt, ConventionNames opt2,Convention opt3,float probeRad,float probeEmax ,int withH,bool center_hyd,float rm,float f1,float f2)
{

	vlVolume *vol=createVDW(voxel_size, opt, opt2, opt3,probeRad,probeEmax ,withH,center_hyd,rm,f1,f2);
	FOPS::writeFile(vol,filename);
	free(vol);
}

vlVolume * Macromolecule::createVDW(float voxel_size, int opt, ConventionNames opt2,Convention opt3,float probeRad,float probeEmax ,int withH,bool center_hyd,float rm,float f1,float f2)
{
  PDB_Container * fath;
  Atom * a;
  Element *e;
  int i;
  int aminoacid,atom;
  char name[5];
  vlPoint3f origin,*position;
  vlPoint3ui *pos_aux;
  Tcoor posA;
  float max_dist;
  float max_dist_voxels;
  int limit_max[3],limit_min[3];
  int ii,jj,kk;
  float dist;
  int dist_down;
  float dist_pos;
  float value;
  int cont_mols;
  bool atom_found, primera_vez;

  primera_vez=true;

 // init_aminoacids(opt3,opt2);
 // nucleotids();
 // aa_pdb2pqr();
 // rename_residues();

  switch(opt)
  {
    case 0: init_vdw_c(probeRad,probeEmax,rm,f1,f2);
            break;
    case 1: init_vdw_h();
            break;
  }


  //for(i=0;i<MAX_VDW;i++)
  //printf("%d %f %f\n",i,VDW_C[0][i],i*VDW_STEP);
  //exit(1);



  max_dist=(MAX_VDW*VDW_STEP);
  max_dist_voxels=max_dist/(voxel_size*2.0);
  position= new vlPoint3f();
  vlVolume *vol_aux=createVolumeCentered(voxel_size,center_hyd);
  vlVolume *vol=FOPS::padVolume( vol_aux, vlDim(max_dist_voxels,max_dist_voxels,max_dist_voxels));
  vol_aux->~vlVolume();
  vol->clear(0);
  vol->getPosition(&origin);
 	//printf("Dimensiones del volumen: %d %d %d\nmaxima distancia:%f\n",vol->dim().x(),vol->dim().y(),vol->dim().z(),max_dist);




 	Macromolecule *mol[2];
	float mul[2];
	Conditions *conds= new Conditions();
	  Condition *condition= new Condition(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
  	condition->add(Nstr2);
  	condition->add(CAstr2);
  	condition->add(Cstr2);
  	condition->add(Ostr2);
  	condition->add(CBstr2);
  	condition->add(CGstr2);
  	condition->add(CDstr2);
  	condition->add(CEstr2);
  	condition->add(CZstr2);

//  	condition->add(" HA ");
//  	condition->add(" H  ");
  	conds->add(condition);

  mol[0]=this->select_cpy(conds);
  mul[0]=1.0;
	mol[1]=this->select_cpy(conds, false);
 	mul[1]=1.0;

 	//mol[0]->writePDB("kk1.pdb");
 	// mol[1]->writePDB("kk2.pdb");



 	for(cont_mols=0;cont_mols<2;cont_mols++)
 	{
 		mol[cont_mols]->initAll();

 		if(mol[cont_mols]->get_num_atoms()>0)
 		do
 	  	{
 	  		a=mol[cont_mols]->getCurrentAtom();
 	  		e=a->getElement();

 	  		fath=(PDB_Container*)a->getFather();
 	  		if(e->symbol!=H || withH )
 	  		{
 	  			//search Father
 	  			if(fath->getClass()==pdb_nucleotide || fath->getClass()==pdb_residue)
 	  			{

 	  				for(i=0;i<N_AMINO;i++)
 	  				{
 	  					if(strcmp(fath->getName(),AA[i].aa_name3)==0)
 	  					{
 	  						aminoacid=i;
 	  						//printf("%s en %s\n",fath->getName(),AA[i].aa_name3);
 	  						i=N_AMINO+10;
 	  						//getchar();
 	  					}
 	  				}
 	  				if(i==N_AMINO+1)
 	  				{
 	  					//aminoacid=GLY;
						aminoacid=-1;
						fprintf(stderr, "aa %s not found\n",fath->getName());
						primera_vez=true;
 	  				}
 	  			}
 	  			else
 	  			{
 	  				aminoacid=-1;
 	  			}


 	  			//search son
 	  			//fprintf(stderr,"Nombre %s\n",a->getName());
 	  			strcpy(name,a->getName());
 	  			name[4]='\0';
 	  			if(aminoacid>=0)
 	  			{
 	  				atom=21;
 	  				atom_found=false;
 	  				for(i=0;i<AA[aminoacid].natoms+3;i++)
 	  				{

 	  					if(strcmp(name,AA[aminoacid].atom[i].atom_name)==0)
 	  					{
 	  						atom=i;
 	  						atom_found=true;
 	  						atom=AA[aminoacid].atom[i].fullatom_type-1;
 	  						i=1000;
 	  						//fprintf(stderr,"found\n");
 	  					}
 	  				}
 	  			}

 	  			if(aminoacid<0 || !atom_found)
 	  			{
 	  				//if(!atom_found) {
 	  			     // fprintf(stderr,"Warning at createVDW Atom [%s] not found  %s %d\n",name,fath->getName(),a->getPdbSerial());
 	  				//getchar();
 	  				//}
 	  				atom_found=true;
 	  				atom=get_element_type(e);
 	  				if(atom==-1)
 	  				{
 	  					atom_found=false;
 	  				 	if(num_atom_type==54)
 	  				 	 	atom=21;
 	  				 	if(num_atom_type==45)
 	  				 	 	atom=14;

 	  				}
 	   	  		}

 	  			if(!atom_found)
 	  			{
 	  				if (primera_vez) {
 	  					fprintf(stderr,"CreateVDW skipped atoms\n");
 	  					primera_vez=false;
 	  				}

 	  				fprintf(stderr,"%s[%s] %d;", fath->getName(),name, a->getPdbSerial());
 	  			}
 	  			else
 	  			{
					//Operation
					//cout<<"Elemento: "<<(a->getElement())->symbol<<endl;
					//Exploracion de voxels vecinos
					a->getPosition(posA);
					//printf("%s %s %d %d\n",a->getPdbName(),fath->getName(),aminoacid,atom);

					position->x( ((posA[0]-origin.x())/vol->units().x()) );
					position->y( ((posA[1]-origin.y())/vol->units().y()) );
					position->z( ((posA[2]-origin.z())/vol->units().z()) );

					limit_max[0]=(int)(ceil(position->x())+max_dist);
					limit_max[1]=(int)(ceil(position->y())+max_dist);
					limit_max[2]=(int)(ceil(position->z())+max_dist);
					limit_min[0]=(int)(floor(position->x())-max_dist);
					limit_min[1]=(int)(floor(position->y())-max_dist);
					limit_min[2]=(int)(floor(position->z())-max_dist);

					if(limit_max[0]>=vol->dim().x())
						limit_max[0]=vol->dim().x()-1;
					if(limit_max[1]>=vol->dim().y())
						limit_max[1]=vol->dim().y()-1;
					if(limit_max[2]>=vol->dim().z())
						limit_max[2]=vol->dim().z()-1;

					if(limit_min[0]<0)
						limit_min[0]=0;
					if(limit_min[1]<0)
						limit_min[1]=0;
					if(limit_min[2]<0)
						limit_min[2]=0;

					//printf("limites de influencia\nx: %d/%d\ny: %d/%d\nz: %d/%d\n",
					//      limit_min[0],limit_max[0],
					//      limit_min[1],limit_max[1],
					//      limit_min[2],limit_max[2]);


					for(ii=limit_min[0];ii<limit_max[0];ii++)
						for(jj=limit_min[1];jj<limit_max[1];jj++)
							for(kk=limit_min[2];kk<limit_max[2];kk++)
							{
								dist=sqrt(((position->x()-ii)*(position->x()-ii))+
										((position->y()-jj)*(position->y()-jj))+
										((position->z()-kk)*(position->z()-kk)));
								dist*=vol->units().x()/VDW_STEP;
								dist_down=(int)floor(dist);
								if(dist_down<MAX_VDW-1)
								{
									dist_pos=dist-(float)dist_down;
									pos_aux = new vlPoint3ui( ii, jj, kk );
									vol->getVoxel( *pos_aux, value );
									switch(opt)
									{
										case 0:
											value +=(VDW_C[atom][dist_down]*(1-dist_pos)+VDW_C[atom][dist_down+1]*(dist_pos))*mul[cont_mols];
										break;
										//case 1: value +=(VDW_H[atom][dist_down]*(1-dist_pos)+VDW_H[atom][dist_down+1]*(dist_pos))*mul[cont_mols];
										//break;
									}
									vol->setVoxel( *pos_aux, value );
									delete( pos_aux );

								}
							}
							//printf("Voxels afectados:%d\n",cont);
 	  			}
 	  		}//if hyd
 	  	}while (mol[cont_mols]->nextAtom());
 	}//for end


//  FOPS::thresholdUp(vol, 3);

	//FOPS::threshold(vol,-3);
 // FOPS::thresholdUp(vol,15);
 delete conds;
 delete mol[0];
 delete mol[1];
 if (!primera_vez) fprintf(stderr,"\n");
 return vol;
}

// MPI version...
vlVolume * Macromolecule::createVDW(int start, int end, float voxel_size, int opt, ConventionNames opt2,Convention opt3,float probeRad,float probeEmax ,int withH,bool center_hyd,float rm,float f1,float f2)
{
	PDB_Container *fath;
	Atom *a;
	Element *e;
	int i;
	int aminoacid,atom;
	char name[5];
	vlPoint3f origin,*position;
	vlPoint3ui *pos_aux;
	Tcoor posA;
	float max_dist;
	float max_dist_voxels;
	int limit_max[3],limit_min[3];
	int ii,jj,kk;
	float dist;
	int dist_down;
	float dist_pos;
	float value;
	bool atom_found, primera_vez;

	primera_vez=true;

	// Initialize VdW look-up table...
	init_vdw_c(probeRad,probeEmax,rm,f1,f2);
	//for(i=0;i<MAX_VDW;i++)
	//printf("%d %f %f\n",i,VDW_C[0][i],i*VDW_STEP);
	//exit(1);

	// max_dist=(MAX_VDW*VDW_STEP)/2; // Mon: factor "2" reduces the effective VdW range to 15 A instead of 30 A.
	max_dist=20;
	max_dist_voxels=max_dist/(voxel_size*2.0);
	position= new vlPoint3f();
	vlVolume *vol_aux=createVolumeCentered(voxel_size,center_hyd);
	vlVolume *vol=FOPS::padVolume( vol_aux, vlDim(max_dist_voxels,max_dist_voxels,max_dist_voxels));
	vol_aux->~vlVolume();
	vol->clear(0.0);
	vol->getPosition(&origin);
	//printf("Dimensiones del volumen: %d %d %d\nmaxima distancia:%f\n",vol->dim().x(),vol->dim().y(),vol->dim().z(),max_dist);

	// Mon: Watch out! "this" should be a full copy...?
	// this->initAll();

	pdbIter *iter = new pdbIter( this ); // iter to screen atoms
	// pdbIter *iter = new pdbIter( this, true, true, true, true ); // iter to screen atoms

	for ( iter->pos_atom = start; iter->pos_atom < end; iter->next_atom() ) // screens atoms of selected chunk
	{
		a = iter->get_atom();
		e = a->getElement();

		if(e->symbol!=H || withH )
		{
			aminoacid = -1;  // Aminoacid not found by default

			// get father (aminoacid)
			fath = (PDB_Container*) a->getFather();
			if(fath->getClass()==pdb_nucleotide || fath->getClass()==pdb_residue)
			{
				// search father
				for(i=0;i<N_AMINO;i++)
				{
					if(strcmp(fath->getName(),AA[i].aa_name3)==0)
					{
						aminoacid = i;
						//printf("%s en %s\n",fath->getName(),AA[i].aa_name3);
						break; // Aminoacid found
					}
				}
			}

			if(aminoacid < 0) // Aminoacid not found
			{
				primera_vez=true;
				fprintf(stderr, "aa %s not found\n",fath->getName());
			}

			// Search son (atom) for current aminoacid
			strcpy(name,a->getName());
			name[4]='\0';
			//fprintf(stderr,"Nombre %s\n",a->getName());
			if(aminoacid >= 0) // aminoacid found
			{
				atom=21; // H (in non-sybil)
				atom_found=false;
				for(i=0;i<AA[aminoacid].natoms+3;i++)
				{
					if(strcmp(name,AA[aminoacid].atom[i].atom_name)==0)
					{
						atom=i;
						atom_found=true;
						atom=AA[aminoacid].atom[i].fullatom_type-1;
						break;
						//fprintf(stderr,"found\n");
					}
				}
			}

			if(aminoacid<0 || !atom_found) // atom not found for current aminoacid
			{
				//if(!atom_found) {
				// fprintf(stderr,"Warning at createVDW Atom [%s] not found  %s %d\n",name,fath->getName(),a->getPdbSerial());
				//getchar();
				//}
				atom_found=true;
				atom = get_element_type(e); // Mon: some kind of generic type???

				if(atom==-1) // if element type not found, then use H data
				{
					atom_found=false;
					if(num_atom_type==54)
						atom=21; // H not-sybil
					if(num_atom_type==45)
						atom=14; // H sybil
				}
			}

			if(!atom_found) // Atom finally not found!
			{
				if (primera_vez)
				{
					fprintf(stderr,"CreateVDW skipped atoms\n");
					primera_vez=false;
				}

				fprintf(stderr,"%s[%s] %d;", fath->getName(),name, a->getPdbSerial());
			}
			else
			{
				//Operation
				//cout<<"Elemento: "<<(a->getElement())->symbol<<endl;
				//Exploracion de voxels vecinos
				a->getPosition(posA);
				//printf("%s %s %d %d\n",a->getPdbName(),fath->getName(),aminoacid,atom);

				position->x( ((posA[0]-origin.x())/vol->units().x()) );
				position->y( ((posA[1]-origin.y())/vol->units().y()) );
				position->z( ((posA[2]-origin.z())/vol->units().z()) );

				limit_max[0]=(int)(ceil(position->x())+max_dist);
				limit_max[1]=(int)(ceil(position->y())+max_dist);
				limit_max[2]=(int)(ceil(position->z())+max_dist);
				limit_min[0]=(int)(floor(position->x())-max_dist);
				limit_min[1]=(int)(floor(position->y())-max_dist);
				limit_min[2]=(int)(floor(position->z())-max_dist);

				if(limit_max[0]>=vol->dim().x())
					limit_max[0]=vol->dim().x()-1;
				if(limit_max[1]>=vol->dim().y())
					limit_max[1]=vol->dim().y()-1;
				if(limit_max[2]>=vol->dim().z())
					limit_max[2]=vol->dim().z()-1;

				if(limit_min[0]<0)
					limit_min[0]=0;
				if(limit_min[1]<0)
					limit_min[1]=0;
				if(limit_min[2]<0)
					limit_min[2]=0;

				//printf("limites de influencia\nx: %d/%d\ny: %d/%d\nz: %d/%d\n",
				//      limit_min[0],limit_max[0],
				//      limit_min[1],limit_max[1],
				//      limit_min[2],limit_max[2]);

				for(ii=limit_min[0];ii<limit_max[0];ii++)
					for(jj=limit_min[1];jj<limit_max[1];jj++)
						for(kk=limit_min[2];kk<limit_max[2];kk++)
						{
							dist=sqrt(((position->x()-ii)*(position->x()-ii))+
									((position->y()-jj)*(position->y()-jj))+
									((position->z()-kk)*(position->z()-kk)));
							dist*=vol->units().x()/VDW_STEP;
							dist_down=(int)floor(dist);
							if(dist_down<MAX_VDW-1)
							{
								dist_pos=dist-(float)dist_down;
								pos_aux = new vlPoint3ui( ii, jj, kk );
								vol->getVoxel( *pos_aux, value );
								switch(opt)
								{
								case 0:
									//value +=(VDW_C[atom][dist_down]*(1-dist_pos)+VDW_C[atom][dist_down+1]*(dist_pos))*mul[cont_mols];
									value +=(VDW_C[atom][dist_down]*(1-dist_pos)+VDW_C[atom][dist_down+1]*(dist_pos));
									break;
									//case 1: value +=(VDW_H[atom][dist_down]*(1-dist_pos)+VDW_H[atom][dist_down+1]*(dist_pos))*mul[cont_mols];
									//break;
								}
								vol->setVoxel( *pos_aux, value );
								delete( pos_aux );

							}
						}
				//printf("Voxels afectados:%d\n",cont);
			}
		}//if hyd

	}
	delete(iter);
	if (!primera_vez)
		fprintf(stderr,"\n");
	return vol;
}




void Macromolecule::writeASA_c(char *filename, float voxel_size,int withH)
{

  vlVolume *vol=createVolumeCentered(voxel_size);
  //vol=FOPS::padVolume(vol,vlDim(5,5,5),true);
 	projectASA(vol,withH);
	FOPS::mul(vol,-1.0);
  FOPS::writeFile(vol,filename);
	printf("total mass=%f\n",FOPS::calc_total_posneg(vol));
	free (vol);
}

void Macromolecule::writeCharge(char *filename, float voxel_size,int withH)
{

  vlVolume *vol=createVolumeCentered(voxel_size);
 	projectCharge(vol,withH);

  FOPS::writeFile(vol,filename);
	free(vol);
}

void Macromolecule::projectASA(vlVolume *vol,int withH)
{
  vlPoint3f pos;
  Tcoor atomCoor;
	Atom *at;
	bool end=true;
	float weight;
  vol->clear(0.0);
  vol->getPosition(&pos);
  float gx, gy, gz;

Tcoor volcenter;
	FOPS::center_vol(vol,volcenter);

	initAll();
  while ( end != false )
  {
    /* Obtengo un atomo */
    at = getCurrentAtom();
    /* Obtencion de la posicion real del atomo en el PDb en amstrongs */
    at->getPosition(atomCoor);

    /* Calculo de la posicion relativa del atomo en el volumen en numero de voxels */
    gx = ((atomCoor[0]- pos.x()) / vol->units().x());
    gy = ((atomCoor[1]- pos.y()) / vol->units().y());
    gz = ((atomCoor[2]- pos.z()) / vol->units().z());

    /* Peso atomico del atomo */
    //if(at->getPdbocc()>15.0)
    //{
    	weight = ( float )( at->getPdbocc()/100.0);
			//weight=1.0;

    if(weight>0 )
    {
  			vol->interpolateVoxel(vlPoint3f(gx,gy,gz), weight);

    	  //vol->setVoxel(vlPoint3ui(gx,gy,gz), weight);
    	//}
    	//return;
    }
    end = nextAtom();
  }



}

void Macromolecule::projectCharge(vlVolume *vol,int withH)
{
  vlPoint3f pos;
  Tcoor atomCoor;
	Atom *at;
	bool end=true;
	float weight;
  vol->clear(0.0);
  vol->getPosition(&pos);
  float gx, gy, gz;

	initAll();
  while ( end != false )
  {
    /* Obtengo un atomo */
    at = getCurrentAtom();
    /* Obtencion de la posicion real del atomo en el PDb en amstrongs */
    if((at->getElement())->symbol!=H || withH)
    {
    	at->getPosition(atomCoor);

    	/* Calculo de la posicion relativa del atomo en el volumen en numero de voxels */
    	gx = ((atomCoor[0]- pos.x()) / vol->units().x()) ;
    	gy = ((atomCoor[1]- pos.y()) / vol->units().y()) ;
    	gz = ((atomCoor[2]- pos.z()) / vol->units().z()) ;

    	/* Peso atomico del atomo */
    	//if(at->getPdbocc()>15.0)
    	//{
    	weight = ( float )( at->getPdbfact());
			//weight=1.0;
  			vol->interpolateVoxel(vlPoint3f(gx,gy,gz), weight);
    }
    end = nextAtom();
  }


}



void Macromolecule::writeELE(char *filename, float voxel_size, ConventionNames opt2,Convention opt3, float threshold,int type,int withH)
{
	vlVolume *vol= createELE(voxel_size, threshold, type,withH);
	FOPS::writeFile(vol,filename);
	free(vol);
}
void Macromolecule::writeELE(char *filename, float voxel_size, int type, bool withH)
{
	vlVolume *vol= createELE_file(voxel_size, type,  withH);
	FOPS::writeFile(vol,filename);
	free(vol);
}

vlVolume *Macromolecule::createELE(float voxel_size, float threshold,int type,int withH)
{
  PDB_Container * fath;
  Atom * a;
  Element *e;
  int i;
  int aminoacid;
  char name[5];
  vlPoint3f origin,*position;
  vlPoint3ui *pos_aux;
  Tcoor posA;
  float max_dist;
  float max_dist_voxels;
  int limit_max[3],limit_min[3];
  int ii,jj,kk;
  float dist;
  float value;
  bool atom_found;
//	float in300=(332.0)/((MAX_ELE*ELE_STEP*MAX_ELE*ELE_STEP*4.0));
	float radius;
//	float e0=8.8542e-12;
	float ron=8.0,roff=30.0;
//	float elec_B=78.4+8.5525;
//	float er;
	float charge;

    // Change default charges....
    aa_iupac2pqr();
    init_ele();
  	change_charges_withoutH();

  	//Dimensiones del mapa
  	max_dist=(300.0*ELE_STEP);
  	max_dist_voxels=max_dist/(voxel_size);

  	position= new vlPoint3f();
  	vlVolume *vol_aux=createVolumeCentered(voxel_size);
  	vlVolume *vol=FOPS::padVolume( vol_aux, vlDim(max_dist_voxels/2.0,max_dist_voxels/2.0,max_dist_voxels/2.0));

  	//Dimensiones de la busqueda
  	max_dist=roff+1.0;//(MAX_ELE*ELE_STEP);


  	//vlVolume *vol=FOPS::padVolume( vol_aux, vlDim((10.0/voxel_size),(10.0/voxel_size),(10.0/voxel_size)));
  	delete vol_aux;
  	vol->clear(0);
  	vol->getPosition(&origin);

 		initAll();

 		if(get_num_atoms()>0)
 	  	do
 	  	{
 	  		a=getCurrentAtom();
 	  		e=a->getElement();
 	  		fath=(PDB_Container*)a->getFather();
 	  		if( (e->symbol!=H || withH) && fath->getClass()!= pdb_smol)
 	  		{
 	  			//search Father
 	  			if(fath->getClass()==pdb_nucleotide || fath->getClass()==pdb_residue)
 	  			{

 	  				for(i=0;i<N_AMINO;i++)
 	  					if(strcmp(fath->getName(),AA[i].aa_name3)==0)
 	  					{
 	  						aminoacid=i; break;
 	  						//printf("%s en %s\n",res->getName(),AA[i].aa_name3);
 	  						//i=N_AMINO+10;
 	  					}

 	  				if(i==N_AMINO+1)
 	  				{	aminoacid=GLY;
						aminoacid=-1;
 	  					printf("Alert: Aminoacid not found %s. The atoms will be evaluated by its element.\n",fath->getName());
 	  				}
 	  			}
 	  			else
 	  			{
 	  				aminoacid=-1;
 	  			}


 	  			//search son
 	  			//fprintf(stderr,"Nombre %s\n",a->getName());
 	  			strcpy(name,a->getName());
 	  			name[4]='\0';
 	  			if(aminoacid>=0)
 	  			{
 	  				atom_found=false;
 	  				for(i=0;i<AA[aminoacid].natoms+3;i++)
 	  				{
 	  					if(strcmp(name,AA[aminoacid].atom[i].atom_name)==0)
 	  					{
 	  						radius= atom_types[AA[aminoacid].atom[i].fullatom_type-1].vdw;
 	  						charge=AA[aminoacid].atom[i].charge;
 	  						//fprintf(stderr,"RES=%s Atom=%s charge=%f\n",AA[aminoacid].aa_name3,AA[aminoacid].atom[i].atom_name,charge);
 	  						atom_found=true;
 	  						break;
 	  					}
 	  				}

 	  			}
 	  			else
 	  			{
 	  				atom_found=true;
 	  				// atom not found
 	  				charge=0;
 	  				radius=1;
 	  			}

 	  			if(!atom_found)
 	  			{
 	  				fprintf(stderr,"Alert: Atom type not found. Treated as Hydrogen: %s %d\n",name,a->getPdbSerial());
 	  			}

 	  			//Operation
 	  			//cout<<"Elemento: "<<(a->getElement())->symbol<<endl;
 	  			//Exploracion de voxels vecinos
 	  			a->getPosition(posA);

 	  			position->x( ((posA[0]-origin.x())/vol->units().x()) );
 	  			position->y( ((posA[1]-origin.y())/vol->units().y()) );
 	  			position->z( ((posA[2]-origin.z())/vol->units().z()) );

 	  			limit_max[0]=(int)(ceil(position->x())+max_dist);
 	  			limit_max[1]=(int)(ceil(position->y())+max_dist);
 	  			limit_max[2]=(int)(ceil(position->z())+max_dist);
 	  			limit_min[0]=(int)(floor(position->x())-max_dist);
 	  			limit_min[1]=(int)(floor(position->y())-max_dist);
 	  			limit_min[2]=(int)(floor(position->z())-max_dist);

 	  			if(limit_max[0]>=vol->dim().x())
 	  				limit_max[0]=vol->dim().x()-1;
 	  			if(limit_max[1]>=vol->dim().y())
 	  				limit_max[1]=vol->dim().y()-1;
				if(limit_max[2]>=vol->dim().z())
					limit_max[2]=vol->dim().z()-1;

				if(limit_min[0]<0)
					limit_min[0]=0;
				if(limit_min[1]<0)
					limit_min[1]=0;
				if(limit_min[2]<0)
					limit_min[2]=0;

				//printf("limites de influencia\nx: %d/%d\ny: %d/%d\nz: %d/%d\n",
				//      limit_min[0],limit_max[0],
				//      limit_min[1],limit_max[1],
				//      limit_min[2],limit_max[2]);


				for(ii=limit_min[0];ii<limit_max[0];ii++)
					for(jj=limit_min[1];jj<limit_max[1];jj++)
						for(kk=limit_min[2];kk<limit_max[2];kk++)
						{
							dist=sqrt(((position->x()-ii)*(position->x()-ii))+
									((position->y()-jj)*(position->y()-jj))+
									((position->z()-kk)*(position->z()-kk)));
							dist*=vol->units().x();

							if(dist<roff )
							{
								pos_aux = new vlPoint3ui( ii, jj, kk );
								vol->getVoxel( *pos_aux, value );

								value += autodock_ele (dist, ron, roff,radius,charge, type );

								vol->setVoxel( *pos_aux, value );
								delete( pos_aux );

							}

						}
      		            //printf("Voxels afectados:%d\n",cont);
 	  		}//if hyd
 	  	}while (nextAtom());

//  FOPS::thresholdUp(vol, 3);

	//FOPS::threshold(vol,-3);
 // FOPS::thresholdUp(vol,15);
  return vol;

}

vlVolume *Macromolecule::createELE(int start, int end, float voxel_size, float threshold,int type,int withH)
{
	PDB_Container * fath;
	Atom * a;
	Element *e;
	int i;
	int aminoacid;
	char name[5];
	vlPoint3f origin,*position;
	vlPoint3ui *pos_aux;
	Tcoor posA;
	float max_dist;
	float max_dist_voxels;
	int limit_max[3],limit_min[3];
	int ii,jj,kk;
	float dist;
	float value;
	bool atom_found;
	float radius;
	float ron=8.0,roff=30.0;
	float charge;

	// Change default charges....
	aa_iupac2pqr();
	init_ele();
	change_charges_withoutH();

	//Dimensiones del mapa
	max_dist=(300.0*ELE_STEP);
	max_dist_voxels=max_dist/(voxel_size);

	position= new vlPoint3f();
	vlVolume *vol_aux=createVolumeCentered(voxel_size);
	vlVolume *vol=FOPS::padVolume( vol_aux, vlDim(max_dist_voxels/2.0,max_dist_voxels/2.0,max_dist_voxels/2.0));

	//Dimensiones de la busqueda
	max_dist=roff+1.0;//(MAX_ELE*ELE_STEP);


	//vlVolume *vol=FOPS::padVolume( vol_aux, vlDim((10.0/voxel_size),(10.0/voxel_size),(10.0/voxel_size)));
	delete vol_aux;
	vol->clear(0);
	vol->getPosition(&origin);

	initAll();

	pdbIter *iter = new pdbIter( this ); // iter to screen atoms
	// pdbIter *iter = new pdbIter( this, true, true, true, true ); // iter to screen atoms

	for ( iter->pos_atom = start; iter->pos_atom < end; iter->next_atom() ) // screens atoms of selected chunk
	{
		a = iter->get_atom();
		e = a->getElement();

		fath = (PDB_Container*) a->getFather();
		if( (e->symbol!=H || withH) && fath->getClass()!= pdb_smol)
		{
			//search Father
			if(fath->getClass()==pdb_nucleotide || fath->getClass()==pdb_residue)
			{

				for(i=0;i<N_AMINO;i++)
					if(strcmp(fath->getName(),AA[i].aa_name3)==0)
					{
						aminoacid=i;
						break;
						//printf("%s en %s\n",res->getName(),AA[i].aa_name3);
						//i=N_AMINO+10;
					}

				if(i==N_AMINO+1)
				{
					// aminoacid=GLY;
					aminoacid=-1;
					printf("Alert: Aminoacid not found %s. The atoms will be evaluated by its element.\n",fath->getName());
				}
			}
			else
			{
				aminoacid=-1;
			}


			//search son
			//fprintf(stderr,"Nombre %s\n",a->getName());
			strcpy(name,a->getName());
			name[4]='\0';
			if(aminoacid>=0)
			{
				atom_found=false;
				for(i=0;i<AA[aminoacid].natoms+3;i++)
				{
					if(strcmp(name,AA[aminoacid].atom[i].atom_name)==0)
					{
						radius= atom_types[AA[aminoacid].atom[i].fullatom_type-1].vdw;
						charge=AA[aminoacid].atom[i].charge;
						//fprintf(stderr,"RES=%s Atom=%s charge=%f\n",AA[aminoacid].aa_name3,AA[aminoacid].atom[i].atom_name,charge);
						atom_found=true;
						break;
					}
				}

			}
			else
			{
				atom_found=true;
				// atom not found
				charge=0;
				radius=1;
			}

			if(!atom_found)
			{
				fprintf(stderr,"Alert: Atom type not found. Treated as Hydrogen: %s %d\n",name,a->getPdbSerial());
			}

			//Operation
			//cout<<"Elemento: "<<(a->getElement())->symbol<<endl;
			//Exploracion de voxels vecinos
			a->getPosition(posA);

			position->x( ((posA[0]-origin.x())/vol->units().x()) );
			position->y( ((posA[1]-origin.y())/vol->units().y()) );
			position->z( ((posA[2]-origin.z())/vol->units().z()) );

			limit_max[0]=(int)(ceil(position->x())+max_dist);
			limit_max[1]=(int)(ceil(position->y())+max_dist);
			limit_max[2]=(int)(ceil(position->z())+max_dist);
			limit_min[0]=(int)(floor(position->x())-max_dist);
			limit_min[1]=(int)(floor(position->y())-max_dist);
			limit_min[2]=(int)(floor(position->z())-max_dist);

			if(limit_max[0]>=vol->dim().x())
				limit_max[0]=vol->dim().x()-1;
			if(limit_max[1]>=vol->dim().y())
				limit_max[1]=vol->dim().y()-1;
			if(limit_max[2]>=vol->dim().z())
				limit_max[2]=vol->dim().z()-1;

			if(limit_min[0]<0)
				limit_min[0]=0;
			if(limit_min[1]<0)
				limit_min[1]=0;
			if(limit_min[2]<0)
				limit_min[2]=0;

			for(ii=limit_min[0];ii<limit_max[0];ii++)
				for(jj=limit_min[1];jj<limit_max[1];jj++)
					for(kk=limit_min[2];kk<limit_max[2];kk++)
					{
						dist=sqrt(((position->x()-ii)*(position->x()-ii))+
								((position->y()-jj)*(position->y()-jj))+
								((position->z()-kk)*(position->z()-kk)));
						dist*=vol->units().x();

						if(dist<roff )
						{
							pos_aux = new vlPoint3ui( ii, jj, kk );
							vol->getVoxel( *pos_aux, value );

							value += autodock_ele (dist, ron, roff,radius,charge, type );

							vol->setVoxel( *pos_aux, value );
							delete( pos_aux );

						}

					}
			//printf("Voxels afectados:%d\n",cont);
		}//if hyd
	}

	return vol;

}


// Mon: Warning, this function has not been parallelized yet...
vlVolume *Macromolecule::createELE_file(float voxel_size, int type, int withH)
{
  Atom * a;
  Element *e;
  vlPoint3f origin,*position;
  vlPoint3ui *pos_aux;
  Tcoor posA;
  float max_dist;
  float max_dist_voxels;
  int limit_max[3],limit_min[3];
  int ii,jj,kk;
  float dist;
  float value;
	float radius;
	float ron=8.0,roff=30.0;
	float charge;

  	//Dimensiones del mapa
  	max_dist=(300.0*ELE_STEP);
  	max_dist_voxels=max_dist/(voxel_size);

  	position= new vlPoint3f();
  	vlVolume *vol_aux=createVolumeCentered(voxel_size);
  	vlVolume *vol=FOPS::padVolume( vol_aux, vlDim(max_dist_voxels/2.0,max_dist_voxels/2.0,max_dist_voxels/2.0));

  	//Dimensiones de la busqueda
  	max_dist=roff+1.0;//(MAX_ELE*ELE_STEP);


  	//vlVolume *vol=FOPS::padVolume( vol_aux, vlDim((10.0/voxel_size),(10.0/voxel_size),(10.0/voxel_size)));
  	delete vol_aux;
  	vol->clear(0);
  	vol->getPosition(&origin);

 		initAll();

 		if(get_num_atoms()>0)
 	  	do
 	  	{
 	  		a=getCurrentAtom();
 	  		e=a->getElement();
 	  		if( (e->symbol!=H || withH) )
 	  		{


 	  			//Operation
 	  			//cout<<"Elemento: "<<(a->getElement())->symbol<<endl;
 	  			//Exploracion de voxels vecinos
 	  			a->getPosition(posA);
 	  			charge=a->getPdbfact();
 	  			radius=a->getPdbocc();
 	  			if(radius==0)
 	  			{
 	  				fprintf(stderr,"Alert: Radius 0.0??? changed to 1.2\n");
 	  				radius=1.2;
 	  			}
 	  			position->x( ((posA[0]-origin.x())/vol->units().x()) );
 	  			position->y( ((posA[1]-origin.y())/vol->units().y()) );
 	  			position->z( ((posA[2]-origin.z())/vol->units().z()) );

 	  			limit_max[0]=(int)(ceil(position->x())+max_dist);
 	  			limit_max[1]=(int)(ceil(position->y())+max_dist);
 	  			limit_max[2]=(int)(ceil(position->z())+max_dist);
 	  			limit_min[0]=(int)(floor(position->x())-max_dist);
 	  			limit_min[1]=(int)(floor(position->y())-max_dist);
 	  			limit_min[2]=(int)(floor(position->z())-max_dist);

 	  			if(limit_max[0]>=vol->dim().x())
 	  				limit_max[0]=vol->dim().x()-1;
 	  			if(limit_max[1]>=vol->dim().y())
 	  				limit_max[1]=vol->dim().y()-1;
				if(limit_max[2]>=vol->dim().z())
					limit_max[2]=vol->dim().z()-1;

				if(limit_min[0]<0)
					limit_min[0]=0;
				if(limit_min[1]<0)
					limit_min[1]=0;
				if(limit_min[2]<0)
					limit_min[2]=0;


				for(ii=limit_min[0];ii<limit_max[0];ii++)
					for(jj=limit_min[1];jj<limit_max[1];jj++)
						for(kk=limit_min[2];kk<limit_max[2];kk++)
						{
							dist=sqrt(((position->x()-ii)*(position->x()-ii))+
									((position->y()-jj)*(position->y()-jj))+
									((position->z()-kk)*(position->z()-kk)));
							dist*=vol->units().x();
							if(dist<roff )
							{
								pos_aux = new vlPoint3ui( ii, jj, kk );
								vol->getVoxel( *pos_aux, value );
								value += autodock_ele (dist, ron, roff,radius,charge, type );
								vol->setVoxel( *pos_aux, value );
								delete( pos_aux );

							}

						}
      		            //printf("Voxels afectados:%d\n",cont);
 	  		}//if hyd
 	  	}while (nextAtom());


  return vol;

}

int roundFloat( float x )
{

  int aux1, aux2;
  float dist1, dist2;
  aux1 = ( int )ceil( x );
  aux2 = ( int )floor( x );
  dist1 = x - aux1;
  if ( dist1 < 0 ) dist1 *= -1;
  dist2 = x - aux2;
  if ( dist2 < 0 ) dist2 *= -1;

  if ( dist2 < dist1 ) return aux2;
  else
    return aux1;
}


vlVolume * Macromolecule::project_radiusVDW(float unit, int withH,bool center_hyd, float factor_radius)

{
  bool end = true;
  Tcoor atomCoor;

  Atom * at;
  Element *e;
  vlVolume * vlm;
  vlPoint3ui * position;
  float vdwr, value;

  float shiftX, shiftY, shiftZ;
  float centerX, centerY, centerZ, ii, jj, kk;
  float gx, gy, gz;

  int i, j, k;


  /* create paded volume size + kernel size */
  vlm = createVolume_pad( unit, 3 );
	Tcoor volcenter;
  FOPS::center_vol(vlm,volcenter);
  //printf("  grid size %f %f\n",unit, vlm->units().x());
  //printf(" Protein box %d %d %d \n",box[0],box[1],box[2]);
  //printf(" Protein center %f %f %f\n",macromoleculeDimension.geometricCenter[0],
  //       macromoleculeDimension.geometricCenter[1],macromoleculeDimension.geometricCenter[2]);

  /* shift from the PDB geometric center to center of the map */
  shiftX = volcenter[0] - macromoleculeDimension.geometricCenter[0] / vlm->units().x();
  shiftY = volcenter[1] - macromoleculeDimension.geometricCenter[1] / vlm->units().y();
  shiftZ = volcenter[2] - macromoleculeDimension.geometricCenter[2] / vlm->units().z();

  // get geometric center of the map in amstrons
  centerX = ( float ) - unit * volcenter[0] + macromoleculeDimension.geometricCenter[0];
  centerY = ( float ) - unit * volcenter[1] + macromoleculeDimension.geometricCenter[1];
  centerZ = ( float ) - unit * volcenter[2] + macromoleculeDimension.geometricCenter[2];
  vlm->setPosition( vlPoint3f( centerX, centerY, centerZ ) );
  //printf(" Protein center %f %f %f    %f %f %f \n", centerX,centerY,centerX,shiftX, shiftY, shiftZ);


  /* gaussian expansion of pdb */
  initAll();
  position = new vlPoint3ui( 0, 0, 0 );

  value=1;
  while ( end != false )
  {
	  /* get atom */
	  at = getCurrentAtom();
	  /* getelectronic density */
	  e=at->getElement();
	  // 			fath=(PDB_Container*)a->getFather();
	  if(e->symbol!=H || withH )
	  {

		  // vdwr=( at->getElement() )->vdw*factor_radius;
		  vdwr=e->vdw*factor_radius;
		  //vdwr=1.5*2.0*factor_radius;
		  /* get position amstrons */
		  at->getPosition(atomCoor);

		  /* get relative map position */
		  gx = ( float )( atomCoor[0] / unit + shiftX );
		  gy = ( float )( atomCoor[1] / unit + shiftY );
		  gz = ( float )( atomCoor[2] / unit + shiftZ );

		  /* calculate grid pos around point <vdwr*/
		  float li;
		  li=roundFloat(vdwr/unit)*unit;

		  for(ii=-li;ii<=li;ii+=unit)
			  for(jj=-li;jj<=li;jj+=unit)
				  for(kk=-li;kk<=li;kk+=unit)


					  //if (dist=sqrt(ii*ii+jj*jj+kk*kk)<=vdwr)  	 {
					  if ((ii*ii+jj*jj+kk*kk)<=vdwr*vdwr)
					  {


						  i = ( int )roundFloat( gx + ii/unit);
						  j = ( int )roundFloat( gy + jj/unit);
						  k = ( int )roundFloat( gz + kk/unit);

						  // printf("%s en %f %d %d %d  %f %f %f\n",at->getName(),vdwr, i, j, k, ii, jj, kk);

						  //value=0.0001+vdwr-dist;
						  value=1.0;
						  /* set position */
						  position->x( i );
						  position->y( j );
						  position->z( k );
						  vlm->setVoxel( * position, value );
					  }
	  }
	  end = nextAtom();
  }
  return vlm;
}


// // Mon's version...
//vlVolume * Macromolecule::project_radiusVDW(float voxel_size, int withH,bool center_hyd, float factor_radius)
//{
//  PDB_Container *fath;
//  Atom * a;
//  Element *e;
//   vlPoint3f origin, *position;
//	int cont_mols;
//	int i, aminoacid, atom;
//	char name[5];
//	Tcoor posA;
//  int limit_max[3],limit_min[3];
//  float vdwr, value=1,dist;
//  int ii,jj,kk;
//  bool atom_found;
//  int cont=0;
//
//  vlPoint3ui *pos_aux;
//  // init_aminoacids(opt3,opt2);
//  //aa_iupac2pqr();
//  // rename_residues();
//  vlVolume *vol=createVolumeCentered(voxel_size,center_hyd);
//  vol->clear(0);
//  vol->getPosition(&origin);
//  position= new vlPoint3f();
//
//  if (factor_radius<=0) factor_radius=1;
//
// 	Macromolecule *mol[2];
//	float mul[2];
//	Conditions *conds= new Conditions();
//	  Condition *condition= new Condition(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);
//  	condition->add(" C  ");
//  	condition->add(" CA ");
//  	condition->add(" N  ");
//  	condition->add(" O  ");
//  	condition->add(" CB ");
//  	condition->add(" HA ");
//  	condition->add(" H  ");
//  	conds->add(condition);
//
//  mol[0]=this->select_cpy(conds);
//  mul[0]=1.0;
//	mol[1]=this->select_cpy(conds, false);
// 	mul[1]=1.0;
//
// 	//mol[0]->writePDB("kk1.pdb");
// 	//mol[1]->writePDB("kk2.pdb");
//
//
// 	for(cont_mols=0;cont_mols<2;cont_mols++)
// 	{
// 		mol[cont_mols]->initAll();
//
// 		if(mol[cont_mols]->get_num_atoms()>0)
// 		do
// 		{
// 			a=mol[cont_mols]->getCurrentAtom();
// 			e=a->getElement();
// 			fath=(PDB_Container*)a->getFather();
// 			if(e->symbol!=H || withH )
// 			{
//
// 				//fprintf(stderr,"Nombre %s\n",a->getName());
// 			 	//strcpy(name,a->getName());
// 				vdwr=e->vdw*factor_radius;
//
// 				//Exploracion de voxels vecinos
// 				a->getPosition(posA);
//
// 				position->x( ((posA[0]-origin.x())/vol->units().x()) );
// 				position->y( ((posA[1]-origin.y())/vol->units().y()) );
// 				position->z( ((posA[2]-origin.z())/vol->units().z()) );
//
// 				limit_max[0]=(int)(ceil(position->x())+vdwr);
// 				limit_max[1]=(int)(ceil(position->y())+vdwr);
// 				limit_max[2]=(int)(ceil(position->z())+vdwr);
// 				limit_min[0]=(int)(floor(position->x())-vdwr);
// 				limit_min[1]=(int)(floor(position->y())-vdwr);
// 				limit_min[2]=(int)(floor(position->z())-vdwr);
//
// 				if(limit_max[0]>=vol->dim().x())
// 					limit_max[0]=vol->dim().x()-1;
// 				if(limit_max[1]>=vol->dim().y())
// 					limit_max[1]=vol->dim().y()-1;
// 				if(limit_max[2]>=vol->dim().z())
// 					limit_max[2]=vol->dim().z()-1;
//
// 				if(limit_min[0]<0)
// 					limit_min[0]=0;
// 				if(limit_min[1]<0)
// 					limit_min[1]=0;
// 				if(limit_min[2]<0)
// 					limit_min[2]=0;
//
// 				for(ii=limit_min[0];ii<limit_max[0];ii++)
// 					for(jj=limit_min[1];jj<limit_max[1];jj++)
// 						for(kk=limit_min[2];kk<limit_max[2];kk++)
// 						{
// 							dist=sqrt(((position->x()-ii)*(position->x()-ii))+
// 									((position->y()-jj)*(position->y()-jj))+
// 									((position->z()-kk)*(position->z()-kk)));
// 							dist*=vol->units().x();
//
// 							if(dist<vdwr)
// 							{
// 								pos_aux = new vlPoint3ui( ii, jj, kk );
// 								vol->setVoxel( *pos_aux, value );
// 								delete (pos_aux);
// 								//cont++;
// 							}
// 						}
// 				//printf("Voxels afectados:%d\n",cont);
//				}
// 			}while (mol[cont_mols]->nextAtom());
// 		}
//
//
//
//
//  delete (position);
//  delete conds;
//  delete mol[0];
//  delete mol[1];
//  return vol;
//
//}

vlVolume * Macromolecule::project_Hydrophoby(float voxel_size,  ConventionNames opt2,Convention opt3,float ideal, float max,float sigma,int withH,bool center_hyd)
{
	PDB_Container *fath;
	Atom * a;
	Element *e;
	vlPoint3f origin, *position;
	int i, aminoacid, atom;
	char name[5];
	Tcoor posA;
	int limit_max[3],limit_min[3];
	float dist;
	int ii,jj,kk;
	bool atom_found;
	float max_voxel;
	float value,value2;
	int polarity;
	float in_ideal,in_max;
	float vdwr;

	vlPoint3ui *pos_aux;

	vlVolume *vol=createVolumeCentered(voxel_size,center_hyd);
	vol->clear(0);
	vol->getPosition(&origin);
	position= new vlPoint3f();
	max_voxel=(max/vol->units().x())+1.0;

	Macromolecule *mol[2];
	float mul[2];
	Conditions *conds= new Conditions();
	Condition *condition= new Condition(-1,-1,-1,-1,-1,-1,-1,-1,-1,-1);


	condition->add(Nstr2);
  	condition->add(CAstr2);
  	condition->add(Cstr2);
  	condition->add(Ostr2);
  	condition->add(CBstr2);
  	condition->add(HAstr2);
  	condition->add(Hstr2);
	conds->add(condition);

	mol[0]=this->select_cpy(conds);
	mul[0]=1.0;
	mol[1]=this->select_cpy(conds, false);
	mul[1]=1.0;

	//mol[0]->writePDB("kk1.pdb");
	//mol[1]->writePDB("kk2.pdb");


	// 	for(cont_mols=0;cont_mols<2;cont_mols++)
	// 	{
	// 		mol[cont_mols]->initAll();
	initAll();

	// 		if(mol[cont_mols]->get_num_atoms()>0)
	if(get_num_atoms()>0)
		do
		{
			// 			a=mol[cont_mols]->getCurrentAtom();
			a=getCurrentAtom();
			e=a->getElement();
			fath=(PDB_Container*)a->getFather();
			if (a->getPdbocc()>=5.0)
			{
				if(e->symbol!=H || withH )
				{
					//search Father
					if(fath->getClass()==pdb_nucleotide || fath->getClass()==pdb_residue)
					{

						for(i=0;i<N_AMINO;i++)
						{
							if(strcmp(fath->getName(),AA[i].aa_name3)==0)
							{
								aminoacid=i;
								//printf("%s en %s\n",res->getName(),AA[i].aa_name3);
								i=N_AMINO+10;
							}
						}
						if(i==N_AMINO+1)
						{	aminoacid=GLY;
						printf("Alert: Aminoacid not found %s\n",fath->getName());
						}
					}
					else
					{
						aminoacid=-1;
					}

					//search son
					//fprintf(stderr,"Nombre %s\n",a->getName());
					strcpy(name,a->getName());
					name[4]='\0';
					polarity=true;
					if(aminoacid>0)
					{
						atom=21;
						atom_found=false;
						for(i=0;i<AA[aminoacid].natoms+3;i++)
						{
							if(strcmp(name,AA[aminoacid].atom[i].atom_name)==0)
							{
								atom=i;
								atom_found=true;
								atom=AA[aminoacid].atom[i].fullatom_type-1;
								polarity=atom_types[atom].polar;
								vdwr=atom_types[atom].vdw+ideal;
								i=1000;
							}
						}
					}

					if(!atom_found)
					{
						atom_found=true;
						atom=get_element_type(e);
						if(atom==-1)
						{
							atom_found=false;
							if(num_atom_type==54)
								atom=21;
							if(num_atom_type==45)
								atom=14;

						}
						polarity=atom_types[atom].polar;
						vdwr=atom_types[atom].vdw+ideal;
					}
					in_max=vdwr+1.0;
					max_voxel=(in_max/vol->units().x())+1.0;
					in_ideal=vdwr+0.6;
					//in_ideal=ideal;
					//in_max=max;

					if(!atom_found)
					{
						fprintf(stderr,"Alert hydropho: Atom type not found. Treated as Hydrogen: %s %d %s\n",name,a->getPdbSerial(),e->name);
					}

					//Exploracion de voxels vecinos
					if(atom!=0 && polarity!=POLAR && (e->symbol==C || e->symbol==F || e->symbol==BR || e->symbol==I || e->symbol==CL))
					{
						//fprintf(stderr,"Elemento: %s\n",(a->getElement())->sym);
						a->getPosition(posA);

						position->x( ((posA[0]-origin.x())/vol->units().x()) );
						position->y( ((posA[1]-origin.y())/vol->units().y()) );
						position->z( ((posA[2]-origin.z())/vol->units().z()) );

						limit_max[0]=(int)(ceil(position->x())+max_voxel);
						limit_max[1]=(int)(ceil(position->y())+max_voxel);
						limit_max[2]=(int)(ceil(position->z())+max_voxel);
						limit_min[0]=(int)(floor(position->x())-max_voxel);
						limit_min[1]=(int)(floor(position->y())-max_voxel);
						limit_min[2]=(int)(floor(position->z())-max_voxel);

						if(limit_max[0]>=vol->dim().x())
							limit_max[0]=vol->dim().x()-1;
						if(limit_max[1]>=vol->dim().y())
							limit_max[1]=vol->dim().y()-1;
						if(limit_max[2]>=vol->dim().z())
							limit_max[2]=vol->dim().z()-1;

						if(limit_min[0]<0)
							limit_min[0]=0;
						if(limit_min[1]<0)
							limit_min[1]=0;
						if(limit_min[2]<0)
							limit_min[2]=0;

						for(ii=limit_min[0];ii<limit_max[0];ii++)
							for(jj=limit_min[1];jj<limit_max[1];jj++)
								for(kk=limit_min[2];kk<limit_max[2];kk++)
								{
									dist=sqrt(((position->x()-ii)*(position->x()-ii))+
											((position->y()-jj)*(position->y()-jj))+
											((position->z()-kk)*(position->z()-kk)));
									dist*=vol->units().x();

//
//									if (dist<vdwr)
//										value==0;
//									if((dist>vdwr) && (dist<in_ideal))
//										value=2.0*(dist-vdwr);
//									if(dist>=in_ideal && dist<in_max)
//										value=1.0*(vdwr-dist+1.5);
//									//value=Bgaussian(dist,in_ideal,in_max,sigma);

									if (dist<vdwr)
										value=0;
									if ((dist>vdwr) && (dist<=in_ideal))
										value=1.0;
									if (dist>=in_max)
										value=0;
									if ((dist>in_ideal) && (dist<in_max))
										value= (1/0.4)*(vdwr+1.0-dist);


									if(dist<in_max && dist>vdwr)
									{
										pos_aux = new vlPoint3ui( ii, jj, kk );
										vol->getVoxel( *pos_aux, value2 );
										//value+=value2;
										//if(value2<value)
											vol->setVoxel( *pos_aux, value );
										delete (pos_aux);
										//cont++;
									}
								}
					}
					//printf("Voxels afectados:%d\n",cont);
				}
			}
		}while (nextAtom());
	//		}





	delete conds;
	delete mol[0];
	delete mol[1];
	delete (position);
	return vol;

}


//Added by Pablo
/// <> pdb proyection using simple gaussian expansion at a given resolution (res)
// the constant a for a Gaussian y=exp(-a*x^2) is
//        a=ln(2)/res^2   for res defined at y=0.5
//        a=1/res^2       for res defined at y=1/e
//
//  here we use the res defined in Fourier space
//  eg. Fourier transform of the above function: Y=exp(-pi^2*k^2/a)
//   a=ln(2)*pi^2/res^2   for res defined at Y=0.5
//   a=pi^2/res^2         for res defined at Y=1/e
//                          reciprocal of the 1/2 width of a Gaussian in Fourier space
//
//   Notes:  the gaussian exp 3/2 factor vanishes due to
//                          sigma_3D=sigma_1D*sqrt(3)
//                          sigma1D=(res*1/2)/width
//
//           res_0.5=1/sqrt(ln2)*res_1/e = 1.20112240878*res_1/e
//

vlVolume * Macromolecule::project_radiusVDW_loop(float unit, int withH,bool center_hyd, float factor_radius)
{
  bool end = true;
  Tcoor atomCoor;
  Element *e;

  Atom * at;
  vlVolume * vlm;
  vlPoint3ui * position;
  float vdwr, value;

  float shiftX, shiftY, shiftZ;
  float centerX, centerY, centerZ, ii, jj, kk;
  float gx, gy, gz;

  int i, j, k;


  /* create paded volume size + kernel size */
  vlm = createVolume_pad( unit, 3 );
	Tcoor volcenter;
  FOPS::center_vol(vlm,volcenter);
  //printf("  grid size %f %f\n",unit, vlm->units().x());
  //printf(" Protein box %d %d %d \n",box[0],box[1],box[2]);
  //printf(" Protein center %f %f %f\n",macromoleculeDimension.geometricCenter[0],
  //       macromoleculeDimension.geometricCenter[1],macromoleculeDimension.geometricCenter[2]);

  /* shift from the PDB geometric center to center of the map */
  shiftX = volcenter[0] - macromoleculeDimension.geometricCenter[0] / vlm->units().x();
  shiftY = volcenter[1] - macromoleculeDimension.geometricCenter[1] / vlm->units().y();
  shiftZ = volcenter[2] - macromoleculeDimension.geometricCenter[2] / vlm->units().z();

  // get geometric center of the map in amstrons
  centerX = ( float ) - unit * volcenter[0] + macromoleculeDimension.geometricCenter[0];
  centerY = ( float ) - unit * volcenter[1] + macromoleculeDimension.geometricCenter[1];
  centerZ = ( float ) - unit * volcenter[2] + macromoleculeDimension.geometricCenter[2];
  vlm->setPosition( vlPoint3f( centerX, centerY, centerZ ) );
  //printf(" Protein center %f %f %f    %f %f %f \n", centerX,centerY,centerX,shiftX, shiftY, shiftZ);


  /* gaussian expansion of pdb */
  initAll();
  position = new vlPoint3ui( 0, 0, 0 );

  value=1;
  while ( end != false )
  {
    /* get atom */
    at = getCurrentAtom();
	e=at->getElement();

    /* getelectronic density */
    if(e->symbol!=H || withH )
        {
    //vdwr=( at->getElement() )->vdw*factor_radius;
    vdwr=1.5*2.0*factor_radius;
    /* get position amstrons */
    at->getPosition(atomCoor);

    /* get relative map position */
    gx = ( float )( atomCoor[0] / unit + shiftX );
    gy = ( float )( atomCoor[1] / unit + shiftY );
    gz = ( float )( atomCoor[2] / unit + shiftZ );

    /* calculate grid pos around point <vdwr*/
    float li;
    li=roundFloat(vdwr/unit)*unit;

	for(ii=-li;ii<=li;ii+=unit)
		for(jj=-li;jj<=li;jj+=unit)
			for(kk=-li;kk<=li;kk+=unit)


			  //if (dist=sqrt(ii*ii+jj*jj+kk*kk)<=vdwr)  	 {
				  if ((ii*ii+jj*jj+kk*kk)<=vdwr*vdwr)  	 {


    i = ( int )roundFloat( gx + ii/unit);
    j = ( int )roundFloat( gy + jj/unit);
    k = ( int )roundFloat( gz + kk/unit);

   // printf("%s en %f %d %d %d  %f %f %f\n",at->getName(),vdwr, i, j, k, ii, jj, kk);

    //value=0.0001+vdwr-dist;
    value=1.0;
    /* set position */
    position->x( i );
    position->y( j );
    position->z( k );
    vlm->setVoxel( * position, value );
			}
        }
    end = nextAtom();
  }

  return vlm;
}



vlVolume * Macromolecule::pdb2map_real( float res, float unit )
{
  bool end = true;
  Tcoor atomCoor;


  Atom * at;
  vlVolume * vlm;
  vlPoint3ui * position;
  float eweight, value;

  float shiftX, shiftY, shiftZ;
  float centerX, centerY, centerZ;

  float g_cte, kg;
  int g_width, x[2], y[2], z[2];
  float gx, gy, gz, dist;
  int i, j, k;

  g_cte = res / unit;
  g_width = ( int )( roundFloat( g_cte * 3.0 ) );

  if ( g_cte / 2.0 * sqrt( 2.0 * log( 2.0 ) ) / sqrt( ( float )3 ) < 1.0 )
  {
    printf( "Increase the sampling and try again\n" );
    exit( 0 );
  }

  g_cte = ( PDB_PI / g_cte ) * ( PDB_PI / g_cte ); // constant for a gaussian cutoff
  // for res defined at y=1/2 Gausian width  g_cte=log(2)*(M_PI/g_cte)*(M_PI/g_cte);
  kg = pow( g_cte / PDB_PI, 1.5 ); // normalization constant

  /* create paded volume size + kernel size */
  vlm = createVolume_pad( unit, g_width );
	Tcoor volcenter;
  FOPS::center_vol(vlm,volcenter);
  //printf(" Gaussian width %d, grid size %f, resolution %f PDB_PI:%f\n",g_width,unit,res,PDB_PI);
  //printf(" Gaussian cutoff %f Norm cte %f\n", g_cte,kg);
  //printf(" Protein box %d %d %d \n",box[0],box[1],box[2]);
  //printf(" Protein center %f %f %f\n",macromoleculeDimension.geometricCenter[0],
  //        macromoleculeDimension.geometricCenter[1],macromoleculeDimension.geometricCenter[2]);

  /* shift from the PDB geometric center to center of the map */
  shiftX = volcenter[0] - macromoleculeDimension.geometricCenter[0] / vlm->units().x();
  shiftY = volcenter[1] - macromoleculeDimension.geometricCenter[1] / vlm->units().y();
  shiftZ = volcenter[2] - macromoleculeDimension.geometricCenter[2] / vlm->units().z();

  // get geometric center of the map in amstrons
  centerX = ( float ) - unit * volcenter[0] + macromoleculeDimension.geometricCenter[0];
  centerY = ( float ) - unit * volcenter[1] + macromoleculeDimension.geometricCenter[1];
  centerZ = ( float ) - unit * volcenter[2] + macromoleculeDimension.geometricCenter[2];
  vlm->setPosition( vlPoint3f( centerX, centerY, centerZ ) );


  /* gaussian expansion of pdb */
  initAll();
  position = new vlPoint3ui( 0, 0, 0 );

  while ( end != false )
  {
    /* get atom */
    at = getCurrentAtom();
    /* getelectronic density */
    eweight = ( at->getElement() )->number;
    /* get position amstrons */
    at->getPosition(atomCoor);
    /* get relative map position */
    gx = ( float )( atomCoor[0] / unit + shiftX );
    gy = ( float )( atomCoor[1] / unit + shiftY );
    gz = ( float )( atomCoor[2] / unit + shiftZ );

    /* calculate kernel gaussian size in grid units */
    x[0] = ( int )roundFloat( gx - g_width ); x[1] = ( int )roundFloat( gx + g_width );
    y[0] = ( int )roundFloat( gy - g_width ); y[1] = ( int )roundFloat( gy + g_width );
    z[0] = ( int )roundFloat( gz - g_width ); z[1] = ( int )roundFloat( gz + g_width );

    /* gaussian expansion */
    for ( k = z[0]; k < z[1]; k++ )
      for ( j = y[0]; j < y[1]; j++ )
        for ( i = x[0]; i < x[1]; i++ )
        {
          /* compute distance ... */
          dist = ( ( i - gx ) * ( i - gx ) + ( j - gy ) * ( j - gy ) + ( k - gz ) * ( k - gz ) );
          /* get position */
          position->x( i );
          position->y( j );
          position->z( k );
          vlm->getVoxel( * position, value );
          /* set position */
          value += kg * eweight * exp( -dist * g_cte );
          vlm->setVoxel( * position, value );
        }

    end = nextAtom();
  }

  return vlm;
}

vlVolume * Macromolecule::pdb2map_real_3BB2R( float res, float unit )
{
  bool end = true;
  Tcoor atomCoor;


  Atom * at;
  vlVolume * vlm;
  vlPoint3ui * position;
  float eweight, value;

  float shiftX, shiftY, shiftZ;
  float centerX, centerY, centerZ;

  float g_cte, kg;
  int g_width, x[2], y[2], z[2];
  float gx, gy, gz, dist;
  int i, j, k;

  g_cte = res / unit;
  g_width = ( int )( roundFloat( g_cte * 3.0 ) );

  if ( g_cte / 2.0 * sqrt( 2.0 * log( 2.0 ) ) / sqrt( ( float )3 ) < 1.0 )
  {
    printf( "Increase the sampling and try again\n" );
    exit( 0 );
  }

  g_cte = ( PDB_PI / g_cte ) * ( PDB_PI / g_cte ); // constant for a gaussian cutoff
  // for res defined at y=1/2 Gausian width  g_cte=log(2)*(M_PI/g_cte)*(M_PI/g_cte);
  kg = pow( g_cte / PDB_PI, 1.5 ); // normalization constant

  /* create paded volume size + kernel size */
  vlm = createVolume_pad( unit, g_width );
	Tcoor volcenter;
  FOPS::center_vol(vlm,volcenter);

  //printf(" Gaussian width %d, grid size %f, resolution %f PDB_PI:%f\n",g_width,unit,res,PDB_PI);
  //printf(" Gaussian cutoff %f Norm cte %f\n", g_cte,kg);
  //printf(" Protein box %d %d %d \n",box[0],box[1],box[2]);
  //printf(" Protein center %f %f %f\n",macromoleculeDimension.geometricCenter[0],
  //        macromoleculeDimension.geometricCenter[1],macromoleculeDimension.geometricCenter[2]);

  /* shift from the PDB geometric center to center of the map */
  shiftX = volcenter[0] - macromoleculeDimension.geometricCenter[0] / vlm->units().x();
  shiftY = volcenter[1] - macromoleculeDimension.geometricCenter[1] / vlm->units().y();
  shiftZ = volcenter[2] - macromoleculeDimension.geometricCenter[2] / vlm->units().z();

  // get geometric center of the map in amstrons
  centerX = ( float ) - unit * volcenter[0] + macromoleculeDimension.geometricCenter[0];
  centerY = ( float ) - unit * volcenter[1] + macromoleculeDimension.geometricCenter[1];
  centerZ = ( float ) - unit * volcenter[2] + macromoleculeDimension.geometricCenter[2];
  vlm->setPosition( vlPoint3f( centerX, centerY, centerZ ) );


  /* gaussian expansion of pdb */
  initAll();
  position = new vlPoint3ui( 0, 0, 0 );

  while ( end != false )
  {
    /* get atom */
    at = getCurrentAtom();
    /* getelectronic density */
//	eweight = at->getPdbocc(); // Mass from the 3BB2R reduced model
	eweight = at->getPdbfact(); // Charge from the 3BB2R reduced model
//    eweight = ( at->getElement() )->number;
    /* get position amstrons */
    at->getPosition(atomCoor);
    /* get relative map position */
    gx = ( float )( atomCoor[0] / unit + shiftX );
    gy = ( float )( atomCoor[1] / unit + shiftY );
    gz = ( float )( atomCoor[2] / unit + shiftZ );

    /* calculate kernel gaussian size in grid units */
    x[0] = ( int )roundFloat( gx - g_width ); x[1] = ( int )roundFloat( gx + g_width );
    y[0] = ( int )roundFloat( gy - g_width ); y[1] = ( int )roundFloat( gy + g_width );
    z[0] = ( int )roundFloat( gz - g_width ); z[1] = ( int )roundFloat( gz + g_width );

    /* gaussian expansion */
    for ( k = z[0]; k < z[1]; k++ )
      for ( j = y[0]; j < y[1]; j++ )
        for ( i = x[0]; i < x[1]; i++ )
        {
          /* compute distance ... */
          dist = ( ( i - gx ) * ( i - gx ) + ( j - gy ) * ( j - gy ) + ( k - gz ) * ( k - gz ) );
          /* get position */
          position->x( i );
          position->y( j );
          position->z( k );
          vlm->getVoxel( * position, value );
          /* set position */
          value += kg * eweight * exp( -dist * g_cte );
          vlm->setVoxel( * position, value );
        }

    end = nextAtom();
  }

  return vlm;
}




//Auxiliar function
/*float elect_limited(float point,float limit, float perc, float in300)
{
	float max=(332.0/(limit*limit*4.0))-in300;
	float init=limit + limit*perc;
	float min=(332.0/(init*init*4.0))-in300;

	float pend=log(max-min)/init;

	return (max-exp(pend*point));
}*/




float elect_limited(float point,float limit, float perc, float in300)
{
	float min=(332.0/(limit*limit*4.0))-in300;
	float init=limit - limit*perc;
	float max=(332.0/(init*init*4.0))-in300;

	float pend=log(max-min)/limit;

	return (max-exp(pend*point));
}



int get_element_type(Element *e)
{
	int atom=-1;

	if(num_atom_type==54)
 	{



 		switch(e->symbol)
		{
			case C:  atom=2; break;
			case H:  atom=21;break;
			case O:  atom=12;break;
			case N:  atom=7; break;
			case S:  atom=15;break;
			case P:  atom=20;break;
			case F:  atom=26;break;
			case CL: atom=27;break;
			case BR: atom=28;break;
			case I:  atom=29;break;
			case ZN: atom=30;break;
			case FE: atom=31;break;
			case MG: atom=33;break;
			case MN: atom=33;break; // MN por MG
			case CA: atom=34;break;
			case NA: atom=35;break;
			case K:  atom=36;break;
			default: atom=-1;
		}
 	}
 	else
 	if(num_atom_type==45)  // sybil type
 	{

  		switch(e->symbol)
		{
			case C: atom=0;break;
			case H: atom=14;break;
			case O: atom=9;break;
			case N: atom=6; break;
			case S: atom=12;break;
			case P: atom=13;break;
			case F:  atom=18;break;
			case CL: atom=19;break;
			case BR: atom=20;break;
			case I:  atom=21;break;
			case ZN: atom=22;break;
			case FE: atom=23;break;
			case MG: atom=24;break;
			case CA: atom=25;break;
			case NA: atom=26;break;
			case K:  atom=27;break;
			default: atom=-1;
		}
	}

	return atom;
}

float original_ele(float dist, float ron, float roff, float radius, float charge,float charge2)
{
	float coeff, coeff1, value;

	if (dist<ron)
	{
		coeff1=1.0;
	}
	else
		if (dist>roff)
		{
			//coeff=0.0;
			return 0.0;
		}

		else
			coeff1=((roff-dist)*(roff-dist)*(roff+2*dist-3*ron))/((roff-ron)*(roff-ron)*(roff-ron));

//	if (dist<5.0)
//		charge2=1.0;
//	else
//		charge=1.0;

	if (dist<4.5)
		coeff=1.0;
	else
		coeff=((roff-dist)*(roff-dist)*(roff+2*dist-3*ron))/((roff-ron)*(roff-ron)*(roff-ron))/5.0;

	if((dist>=radius) && (dist<5.0))
		value = 1.7*coeff*charge*332.0/(4*dist*dist);
	if((dist>=radius) && (dist>=5.0))
      value = coeff1*charge2*332.0/(4*dist*dist);
	if((dist<radius) && (dist<5.0))
		value = 1.7*coeff*charge*332.0/(4*radius*radius);

	return value;
}

float autodock_ele(float dist, float ron, float roff, float radius, float charge, int type)
{
	float coeff, value, coeff1;
	//ron=radius;
	float A=-8.5525;
	float B=78.4+8.5525;
	float k=7.7839;
	float d=0.003627;
	float V1,V2;

    switch(type)
    {
    case 0:
        V1=0.95;   //V1=3.16666/0.30;//3.333*0.95;
        V2=0.25;   //V2=0.8333/0.30; //3.333*0.25;
        break;
    case 1:
        V1=1.8;    //V1=6/0.30;        //3.3333*1.8;
        V2=0.85;   //V2=2.83333/0.30;//3.3333*0.85;
        break;
    case 2:
        V1=2.8;    //V1=9.33333/0.30; //3.333*2.8;
        V2=1.5;    //V2=5.0/0.30;    // 3.333*1.5;
        break;
    case -1:
    default:
        V1=0;
        V2=0;
        break;
    }

		if (dist>=radius)
		{
	    coeff=A+(B/1+k*exp(-d*B*dist));
		}
	    else
	    {
	    coeff=A+(B/1+k*exp(-d*B*radius));
	    }


//	if (dist<radius)
//		coeff=A+(B/1+k*exp(-d*B*radius));
//	else
		if (dist<ron)
		{
			coeff1=1.0;
		}
		else
		if (dist>roff)
		{
			//coeff=0.0;
			return 0.0;
		}
		else
			coeff1=((roff-dist)*(roff-dist)*(roff+2*dist-3*ron))/((roff-ron)*(roff-ron)*(roff-ron));


	if(dist>=radius)
	{
			if(V1==0 && V2==0)
			{
				value=coeff1*charge*(((332.0)/((dist*dist*4.0))));
			}
			else
			{
				if (dist<ron)
					value = V1*coeff1*332.0*charge/(dist*coeff);

				else
					value = V2*coeff1*332.0*charge/(dist*coeff);
			}
	}
	if(dist<radius)
	{
		if(V1==0 && V2==0)
		{
			value=coeff1*charge*(((332.0)/((radius*radius*4.0))));
		}
		else
		{
			value = V1*coeff1*332.0*charge/(radius*coeff);
		}
	}
	return value;
}






#endif
