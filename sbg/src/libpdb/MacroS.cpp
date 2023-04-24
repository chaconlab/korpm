

#include <libpdb/include/Macromolecule.h>
#include <libpdb/include/pdbIter.h>
#include <libtools/include/Surface.h>



void signDistanceGridMol2(S_Grid g, Macromolecule *mol, float PR)

{

	int a,b,c;
	int icx,icy,icz;
	int r;
	float s = 512.0/g->N;
	float dx,dy,dz,fr;
	pdbIter *iter=new pdbIter(mol);
	Atom *at;
	Element *e;
	Tcoor pos;
	float vdw;

	for(iter->pos_atom=0;!iter->gend_atom();iter->next_atom())
	{


		at=iter->get_atom();
		at->getPosition(pos);

		e=at->getElement();
		vdw= (e->vdw+PR)/s;

		r = (int)(vdw)+1;
		icx = int ((pos[0]+256.0)/s);
		icy = int ((pos[1]+256.0)/s);
		icz = int ((pos[2]+256.0)/s);

		//fprintf(stderr,"icx=%d icy=%d icz=%d s=%f\n",icx,icy,icz,s );
		for (a=icx-r;a<=icx+r;a++)
			for (b=icy-r;b<=icy+r;b++)
				for (c=icz-r;c<=icz+r;c++){
					if(a<g->N && b<g->N && c<g->N && a>=0 && b>=0 && c>=0)
					{
						dx = (float)g->matrix[a][b][c].point.x - pos[0];
						dy = (float)g->matrix[a][b][c].point.y - pos[1];
						dz = (float)g->matrix[a][b][c].point.z - pos[2];
						fr = vdw*s;

						//fprintf(stderr,"%d %hd %f %d\n",a,g->matrix[a][b][c].point.x,pos[0]*s,g->N );
						if (dx*dx + dy*dy + dz*dz <= fr*fr)
						{
							g->matrix[a][b][c].phi = 0;
						}
					}
				}
	}
	delete iter;
}






void signDistanceGridMol(S_Grid g, float voxelSize,Macromolecule *mol, float PR)
{
	int a,b,c;
	int icx,icy,icz;
	int r;
        float s = 512.0/g->N;
	float dx,dy,dz,fr;
	pdbIter *iter=new pdbIter(mol);
        Atom *at;
        Element *e;
        Tcoor pos;
        float vdw;

	 for(iter->pos_atom=0;!iter->gend_atom();iter->next_atom())
         {


                at=iter->get_atom();
                at->getPosition(pos);
                pos[0]=pos[0]/voxelSize;
                pos[1]=pos[1]/voxelSize;
                pos[2]=pos[2]/voxelSize;


                e=at->getElement();
                vdw= (e->vdw+PR)/voxelSize;


                r = (int)(vdw)+1;
		icx = (pos[0])+(256.0/s);
		icy = (pos[1])+(256.0/s);
		icz = (pos[2])+(256.0/s);

                //fprintf(stderr,"icx=%d icy=%d icz=%d s=%f\n",icx,icy,icz,s );
                for (a=icx-r;a<=icx+r;a++)
		for (b=icy-r;b<=icy+r;b++)
		for (c=icz-r;c<=icz+r;c++){
                        if(a<g->N && b<g->N && c<g->N && a>=0 && b>=0 && c>=0)
                        {
                          dx = (float)g->matrix[a][b][c].point.x - pos[0]*s;
			  dy = (float)g->matrix[a][b][c].point.y - pos[1]*s;
			  dz = (float)g->matrix[a][b][c].point.z - pos[2]*s;
			  fr = vdw*s;

                          //fprintf(stderr,"%d %hd %f %d\n",a,g->matrix[a][b][c].point.x,pos[0]*s,g->N );
			  if (dx*dx + dy*dy + dz*dz <= fr*fr)
			  {
                            g->matrix[a][b][c].phi = 0;
			  }
                        }
		}
	}
        delete iter;
}


void Macromolecule::lsms( float PR,  bool inner )
{
	Macromolecule *pdb2;
	Tcoor original_center;


  int maxLength;


	pdb2=new Macromolecule(this);
	pdb2->geoCenter( original_center );
    pdb2->geoBox();

  original_center[0] = pdb2->macromoleculeDimension.geometricCenter[0] * ( -1.0 );
  original_center[1] = pdb2->macromoleculeDimension.geometricCenter[1] * ( -1.0 );
  original_center[2] = pdb2->macromoleculeDimension.geometricCenter[2] * ( -1.0 );
  pdb2->moveAll( original_center );


  maxLength=floor (pdb2->maxLength()+2*PR+3);
//  fprintf(stderr,"createGrid maxLength=%f %d\n",pdb2->maxLength(), maxLength );
  S_Grid grid=  createGrid(512);
//  fprintf(stderr,"signDistanceGridMol\n");
  signDistanceGridMol2(grid, pdb2, PR);
//  fprintf(stderr,"shrink\n");
  shrink(grid,PR);
//  fprintf(stderr,"fastMarching\n");
//  fastMarching(grid, inner);
//  int cont=0;
//  for(m=0;m<maxLength;m++)
//      for (n=0;n<maxLength;n++)
//        for (k=0;k<maxLength;k++)
//        {
//          if(grid->matrix[m][n][k].phi==1) cont++;
//        }
//  fprintf(stderr,"Marching Cube%d\n", cont);
  marchingCube(grid);

  fprintf(stdout,"lsms> Saved surf.vmd\n");
  destroyGrid(grid);
  delete pdb2;

}



#ifdef VOLUME_INCLUDED

vlVolume * Macromolecule::projectSurface( float unit, float PR, bool inner )
{
  vlVolume *vol;
  int maxLength= (this->maxLength()+(PR+3)*2)/unit;
  int m,n,k;
  vlDim pad;

  //fprintf(stderr,"maxLength=%f\n",maxLength);
  vol=createVolume( unit );
  pad.x((maxLength-vol->dim().x())/2);
  pad.y((maxLength-vol->dim().y())/2);
  pad.z((maxLength-vol->dim().z())/2);
  vol=FOPS::padVolume(vol, pad,true);

  if(vol->dim().x()!=maxLength || vol->dim().y()!=maxLength || vol->dim().z()!=maxLength)
  {
    pad.x(1);pad.y(1);pad.z(1);
    //fprintf(stderr,"Dimensiones erroneas: %d %d %d maxLength=%d\n"
    //,vol->dim().x(),vol->dim().y(),vol->dim().z(),maxLength);

    vol=FOPS::padVolume(vol, pad,1,true);

    //fprintf(stderr,"Dimensiones erroneas: %d %d %d maxLength=%d\n"
    //,vol->dim().x(),vol->dim().y(),vol->dim().z(),maxLength);
  }


  //fprintf(stderr,"createGrid maxLength=%d\n",maxLength);
  S_Grid grid=  createGrid(maxLength);
  //fprintf(stderr,"signDistanceGridMol\n");
  signDistanceGridMol(grid, unit ,this, PR);
  //fprintf(stderr,"shrink\n");
  shrink(grid,PR);
  //fprintf(stderr,"fastMarching\n");
  fastMarching(grid, inner);


  for(m=0;m<maxLength;m++)
    for (n=0;n<maxLength;n++)
      for (k=0;k<maxLength;k++)
      {
        if(grid->matrix[m][n][k].phi==0)
        {
          vol->setVoxel(vlPoint3ui(m,n,k),(float)1.0);
        }
      }

  destroyGrid(grid);
  return(vol);
}

void Macromolecule::projectSurface( vlVolume *vol, float PR, bool inner )
{
	Macromolecule *pdb2;
	Tcoor original_center;
	float unit=vol->units().x();
  float diff[3];
  //int maxLength= (this->maxLength()+(PR+3)*2)/unit;
  int maxLength;
  int m,n,k;
  vlDim pad;

	pdb2=new Macromolecule(this);
    pdb2->geoCenter( original_center );

//     pdb2->geoBox();
//     original_center[0] = pdb2->macromoleculeDimension.geometricCenter[0] * ( -1.0 );
//     original_center[1] = pdb2->macromoleculeDimension.geometricCenter[1] * ( -1.0 );
//     original_center[2] = pdb2->macromoleculeDimension.geometricCenter[2] * ( -1.0 );
     original_center[0]*=-1.0;
     original_center[1]*=-1.0;
     original_center[2]*=-1.0;
     pdb2->moveAll( original_center );



	maxLength=vol->dim().x();
	if(maxLength<vol->dim().y())
		maxLength=vol->dim().y();
	if(maxLength<vol->dim().z())
		maxLength=vol->dim().z();




 	vol->clear(0.0);

	diff[0]=maxLength-vol->dim().x();
	diff[1]=maxLength-vol->dim().y();
	diff[2]=maxLength-vol->dim().z();

	for(m=0;m<3;m++)
	{
		diff[m]=diff[m]/2.0;
	}


  //fprintf(stderr,"createGrid maxLength=%d\n",maxLength);
  S_Grid grid=  createGrid(maxLength);
  //fprintf(stderr,"signDistanceGridMol\n");
  signDistanceGridMol(grid, unit ,pdb2, PR);
  //fprintf(stderr,"shrink\n");
  shrink(grid,PR);
  //fprintf(stderr,"fastMarching\n");
  fastMarching(grid, inner);

  for(m=0;m<maxLength;m++)
    for (n=0;n<maxLength;n++)
      for (k=0;k<maxLength;k++)
      {
        if(grid->matrix[m][n][k].phi==0)
        {

          vol->interpolateVoxel(vlPoint3f((float)m-diff[0],(float)n-diff[1],(float)k-diff[2]),(float)1.0);
        }
      }

  destroyGrid(grid);
  delete pdb2;

}


#endif
