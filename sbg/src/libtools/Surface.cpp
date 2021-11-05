
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "Surface.h"


S_Grid createGrid(short i)
{
	float s;
	short m,n,k;
	S_Grid g = (S_Grid)malloc(sizeof(struct g));
	g->N = i;

	s = 512.0/(float)i;
	g->stepSize = s;
	g->matrix = (S_GridPoint***)malloc(i*sizeof(S_GridPoint**));
	for (m=0;m<i;m++)
		g->matrix[m] = (S_GridPoint**)malloc(i*sizeof(S_GridPoint*));
	for (m=0;m<i;m++)
		for (n=0;n<i;n++)
			g->matrix[m][n] = (S_GridPoint*)malloc(i*sizeof(S_GridPoint));

	for(m=0;m<i;m++)
		for (n=0;n<i;n++)
			for (k=0;k<i;k++)
	{
		g->matrix[m][n][k].point.x = -256+(short)((float)m*s);
        g->matrix[m][n][k].point.y = -256+(short)((float)n*s);
		g->matrix[m][n][k].point.z = -256+(short)((float)k*s);
		g->matrix[m][n][k].phi = 1; // everything is outside surface
		g->matrix[m][n][k].from=-2;
		g->matrix[m][n][k].dist= 0;
	}


  return g;
}


S_Grid copyGrid(S_Grid g2)
{
        float s;
        short m,n,k;
        S_Grid g = (S_Grid)malloc(sizeof(struct g));
        g->N = g2->N;
        s = 512.0/(float)g2->N;
        g->stepSize = s;
        g->matrix = (S_GridPoint***)malloc(g2->N*sizeof(S_GridPoint**));
        for (m=0;m<g2->N;m++)
                g->matrix[m] = (S_GridPoint**)malloc(g2->N*sizeof(S_GridPoint*));
        for (m=0;m<g2->N;m++)
                for (n=0;n<g2->N;n++)
                        g->matrix[m][n] = (S_GridPoint*)malloc(g2->N*sizeof(S_GridPoint));

        for(m=0;m<g2->N;m++)
                for (n=0;n<g2->N;n++)
                        for (k=0;k<g2->N;k++)
        {
                g->matrix[m][n][k].point.x = g2->matrix[m][n][k].point.x;
                g->matrix[m][n][k].point.y = g2->matrix[m][n][k].point.y;
                g->matrix[m][n][k].point.z = g2->matrix[m][n][k].point.z;
                g->matrix[m][n][k].phi = g2->matrix[m][n][k].phi;
                g->matrix[m][n][k].from= g2->matrix[m][n][k].from;
                g->matrix[m][n][k].dist= g2->matrix[m][n][k].dist;
        }
        return g;
}

void destroyGrid(S_Grid g)
{
        int m,n;
        short i=g->N ;
        for (m=0;m<i;m++)
        {
          for (n=0;n<i;n++)
               free(g->matrix[m][n]);
          free(g->matrix[m]);
        }
        free (g->matrix);
        free (g);
}


int findProbesMol(S_Grid g, float PR)
{
	int i,j,k;
	int a,b,c;
	int co = 0;
	float fx,fy,fz,fr;
	float dx,dy,dz;
	int icx,icy,icz,r;
	int l = g->N;
	float s = 512/g->N;
        float PR2=PR*s;

	co = 0;

	for (i=1;i<l-1;i++)
	for (j=1;j<l-1;j++)
	for (k=1;k<l-1;k++)
		if (g->matrix[i][j][k].phi==1 &&
		   (g->matrix[i][j][k-1].phi==0 ||
		    g->matrix[i][j][k+1].phi==0 ||
		    g->matrix[i][j-1][k].phi==0 ||
		    g->matrix[i][j+1][k].phi==0 ||
		    g->matrix[i+1][j][k].phi==0 ||
		    g->matrix[i-1][j][k].phi==0))
			co++;

	S_Molecule m = (S_Molecule)malloc(co*sizeof(struct S_prot));
	m->npoints = co;
	m->xpoints = (float *)malloc(co*sizeof(float));
	m->ypoints = (float *)malloc(co*sizeof(float));
	m->zpoints = (float *)malloc(co*sizeof(float));
	m->rpoints = (float *)malloc(co*sizeof(float));
	co = 0;
	for (i=1;i<l-1;i++)
	for (j=1;j<l-1;j++)
	for (k=1;k<l-1;k++)
		if (g->matrix[i][j][k].phi==1 &&
		   (g->matrix[i][j][k-1].phi==0 ||
		    g->matrix[i][j][k+1].phi==0 ||
		    g->matrix[i][j-1][k].phi==0 ||
		    g->matrix[i][j+1][k].phi==0 ||
		    g->matrix[i+1][j][k].phi==0 ||
		    g->matrix[i-1][j][k].phi==0)){
			m->xpoints[co] = g->matrix[i][j][k].point.x;
			m->ypoints[co] = g->matrix[i][j][k].point.y;
			m->zpoints[co] = g->matrix[i][j][k].point.z;
			m->rpoints[co] = PR2;
			co++;
		    }

	for (i=0;i<m->npoints;i++){
		fx = m->xpoints[i];
		fy = m->ypoints[i];
		fz = m->zpoints[i];
		fr = m->rpoints[i];
		r = (int)(fr/s)+1;
		icx = (fx+256)/s;
		icy = (fy+256)/s;
		icz = (fz+256)/s;
		for (a=icx-r;a<=icx+r;a++)
		for (b=icy-r;b<=icy+r;b++)
		for (c=icz-r;c<=icz+r;c++){
			if (a>=g->N || b>=g->N || c>=g->N || a<0 || b<0 || c<0)
				continue;
			dx = g->matrix[a][b][c].point.x - fx;
			dy = g->matrix[a][b][c].point.y - fy;
			dz = g->matrix[a][b][c].point.z - fz;
			if (dx*dx + dy*dy + dz*dz <= fr*fr)
				g->matrix[a][b][c].phi = 1;
		}
	}

	return co;
}

void shrink(S_Grid g, float PR)
{
        int xp,xn,yp,yn,zp,zn;
        int i,j,k;
        int co = 0,p;
        float fr;
        float dx,dy,dz;
        int l = g->N;
	float s = 512/g->N;
        float temp_dist;
        S_PPoint* nb_head, *nb_orig, *nb_next_head,*nb_next_orig;
         float PR2=PR*s;
//fprintf(stderr,"FLAG1\n");
        PR2=PR2*PR2;	//be careful
        co = 0;

        for (i=1;i<l-1;i++)
        for (j=1;j<l-1;j++)
        for (k=1;k<l-1;k++)
                if (g->matrix[i][j][k].phi==1 &&
                   (g->matrix[i][j][k-1].phi==0 ||
                    g->matrix[i][j][k+1].phi==0 ||
                    g->matrix[i][j-1][k].phi==0 ||
                    g->matrix[i][j+1][k].phi==0 ||
                    g->matrix[i+1][j][k].phi==0 ||
                    g->matrix[i-1][j][k].phi==0))
                        co++;

        nb_head = (S_PPoint*)malloc(co*sizeof(S_PPoint));	//head line
        nb_orig = (S_PPoint*)malloc(co*sizeof(S_PPoint));	//where does this point come from
//fprintf(stderr,"FLAG2 co=%d l=%d\n",co,l);

        co=0;
        for (i=1;i<l-1;i++)
        for (j=1;j<l-1;j++)
        for (k=1;k<l-1;k++){
                //should we start from first "1" or from last "0"
                if (g->matrix[i][j][k].phi==1){
                        if(g->matrix[i][j][k-1].phi==0 ||
                                    g->matrix[i][j][k+1].phi==0 ||
                                    g->matrix[i][j-1][k].phi==0 ||
                                    g->matrix[i][j+1][k].phi==0 ||
                                    g->matrix[i+1][j][k].phi==0 ||
                                   g->matrix[i-1][j][k].phi==0){
                                        nb_head[co].x = i;
                                        nb_head[co].y = j;
                                        nb_head[co].z = k;
                                        nb_orig[co].x = i;
                                        nb_orig[co].y = j;
                                        nb_orig[co].z = k;
                                        co++;
                        }
                }
        }

//fprintf(stderr,"FLAG2 co=%d\n",co);


        while (co!=0){
//fprintf(stderr,"%d\n",co);
                nb_next_head = (S_PPoint*)malloc(co*6*sizeof(S_PPoint));	//next level (at most 6 times larger)
                nb_next_orig = (S_PPoint*)malloc(co*6*sizeof(S_PPoint));	//it is possible to be larger than nb_head
//fprintf(stderr,"FLAG3\n");
                p = 0;
                for (i=0;i<co;i++){
                        xp=nb_head[i].x+1;
                        xn=nb_head[i].x-1;
                        yp=nb_head[i].y+1;
                        yn=nb_head[i].y-1;
                        zp=nb_head[i].z+1;
                        zn=nb_head[i].z-1;
//fprintf(stderr,"FLAG3.5 xp=%d nheady=%d nheadz=%d  i=%d co=%d l=%d\n",xp,nb_head[i].y,nb_head[i].z,i,co,l);
                        if(xp<l)
                        if (g->matrix[xp][nb_head[i].y][nb_head[i].z].phi==0 ||
                           (g->matrix[xp][nb_head[i].y][nb_head[i].z].phi==1&&g->matrix[xp][nb_head[i].y][nb_head[i].z].from!=-1)){
//fprintf(stderr,"FLAG3.55\n");
                                if(g->matrix[xp][nb_head[i].y][nb_head[i].z].phi==0)	fr=PR2;
                                else	fr=g->matrix[xp][nb_head[i].y][nb_head[i].z].dist;
//fprintf(stderr,"FLAG3.6\n");
                                dx = g->matrix[xp][nb_head[i].y][nb_head[i].z].point.x - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.x;
                                dy = g->matrix[xp][nb_head[i].y][nb_head[i].z].point.y - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.y;
                                dz = g->matrix[xp][nb_head[i].y][nb_head[i].z].point.z - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.z;
                                temp_dist=dx*dx+dy*dy+dz*dz;
//fprintf(stderr,"FLAG3.7\n");
                                if(temp_dist<=fr){
                                        nb_next_head[p].x=xp;
                                        nb_next_head[p].y=nb_head[i].y;
                                        nb_next_head[p].z=nb_head[i].z;
                                        nb_next_orig[p].x=nb_orig[i].x;
                                        nb_next_orig[p].y=nb_orig[i].y;
                                        nb_next_orig[p].z=nb_orig[i].z;
                                        g->matrix[xp][nb_head[i].y][nb_head[i].z].from=i;//index of nb_head
                                        g->matrix[xp][nb_head[i].y][nb_head[i].z].dist=temp_dist;
                                        p++;
//fprintf(stderr,"FLAG3.8\n");
                                        if(g->matrix[xp][nb_head[i].y][nb_head[i].z].phi==1 && temp_dist==fr)
                                                p--;	//doesn't count
                                        g->matrix[xp][nb_head[i].y][nb_head[i].z].phi=1;
                                }
                        }
//fprintf(stderr,"FLAG4\n");

                        if(xn>-1)
                        if (g->matrix[xn][nb_head[i].y][nb_head[i].z].phi==0 ||
                           (g->matrix[xn][nb_head[i].y][nb_head[i].z].phi==1&&g->matrix[xn][nb_head[i].y][nb_head[i].z].from!=-1)){
                                if(g->matrix[xn][nb_head[i].y][nb_head[i].z].phi==0)	fr=PR2;
                                else	fr=g->matrix[xn][nb_head[i].y][nb_head[i].z].dist;

                                dx = g->matrix[xn][nb_head[i].y][nb_head[i].z].point.x - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.x;
                                dy = g->matrix[xn][nb_head[i].y][nb_head[i].z].point.y - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.y;
                                dz = g->matrix[xn][nb_head[i].y][nb_head[i].z].point.z - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.z;
                                temp_dist=dx*dx+dy*dy+dz*dz;

                                if(temp_dist<=fr){
                                        nb_next_head[p].x=xn;
                                        nb_next_head[p].y=nb_head[i].y;
                                        nb_next_head[p].z=nb_head[i].z;
                                        nb_next_orig[p].x=nb_orig[i].x;
                                        nb_next_orig[p].y=nb_orig[i].y;
                                        nb_next_orig[p].z=nb_orig[i].z;
                                        g->matrix[xn][nb_head[i].y][nb_head[i].z].from=i;//index of nb_head
                                        g->matrix[xn][nb_head[i].y][nb_head[i].z].dist=temp_dist;
                                        p++;

                                        if(g->matrix[xn][nb_head[i].y][nb_head[i].z].phi==1 && temp_dist==fr)
                                                p--;
                                        g->matrix[xn][nb_head[i].y][nb_head[i].z].phi=1;
                                }
                        }
//fprintf(stderr,"FLAG5\n");

                        if(yp<l)
                        if (g->matrix[nb_head[i].x][yp][nb_head[i].z].phi==0 ||
                           (g->matrix[nb_head[i].x][yp][nb_head[i].z].phi==1&&g->matrix[nb_head[i].x][yp][nb_head[i].z].from!=-1)){
                                if(g->matrix[nb_head[i].x][yp][nb_head[i].z].phi==0)	fr=PR2;
                                else	fr=g->matrix[nb_head[i].x][yp][nb_head[i].z].dist;

                                dx = g->matrix[nb_head[i].x][yp][nb_head[i].z].point.x - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.x;
                                dy = g->matrix[nb_head[i].x][yp][nb_head[i].z].point.y - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.y;
                                dz = g->matrix[nb_head[i].x][yp][nb_head[i].z].point.z - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.z;
                                temp_dist=dx*dx+dy*dy+dz*dz;

                                if(temp_dist<=fr){
                                        nb_next_head[p].x=nb_head[i].x;
                                        nb_next_head[p].y=yp;
                                        nb_next_head[p].z=nb_head[i].z;
                                        nb_next_orig[p].x=nb_orig[i].x;
                                        nb_next_orig[p].y=nb_orig[i].y;
                                        nb_next_orig[p].z=nb_orig[i].z;
                                        g->matrix[nb_head[i].x][yp][nb_head[i].z].from=i;//index of nb_head
                                        g->matrix[nb_head[i].x][yp][nb_head[i].z].dist=temp_dist;
                                        p++;
                                        if(g->matrix[nb_head[i].x][yp][nb_head[i].z].phi==1 && temp_dist==fr)
                                                p--;
                                        g->matrix[nb_head[i].x][yp][nb_head[i].z].phi=1;
                                }
                        }
//fprintf(stderr,"FLAG6\n");

												if(yn>-1)
                        if (g->matrix[nb_head[i].x][yn][nb_head[i].z].phi==0 ||
                           (g->matrix[nb_head[i].x][yn][nb_head[i].z].phi==1&&g->matrix[nb_head[i].x][yn][nb_head[i].z].from!=-1)){
                                if(g->matrix[nb_head[i].x][yn][nb_head[i].z].phi==0)	fr=PR2;
                                else	fr=g->matrix[nb_head[i].x][yn][nb_head[i].z].dist;

                                dx = g->matrix[nb_head[i].x][yn][nb_head[i].z].point.x - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.x;
                                dy = g->matrix[nb_head[i].x][yn][nb_head[i].z].point.y - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.y;
                                dz = g->matrix[nb_head[i].x][yn][nb_head[i].z].point.z - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.z;
                                temp_dist=dx*dx+dy*dy+dz*dz;

                                if(temp_dist<=fr){
                                        nb_next_head[p].x=nb_head[i].x;
                                        nb_next_head[p].y=yn;
                                        nb_next_head[p].z=nb_head[i].z;
                                        nb_next_orig[p].x=nb_orig[i].x;
                                        nb_next_orig[p].y=nb_orig[i].y;
                                        nb_next_orig[p].z=nb_orig[i].z;
                                        g->matrix[nb_head[i].x][yn][nb_head[i].z].from=i;//index of nb_head
                                        g->matrix[nb_head[i].x][yn][nb_head[i].z].dist=temp_dist;
                                        p++;
                                        if(g->matrix[nb_head[i].x][yn][nb_head[i].z].phi==1 && temp_dist==fr)
                                                p--;
                                        g->matrix[nb_head[i].x][yn][nb_head[i].z].phi=1;
                                }
                        }
//fprintf(stderr,"FLAG7\n");

                        if(zp<l)
                        if (g->matrix[nb_head[i].x][nb_head[i].y][zp].phi==0 ||
                           (g->matrix[nb_head[i].x][nb_head[i].y][zp].phi==1&&g->matrix[nb_head[i].x][nb_head[i].y][zp].from!=-1)){
                                if(g->matrix[nb_head[i].x][nb_head[i].y][zp].phi==0)	fr=PR2;
                                else	fr=g->matrix[nb_head[i].x][nb_head[i].y][zp].dist;

                                dx = g->matrix[nb_head[i].x][nb_head[i].y][zp].point.x - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.x;
                                dy = g->matrix[nb_head[i].x][nb_head[i].y][zp].point.y - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.y;
                                dz = g->matrix[nb_head[i].x][nb_head[i].y][zp].point.z - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.z;
                                temp_dist=dx*dx+dy*dy+dz*dz;

                                if(temp_dist<=fr){
                                        nb_next_head[p].x=nb_head[i].x;
                                        nb_next_head[p].y=nb_head[i].y;
                                        nb_next_head[p].z=zp;
                                        nb_next_orig[p].x=nb_orig[i].x;
                                        nb_next_orig[p].y=nb_orig[i].y;
                                        nb_next_orig[p].z=nb_orig[i].z;
                                        g->matrix[nb_head[i].x][nb_head[i].y][zp].from=i;//index of nb_head
                                        g->matrix[nb_head[i].x][nb_head[i].y][zp].dist=temp_dist;
                                        p++;
                                        if(g->matrix[nb_head[i].x][nb_head[i].y][zp].phi==1 && temp_dist==fr)
                                                p--;
                                        g->matrix[nb_head[i].x][nb_head[i].y][zp].phi=1;
                                }
                        }
//fprintf(stderr,"FLAG8\n");

                        if(zn>-1)
                        if (g->matrix[nb_head[i].x][nb_head[i].y][zn].phi==0 ||
                           (g->matrix[nb_head[i].x][nb_head[i].y][zn].phi==1&&g->matrix[nb_head[i].x][nb_head[i].y][zn].from!=-1)){
//fprintf(stderr,"FLAG8_1\n");
                                if(g->matrix[nb_head[i].x][nb_head[i].y][zn].phi==0)	fr=PR2;
                                else	fr=g->matrix[nb_head[i].x][nb_head[i].y][zn].dist;
//fprintf(stderr,"FLAG8_2\n");
                                dx = g->matrix[nb_head[i].x][nb_head[i].y][zn].point.x - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.x;
                                dy = g->matrix[nb_head[i].x][nb_head[i].y][zn].point.y - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.y;
                                dz = g->matrix[nb_head[i].x][nb_head[i].y][zn].point.z - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.z;
                                temp_dist=dx*dx+dy*dy+dz*dz;
//fprintf(stderr,"FLAG8_3\n");
                                if(temp_dist<=fr){
                                        nb_next_head[p].x=nb_head[i].x;
                                        nb_next_head[p].y=nb_head[i].y;
                                        nb_next_head[p].z=zn;
                                        nb_next_orig[p].x=nb_orig[i].x;
                                        nb_next_orig[p].y=nb_orig[i].y;
                                        nb_next_orig[p].z=nb_orig[i].z;
                                        g->matrix[nb_head[i].x][nb_head[i].y][zn].from=i;//index of nb_head
                                        g->matrix[nb_head[i].x][nb_head[i].y][zn].dist=temp_dist;
                                        p++;
                                        if(g->matrix[nb_head[i].x][nb_head[i].y][zn].phi==1 && temp_dist==fr)
                                                p--;
                                        g->matrix[nb_head[i].x][nb_head[i].y][zn].phi=1;
                                }
//fprintf(stderr,"FLAG8_4\n");
                        }
                }
//fprintf(stderr,"FLAG9\n");

                free(nb_head);
                free(nb_orig);
                nb_head = nb_next_head;
                nb_orig = nb_next_orig;
                co=p;
//fprintf(stderr,"FLAG10\n");

        }
//fprintf(stderr,"FLAG11\n");

        free(nb_head);
        free(nb_orig);
//fprintf(stderr,"FLAG12\n");

}

int fastMarching(S_Grid g, bool inner)
{
        int len = g->N;
        int i;
        int p;
        int x,y,z;
        S_PPoint* nb_head;
        S_PPoint* tmp;
        int co;
        int vol = 0;

        i=0;
        p=len-1;
        //initialize outer surface and narrow band
        for (x=i;x<=p;x++)
                for (y=i;y<=p;y++){
                        g->matrix[x][y][i].phi = 3;
                        g->matrix[x][y][p].phi = 3;
                        g->matrix[x][i][y].phi = 3;
                        g->matrix[x][p][y].phi = 3;
                        g->matrix[i][x][y].phi = 3;
                        g->matrix[p][x][y].phi = 3;
                }
        i++; p--;
//because we use + 1, we directly set a bend of boundary
        co = 0;
        for (x=i;x<=p;x++)
                for (y=i;y<=p;y++){
                        g->matrix[x][y][i].phi = 3;
                        g->matrix[x][y][p].phi = 3;
                        g->matrix[x][i][y].phi = 3;
                        g->matrix[x][p][y].phi = 3;
                        g->matrix[i][x][y].phi = 3;
                        g->matrix[p][x][y].phi = 3;
                        co+=6;
                }

        nb_head = (S_PPoint*)malloc(co*sizeof(S_PPoint));

        co = 0;
        for (x=i;x<=p;x++)
                for (y=i;y<=p;y++){
                        nb_head[co].x = x;
                        nb_head[co].y = y;
                        nb_head[co++].z = i;

                        nb_head[co].x = x;
                        nb_head[co].y = y;
                        nb_head[co++].z = p;

                        nb_head[co].x = x;
                        nb_head[co].y = i;
                        nb_head[co++].z = y;

                        nb_head[co].x = x;
                        nb_head[co].y = p;
                        nb_head[co++].z = y;

                        nb_head[co].x = i;
                        nb_head[co].y = x;
                        nb_head[co++].z = y;

                        nb_head[co].x = p;
                        nb_head[co].y = x;
                        nb_head[co++].z = y;
                }

        while (co!=0){
                p = 0;
                tmp = (S_PPoint*)malloc(co*6*sizeof(S_PPoint));	//at most 6 times larger than nb_head
                for (i=0;i<co;i++){
                        if (g->matrix[nb_head[i].x+1][nb_head[i].y][nb_head[i].z].phi==1){
                                g->matrix[nb_head[i].x+1][nb_head[i].y][nb_head[i].z].phi=3;
                                tmp[p].x=nb_head[i].x+1;
                                tmp[p].y=nb_head[i].y;
                                tmp[p++].z=nb_head[i].z;
                        }
                        if (g->matrix[nb_head[i].x-1][nb_head[i].y][nb_head[i].z].phi==1){
                                g->matrix[nb_head[i].x-1][nb_head[i].y][nb_head[i].z].phi=3;
                                tmp[p].x=nb_head[i].x-1;
                                tmp[p].y=nb_head[i].y;
                                tmp[p++].z=nb_head[i].z;
                        }
                        if (g->matrix[nb_head[i].x][nb_head[i].y+1][nb_head[i].z].phi==1){
                                g->matrix[nb_head[i].x][nb_head[i].y+1][nb_head[i].z].phi=3;
                                tmp[p].x=nb_head[i].x;
                                tmp[p].y=nb_head[i].y+1;
                                tmp[p++].z=nb_head[i].z;
                        }
                        if (g->matrix[nb_head[i].x][nb_head[i].y-1][nb_head[i].z].phi==1){
                                g->matrix[nb_head[i].x][nb_head[i].y-1][nb_head[i].z].phi=3;
                                tmp[p].x=nb_head[i].x;
                                tmp[p].y=nb_head[i].y-1;
                                tmp[p++].z=nb_head[i].z;
                        }
                        if (g->matrix[nb_head[i].x][nb_head[i].y][nb_head[i].z+1].phi==1){
                                g->matrix[nb_head[i].x][nb_head[i].y][nb_head[i].z+1].phi=3;
                                tmp[p].x=nb_head[i].x;
                                tmp[p].y=nb_head[i].y;
                                tmp[p++].z=nb_head[i].z+1;
                        }
                        if (g->matrix[nb_head[i].x][nb_head[i].y][nb_head[i].z-1].phi==1){
                                g->matrix[nb_head[i].x][nb_head[i].y][nb_head[i].z-1].phi=3;
                                tmp[p].x=nb_head[i].x;
                                tmp[p].y=nb_head[i].y;
                                tmp[p++].z=nb_head[i].z-1;
                        }
                }
                free(nb_head);
                nb_head = tmp;
                co = p;
        }
        free(nb_head);

        // outer surface
        if (!inner){
                for (x=0;x<len;x++)
                  for (y=0;y<len;y++)
                    for (z=0;z<len;z++)
                        if (g->matrix[x][y][z].phi==3)
                                g->matrix[x][y][z].phi=1;
                        else{
                                g->matrix[x][y][z].phi=0;
                                vol++;
                        }
        }else{
        // inner cave
                for (x=0;x<len;x++)
                  for (y=0;y<len;y++)
                    for (z=0;z<len;z++)
                        if (g->matrix[x][y][z].phi!=1)
                                g->matrix[x][y][z].phi=1;
                        else{
                                g->matrix[x][y][z].phi=0;
                                vol++;
                        }
        }
        return vol;
}



void expand(S_Grid g, float PR)
{
        int xp,xn,yp,yn,zp,zn;
        int i,j,k;
        int co = 0,p;
        float fr;
        float dx,dy,dz;
        int l = g->N;
	    float s = 512.0/g->N;
        float temp_dist;
        S_PPoint* nb_head, *nb_orig, *nb_next_head,*nb_next_orig;
        float PR2=PR*s;
				//fprintf(stdout,"PR=%f PR2=%f N=%d s=%f\n",PR,PR2,l,s);


        PR2=PR2*PR2;	//be careful
        co = 0;

        for (i=1;i<l-1;i++)
        for (j=1;j<l-1;j++)
        for (k=1;k<l-1;k++)
                if (g->matrix[i][j][k].phi==0 &&
                   (g->matrix[i][j][k-1].phi==1 ||
                    g->matrix[i][j][k+1].phi==1 ||
                    g->matrix[i][j-1][k].phi==1 ||
                    g->matrix[i][j+1][k].phi==1 ||
                    g->matrix[i+1][j][k].phi==1 ||
                    g->matrix[i-1][j][k].phi==1))
                        co++;

        nb_head = (S_PPoint*)malloc(co*sizeof(S_PPoint));	//head line
        nb_orig = (S_PPoint*)malloc(co*sizeof(S_PPoint));	//where does this point come from

        co=0;
        for (i=1;i<l-1;i++)
        for (j=1;j<l-1;j++)
        for (k=1;k<l-1;k++){
                //should we start from first "1" or from last "0"
                if (g->matrix[i][j][k].phi==0){
                        if(g->matrix[i][j][k-1].phi==1 ||
                          g->matrix[i][j][k+1].phi==1 ||
                          g->matrix[i][j-1][k].phi==1 ||
                          g->matrix[i][j+1][k].phi==1 ||
                          g->matrix[i+1][j][k].phi==1 ||
                          g->matrix[i-1][j][k].phi==1){
                            nb_head[co].x = i;
                            nb_head[co].y = j;
                            nb_head[co].z = k;
                            nb_orig[co].x = i;
                            nb_orig[co].y = j;
                            nb_orig[co].z = k;
                            co++;
                        }
                }
        }

        while (co!=0){
//printf("%d\n",co);
                nb_next_head = (S_PPoint*)malloc(co*6*sizeof(S_PPoint));	//next level (at most 6 times larger)
                nb_next_orig = (S_PPoint*)malloc(co*6*sizeof(S_PPoint));	//it is possible to be larger than nb_head

                p = 0;
                for (i=0;i<co;i++){
                        xp=nb_head[i].x+1;
                        xn=nb_head[i].x-1;
                        yp=nb_head[i].y+1;
                        yn=nb_head[i].y-1;
                        zp=nb_head[i].z+1;
                        zn=nb_head[i].z-1;
                        if (g->matrix[xp][nb_head[i].y][nb_head[i].z].phi==1 ||
                           (g->matrix[xp][nb_head[i].y][nb_head[i].z].phi==0&&g->matrix[xp][nb_head[i].y][nb_head[i].z].from!=-1)){
                                if(g->matrix[xp][nb_head[i].y][nb_head[i].z].phi==1)	fr=PR2;
                                else	fr=g->matrix[xp][nb_head[i].y][nb_head[i].z].dist;

                                dx = g->matrix[xp][nb_head[i].y][nb_head[i].z].point.x - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.x;
                                dy = g->matrix[xp][nb_head[i].y][nb_head[i].z].point.y - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.y;
                                dz = g->matrix[xp][nb_head[i].y][nb_head[i].z].point.z - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.z;
                                temp_dist=dx*dx+dy*dy+dz*dz;

                                if((fr-temp_dist)>=0.01){
                                //fprintf(stderr,"1 temp_dist=%f fr=%f dx=%f dy=%f dz=%f nb_head[i].x=%d nb_head[i].y=%d nb_head[i].z=%d nb_orig[i].x=%d nb_orig[i].y=%d nb_orig[i].z=%d\n",temp_dist,fr,dx,dy,dz,nb_head[i].x,nb_head[i].y,nb_head[i].z,nb_orig[i].x,nb_orig[i].y,nb_orig[i].z);
                                        nb_next_head[p].x=xp;
                                        nb_next_head[p].y=nb_head[i].y;
                                        nb_next_head[p].z=nb_head[i].z;
                                        nb_next_orig[p].x=nb_orig[i].x;
                                        nb_next_orig[p].y=nb_orig[i].y;
                                        nb_next_orig[p].z=nb_orig[i].z;
                                        g->matrix[xp][nb_head[i].y][nb_head[i].z].from=i;//index of nb_head
                                        g->matrix[xp][nb_head[i].y][nb_head[i].z].dist=temp_dist;
                                        p++;

                                        if(g->matrix[xp][nb_head[i].y][nb_head[i].z].phi==0 && temp_dist==fr)
                                                p--;	//doesn't count
                                        g->matrix[xp][nb_head[i].y][nb_head[i].z].phi=0;
                                }
                        }

                        if (g->matrix[xn][nb_head[i].y][nb_head[i].z].phi==1 ||
                           (g->matrix[xn][nb_head[i].y][nb_head[i].z].phi==0&&g->matrix[xn][nb_head[i].y][nb_head[i].z].from!=-1)){
                                if(g->matrix[xn][nb_head[i].y][nb_head[i].z].phi==1)	fr=PR2;
                                else	fr=g->matrix[xn][nb_head[i].y][nb_head[i].z].dist;

                                dx = g->matrix[xn][nb_head[i].y][nb_head[i].z].point.x - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.x;
                                dy = g->matrix[xn][nb_head[i].y][nb_head[i].z].point.y - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.y;
                                dz = g->matrix[xn][nb_head[i].y][nb_head[i].z].point.z - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.z;
                                temp_dist=dx*dx+dy*dy+dz*dz;


                                if((fr-temp_dist)>=0.01){
                                //fprintf(stderr,"2 temp_dist=%f fr=%f dx=%f dy=%f dz=%f nb_head[i].x=%d nb_head[i].y=%d nb_head[i].z=%d nb_orig[i].x=%d nb_orig[i].y=%d nb_orig[i].z=%d\n",temp_dist,fr,dx,dy,dz,nb_head[i].x,nb_head[i].y,nb_head[i].z,nb_orig[i].x,nb_orig[i].y,nb_orig[i].z);
                                        nb_next_head[p].x=xn;
                                        nb_next_head[p].y=nb_head[i].y;
                                        nb_next_head[p].z=nb_head[i].z;
                                        nb_next_orig[p].x=nb_orig[i].x;
                                        nb_next_orig[p].y=nb_orig[i].y;
                                        nb_next_orig[p].z=nb_orig[i].z;
                                        g->matrix[xn][nb_head[i].y][nb_head[i].z].from=i;//index of nb_head
                                        g->matrix[xn][nb_head[i].y][nb_head[i].z].dist=temp_dist;
                                        p++;

                                        if(g->matrix[xn][nb_head[i].y][nb_head[i].z].phi==0 && temp_dist==fr)
                                                p--;
                                        g->matrix[xn][nb_head[i].y][nb_head[i].z].phi=0;
                                }
                        }

                        if (g->matrix[nb_head[i].x][yp][nb_head[i].z].phi==1 ||
                           (g->matrix[nb_head[i].x][yp][nb_head[i].z].phi==0&&g->matrix[nb_head[i].x][yp][nb_head[i].z].from!=-1)){
                                if(g->matrix[nb_head[i].x][yp][nb_head[i].z].phi==1)	fr=PR2;
                                else	fr=g->matrix[nb_head[i].x][yp][nb_head[i].z].dist;

                                dx = g->matrix[nb_head[i].x][yp][nb_head[i].z].point.x - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.x;
                                dy = g->matrix[nb_head[i].x][yp][nb_head[i].z].point.y - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.y;
                                dz = g->matrix[nb_head[i].x][yp][nb_head[i].z].point.z - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.z;
                                temp_dist=dx*dx+dy*dy+dz*dz;

                                if((fr-temp_dist)>=0.01){
                                //fprintf(stderr,"3 temp_dist=%f fr=%f dx=%f dy=%f dz=%f nb_head[i].x=%d nb_head[i].y=%d nb_head[i].z=%d nb_orig[i].x=%d nb_orig[i].y=%d nb_orig[i].z=%d\n",temp_dist,fr,dx,dy,dz,nb_head[i].x,nb_head[i].y,nb_head[i].z,nb_orig[i].x,nb_orig[i].y,nb_orig[i].z);
                                        nb_next_head[p].x=nb_head[i].x;
                                        nb_next_head[p].y=yp;
                                        nb_next_head[p].z=nb_head[i].z;
                                        nb_next_orig[p].x=nb_orig[i].x;
                                        nb_next_orig[p].y=nb_orig[i].y;
                                        nb_next_orig[p].z=nb_orig[i].z;
                                        g->matrix[nb_head[i].x][yp][nb_head[i].z].from=i;//index of nb_head
                                        g->matrix[nb_head[i].x][yp][nb_head[i].z].dist=temp_dist;
                                        p++;
                                        if(g->matrix[nb_head[i].x][yp][nb_head[i].z].phi==0 && temp_dist==fr)
                                                p--;
                                        g->matrix[nb_head[i].x][yp][nb_head[i].z].phi=0;
                                }
                        }

                        if (g->matrix[nb_head[i].x][yn][nb_head[i].z].phi==1 ||
                           (g->matrix[nb_head[i].x][yn][nb_head[i].z].phi==0&&g->matrix[nb_head[i].x][yn][nb_head[i].z].from!=-1)){
                                if(g->matrix[nb_head[i].x][yn][nb_head[i].z].phi==1)	fr=PR2;
                                else	fr=g->matrix[nb_head[i].x][yn][nb_head[i].z].dist;

                                dx = g->matrix[nb_head[i].x][yn][nb_head[i].z].point.x - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.x;
                                dy = g->matrix[nb_head[i].x][yn][nb_head[i].z].point.y - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.y;
                                dz = g->matrix[nb_head[i].x][yn][nb_head[i].z].point.z - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.z;
                                temp_dist=dx*dx+dy*dy+dz*dz;

                                if((fr-temp_dist)>=0.01){
                                //fprintf(stderr,"4 temp_dist=%f fr=%f dx=%f dy=%f dz=%f nb_head[i].x=%d nb_head[i].y=%d nb_head[i].z=%d \n",temp_dist,fr,dx,dy,dz,nb_head[i].x,nb_head[i].y,nb_head[i].z);
                                        nb_next_head[p].x=nb_head[i].x;
                                        nb_next_head[p].y=yn;
                                        nb_next_head[p].z=nb_head[i].z;
                                        nb_next_orig[p].x=nb_orig[i].x;
                                        nb_next_orig[p].y=nb_orig[i].y;
                                        nb_next_orig[p].z=nb_orig[i].z;
                                        g->matrix[nb_head[i].x][yn][nb_head[i].z].from=i;//index of nb_head
                                        g->matrix[nb_head[i].x][yn][nb_head[i].z].dist=temp_dist;
                                        p++;
                                        if(g->matrix[nb_head[i].x][yn][nb_head[i].z].phi==0 && temp_dist==fr)
                                                p--;
                                        g->matrix[nb_head[i].x][yn][nb_head[i].z].phi=0;
                                }
                        }

                        if (g->matrix[nb_head[i].x][nb_head[i].y][zp].phi==1 ||
                           (g->matrix[nb_head[i].x][nb_head[i].y][zp].phi==0&&g->matrix[nb_head[i].x][nb_head[i].y][zp].from!=-1)){
                                if(g->matrix[nb_head[i].x][nb_head[i].y][zp].phi==1)	fr=PR2;
                                else	fr=g->matrix[nb_head[i].x][nb_head[i].y][zp].dist;

                                dx = g->matrix[nb_head[i].x][nb_head[i].y][zp].point.x - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.x;
                                dy = g->matrix[nb_head[i].x][nb_head[i].y][zp].point.y - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.y;
                                dz = g->matrix[nb_head[i].x][nb_head[i].y][zp].point.z - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.z;
                                temp_dist=dx*dx+dy*dy+dz*dz;

                                if((fr-temp_dist)>=0.01){
                                //fprintf(stderr,"5 temp_dist=%f fr=%f dx=%f dy=%f dz=%f nb_head[i].x=%d nb_head[i].y=%d nb_head[i].z=%d \n",temp_dist,fr,dx,dy,dz,nb_head[i].x,nb_head[i].y,nb_head[i].z);
                                        nb_next_head[p].x=nb_head[i].x;
                                        nb_next_head[p].y=nb_head[i].y;
                                        nb_next_head[p].z=zp;
                                        nb_next_orig[p].x=nb_orig[i].x;
                                        nb_next_orig[p].y=nb_orig[i].y;
                                        nb_next_orig[p].z=nb_orig[i].z;
                                        g->matrix[nb_head[i].x][nb_head[i].y][zp].from=i;//index of nb_head
                                        g->matrix[nb_head[i].x][nb_head[i].y][zp].dist=temp_dist;
                                        p++;
                                        if(g->matrix[nb_head[i].x][nb_head[i].y][zp].phi==0 && temp_dist==fr)
                                                p--;
                                        g->matrix[nb_head[i].x][nb_head[i].y][zp].phi=0;
                                }
                        }

                        if (g->matrix[nb_head[i].x][nb_head[i].y][zn].phi==1 ||
                           (g->matrix[nb_head[i].x][nb_head[i].y][zn].phi==0&&g->matrix[nb_head[i].x][nb_head[i].y][zn].from!=-1)){
                                if(g->matrix[nb_head[i].x][nb_head[i].y][zn].phi==1)	fr=PR2;
                                else	fr=g->matrix[nb_head[i].x][nb_head[i].y][zn].dist;

                                dx = g->matrix[nb_head[i].x][nb_head[i].y][zn].point.x - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.x;
                                dy = g->matrix[nb_head[i].x][nb_head[i].y][zn].point.y - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.y;
                                dz = g->matrix[nb_head[i].x][nb_head[i].y][zn].point.z - g->matrix[nb_orig[i].x][nb_orig[i].y][nb_orig[i].z].point.z;
                                temp_dist=dx*dx+dy*dy+dz*dz;

                                if((fr-temp_dist)>=0.01){
                                	//fprintf(stderr,"6 temp_dist=%f fr=%f dx=%f dy=%f dz=%f nb_head[i].x=%d nb_head[i].y=%d nb_head[i].z=%d \n",temp_dist,fr,dx,dy,dz,nb_head[i].x,nb_head[i].y,nb_head[i].z);
                                        nb_next_head[p].x=nb_head[i].x;
                                        nb_next_head[p].y=nb_head[i].y;
                                        nb_next_head[p].z=zn;
                                        nb_next_orig[p].x=nb_orig[i].x;
                                        nb_next_orig[p].y=nb_orig[i].y;
                                        nb_next_orig[p].z=nb_orig[i].z;
                                        g->matrix[nb_head[i].x][nb_head[i].y][zn].from=i;//index of nb_head
                                        g->matrix[nb_head[i].x][nb_head[i].y][zn].dist=temp_dist;
                                        p++;
                                        if(g->matrix[nb_head[i].x][nb_head[i].y][zn].phi==0 && temp_dist==fr)
                                                p--;
                                        g->matrix[nb_head[i].x][nb_head[i].y][zn].phi=0;
                                }
                        }
                }
                free(nb_head);
                free(nb_orig);
                nb_head = nb_next_head;
                nb_orig = nb_next_orig;
                co=p;

        }
        free(nb_head);
        free(nb_orig);

}

//a2fVertexOffset lists the positions, relative to vertex0, of each of the 8 vertices of a cube
static const float a2fVertexOffset[8][3] =
{
        {0.0, 0.0, 0.0},{1.0, 0.0, 0.0},{1.0, 1.0, 0.0},{0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},{1.0, 0.0, 1.0},{1.0, 1.0, 1.0},{0.0, 1.0, 1.0}
};

//a2iEdgeConnection lists the index of the endpoint vertices for each of the 12 edges of the cube
static const int a2iEdgeConnection[12][2] =
{
        {0,1}, {1,2}, {2,3}, {3,0},
        {4,5}, {5,6}, {6,7}, {7,4},
        {0,4}, {1,5}, {2,6}, {3,7}
};

//a2fEdgeDirection lists the direction vector (vertex1-vertex0) for each edge in the cube
static const float a2fEdgeDirection[12][3] =
{
        {1.0, 0.0, 0.0},{0.0, 1.0, 0.0},{-1.0, 0.0, 0.0},{0.0, -1.0, 0.0},
        {1.0, 0.0, 0.0},{0.0, 1.0, 0.0},{-1.0, 0.0, 0.0},{0.0, -1.0, 0.0},
        {0.0, 0.0, 1.0},{0.0, 0.0, 1.0},{ 0.0, 0.0, 1.0},{0.0,  0.0, 1.0}
};

// For any edge, if one vertex is inside of the surface and the other is outside of the surface
//  then the edge intersects the surface
// For each of the 8 vertices of the cube can be two possible states : either inside or outside of the surface
// For any cube the are 2^8=256 possible sets of vertex states
// This table lists the edges intersected by the surface for all 256 possible vertex states
// There are 12 edges.  For each entry in the table, if edge #n is intersected, then bit #n is set to 1

int aiCubeEdgeFlags[256]=
{
        0x000, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c, 0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
        0x190, 0x099, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c, 0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
        0x230, 0x339, 0x033, 0x13a, 0x636, 0x73f, 0x435, 0x53c, 0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
        0x3a0, 0x2a9, 0x1a3, 0x0aa, 0x7a6, 0x6af, 0x5a5, 0x4ac, 0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
        0x460, 0x569, 0x663, 0x76a, 0x066, 0x16f, 0x265, 0x36c, 0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
        0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0x0ff, 0x3f5, 0x2fc, 0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
        0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x055, 0x15c, 0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
        0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0x0cc, 0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
        0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc, 0x0cc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
        0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c, 0x15c, 0x055, 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
        0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc, 0x2fc, 0x3f5, 0x0ff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
        0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c, 0x36c, 0x265, 0x16f, 0x066, 0x76a, 0x663, 0x569, 0x460,
        0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac, 0x4ac, 0x5a5, 0x6af, 0x7a6, 0x0aa, 0x1a3, 0x2a9, 0x3a0,
        0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c, 0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x033, 0x339, 0x230,
        0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c, 0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x099, 0x190,
        0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c, 0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x000
};

//  For each of the possible vertex states listed in aiCubeEdgeFlags there is a specific triangulation
//  of the edge intersection points.  a2iTriangleConnectionTable lists all of them in the form of
//  0-5 edge triples with the list terminated by the invalid value -1.
//  For example: a2iTriangleConnectionTable[3] list the 2 triangles formed when corner[0]
//  and corner[1] are inside of the surface, but the rest of the cube is not.
//
//  I found this table in an example program someone wrote long ago.  It was probably generated by hand

int a2iTriangleConnectionTable[256][16] =
{
        {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
        {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
        {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
        {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
        {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
        {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
        {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
        {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
        {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
        {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
        {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
        {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
        {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
        {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
        {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
        {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
        {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
        {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
        {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
        {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
        {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
        {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
        {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
        {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
        {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
        {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
        {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
        {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
        {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
        {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
        {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
        {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
        {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
        {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
        {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
        {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
        {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
        {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
        {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
        {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
        {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
        {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
        {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
        {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
        {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
        {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
        {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
        {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
        {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
        {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
        {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
        {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
        {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
        {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
        {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
        {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
        {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
        {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
        {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
        {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
        {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
        {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
        {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
        {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
        {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
        {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
        {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
        {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
        {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
        {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
        {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
        {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
        {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
        {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
        {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
        {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
        {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
        {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
        {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
        {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
        {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
        {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
        {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
        {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
        {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
        {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
        {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
        {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
        {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
        {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
        {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
        {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
        {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
        {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
        {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
        {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
        {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
        {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
        {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
        {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
        {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
        {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
        {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
        {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
        {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
        {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
        {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
        {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
        {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
        {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
        {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
        {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
        {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
        {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
        {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
        {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
        {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
        {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
        {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
        {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
        {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
        {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
        {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
        {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
        {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
        {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
        {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
        {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
        {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
        {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
        {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
        {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
        {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
        {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
        {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
        {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
        {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
        {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
        {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
        {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
        {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
        {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
        {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
        {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
        {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
        {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
        {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
        {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
        {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
        {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
        {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
        {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
        {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
        {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
        {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
        {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
        {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
        {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
        {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
        {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
        {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
        {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
        {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
        {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
        {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
        {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
        {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
        {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
        {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
        {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
        {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
        {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
        {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
        {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
        {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
        {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
};

//void SceneBuilder::buildCATrace(Molecule& mol)
//{
//	int i;
//	float color[4] = {1.0f,1.0f,1.0f,1.0f};
//	glLineWidth(1.0);
//	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, color);
//	glBegin(GL_LINES);
//		for (i=0;i<mol.npoints-1;i++)
//		{
//			glVertex3f(mol.xpoints[i],mol.ypoints[i],mol.zpoints[i]);
//			glVertex3f(mol.xpoints[i+1],mol.ypoints[i+1],mol.zpoints[i+1]);
//		}
//	glEnd();
//}
//
//void SceneBuilder::buildScene(Grid& grid)
//{
//	float color2[4] = {0.4f,1.0f,0.4f,1.0f};
//	//float color2[4] = {1.0f,1.0f,1.0f,1.0f};
//
//	glColorMaterial(GL_FRONT,GL_DIFFUSE);
//	glEnable(GL_COLOR_MATERIAL);
//	//glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, color2);
//	//glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color2);
//	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, color2);
//	glBegin(GL_TRIANGLES);
//		marchingCube(grid);
//		//showGrid(grid);
//	glEnd();
//}

void vNormalizeVector(GLvector* rfVectorResult, GLvector rfVectorSource)
{
        float fOldLength;
        float fScale;

        fOldLength = sqrt( (rfVectorSource.fX * rfVectorSource.fX) +
                            (rfVectorSource.fY * rfVectorSource.fY) +
                            (rfVectorSource.fZ * rfVectorSource.fZ) );

        if(fOldLength == 0.0)
        {
                rfVectorResult->fX = rfVectorSource.fX;
                rfVectorResult->fY = rfVectorSource.fY;
                rfVectorResult->fZ = rfVectorSource.fZ;
        }
        else
        {
                fScale = 1.0f/fOldLength;
                rfVectorResult->fX = rfVectorSource.fX*fScale;
                rfVectorResult->fY = rfVectorSource.fY*fScale;
                rfVectorResult->fZ = rfVectorSource.fZ*fScale;
        }
}

//vGetNormal() finds the gradient of the scalar field at a point
//This gradient can be used as a very accurate vertx normal for lighting calculations
void vGetNormal(GLvector *rfNormal, S_Grid grid, int x, int y, int z)
{
	int x0,x1,y0,y1,z0,z1;
	x0 = x-1;
	if (x0<0) x0=0;
	x1 = x+1;
	if (x1>=grid->N) x1=(grid->N)-1;
	rfNormal->fX = (grid->matrix[x0][y][z].phi?1.0f:-1.0f) - (grid->matrix[x1][y][z].phi?1.0f:-1.0f);
	//rfNormal->fX = grid.matrix[x0][y][z].phi - grid.matrix[x1][y][z].phi;
	y0 = y-1;
	if (y0<0) y0=0;
	y1 = y+1;
	if (y1>=grid->N) y1=(grid->N)-1;
	rfNormal->fY = (grid->matrix[x][y0][z].phi?1.0f:-1.0f) - (grid->matrix[x][y1][z].phi?1.0f:-1.0f);
	//rfNormal->fY = grid.matrix[x][y0][z].phi - grid.matrix[x][y1][z].phi;
	z0 = z-1;
	if (z0<0) z0=0;
	z1 = z+1;
	if (z1>=grid->N) z1=(grid->N)-1;
	rfNormal->fZ = (grid->matrix[x][y][z0].phi?1.0f:-1.0f) - (grid->matrix[x][y][z1].phi?1.0f:-1.0f);
	//rfNormal->fZ = grid.matrix[x][y][z0].phi - grid.matrix[x][y][z1].phi;

	vNormalizeVector(rfNormal, *rfNormal);
}

//fGetOffset finds the approximate point of intersection of the surface
// between two points with the values fValue1 and fValue2
float fGetOffset(float fValue1, float fValue2, float fValueDesired)
{
        double fDelta = fValue2 - fValue1;

        if(fDelta == 0.0)
        {
                return 0.5;
        }
        return (float)((fValueDesired - fValue1)/fDelta);
}

//vMarchCube performs the Marching Cubes algorithm on a single cube


//marchingCube iterates over the entire dataset, calling vMarchCube on each cube
void marchingCube(S_Grid grid)
{
	int len = grid->N;
	int count = 0;

		       int iCorner, iVertex, iVertexTest, iEdge, iTriangle, iFlagIndex, iEdgeFlags;
				int ix,iy,iz;

		        float fOffset;
		        char afCubeValue[8];
		        GLvector asEdgeVertex[12];
		        GLvector asEdgeNorm[12];
		        char vertice[40],vertice2[512],  lineaV[512], lineaN[512];
		        strcpy(lineaV,""); strcpy(lineaN,""); strcpy(vertice,""); strcpy(vertice2,"");

		        FILE *f;
		        if ( (f=fopen("surf.vmd", "w"))==NULL) {
		        fprintf(stderr, "\n  Error->Cannot open file \n");
		        exit(1);
		        }

	for (int x=0;x<len;x++)
	for (int y=0;y<len;y++)
	for (int z=0;z<len;z++)
        {

			int x0,x1,y0,y1,z0,z1;
			x0 = x-1;
			if (x0<0) x0=0;
			x1 = x+1;
			if (x1>=grid->N) x1=(grid->N)-1;
			y0 = y-1;
			if (y0<0) y0=0;
			y1 = y+1;
			if (y1>=grid->N) y1=(grid->N)-1;
			z0 = z-1;
			if (z0<0) z0=0;
			z1 = z+1;
			if (z1>=grid->N) z1=(grid->N)-1;

			afCubeValue[0] = (grid->matrix[x][y][z].phi?1:-1);
			afCubeValue[1] = (grid->matrix[x1][y][z].phi?1:-1);
			afCubeValue[2] = (grid->matrix[x1][y1][z].phi?1:-1);
			afCubeValue[3] = (grid->matrix[x][y1][z].phi?1:-1);
			afCubeValue[4] = (grid->matrix[x][y][z1].phi?1:-1);
			afCubeValue[5] = (grid->matrix[x1][y][z1].phi?1:-1);
			afCubeValue[6] = (grid->matrix[x1][y1][z1].phi?1:-1);
			afCubeValue[7] = (grid->matrix[x][y1][z1].phi?1:-1);

			/*afCubeValue[0] = grid.matrix[x][y][z].phi;
			afCubeValue[1] = grid.matrix[x1][y][z].phi;
			afCubeValue[2] = grid.matrix[x1][y1][z].phi;
			afCubeValue[3] = grid.matrix[x][y1][z].phi;
			afCubeValue[4] = grid.matrix[x][y][z1].phi;
			afCubeValue[5] = grid.matrix[x1][y][z1].phi;
			afCubeValue[6] = grid.matrix[x1][y1][z1].phi;
			afCubeValue[7] = grid.matrix[x][y1][z1].phi;*/

		        //Find which vertices are inside of the surface and which are outside
		        iFlagIndex = 0;
		        for(iVertexTest = 0; iVertexTest < 8; iVertexTest++)
		        {
		                if(afCubeValue[iVertexTest] <= 0)
		                        iFlagIndex |= 1<<iVertexTest;
		        }

		        //Find which edges are intersected by the surface
		        iEdgeFlags = aiCubeEdgeFlags[iFlagIndex];

		        //If the cube is entirely inside or outside of the surface, then there will be no intersections
		        if(iEdgeFlags != 0)
		        {



		        //Find the point of intersection of the surface with each edge
		        //Then find the normal to the surface at those points
		        for(iEdge = 0; iEdge < 12; iEdge++)
		        {
		                //if there is an intersection on this edge
		                if((iEdgeFlags & (1<<iEdge))!=0)
		                {
		                        fOffset = fGetOffset(afCubeValue[ a2iEdgeConnection[iEdge][0] ], afCubeValue[ a2iEdgeConnection[iEdge][1] ], 0.0f);

		                        asEdgeVertex[iEdge].fX = grid->matrix[x][y][z].point.x + (a2fVertexOffset[ a2iEdgeConnection[iEdge][0] ][0]  +  fOffset * a2fEdgeDirection[iEdge][0]) * grid->stepSize;
		                        asEdgeVertex[iEdge].fY = grid->matrix[x][y][z].point.y + (a2fVertexOffset[ a2iEdgeConnection[iEdge][0] ][1]  +  fOffset * a2fEdgeDirection[iEdge][1]) * grid->stepSize;
		                        asEdgeVertex[iEdge].fZ = grid->matrix[x][y][z].point.z + (a2fVertexOffset[ a2iEdgeConnection[iEdge][0] ][2]  +  fOffset * a2fEdgeDirection[iEdge][2]) * grid->stepSize;

						if (asEdgeVertex[iEdge].fX > (grid->matrix[x][y][z].point.x+grid->matrix[x1][y][z].point.x)/2)
							ix = x1;
						else if (asEdgeVertex[iEdge].fX < (grid->matrix[x][y][z].point.x+grid->matrix[x0][y][z].point.x)/2)
							ix = x0;
						else ix = x;
						if (asEdgeVertex[iEdge].fY > (grid->matrix[x][y][z].point.y+grid->matrix[x][y1][z].point.y)/2)
							iy = y1;
						else if (asEdgeVertex[iEdge].fY < (grid->matrix[x][y][z].point.y+grid->matrix[x][y0][z].point.y)/2)
							iy = y0;
						else iy = y;
						if (asEdgeVertex[iEdge].fZ > (grid->matrix[x][y][z].point.z+grid->matrix[x][y][z1].point.z)/2)
							iz = z1;
						else if (asEdgeVertex[iEdge].fZ < (grid->matrix[x][y][z].point.z+grid->matrix[x][y][z0].point.z)/2)
							iz = z0;
						else iz = z;
		                        vGetNormal(&(asEdgeNorm[iEdge]), grid, ix, iy, iz);
		                }
		        }

		        int countV=0;
		        //Draw the triangles that were found.  There can be up to five per cube
		        for(iTriangle = 0; iTriangle < 5; iTriangle++)
		        {
		                if(a2iTriangleConnectionTable[iFlagIndex][3*iTriangle] < 0)
		                        break;

		                for(iCorner = 2; iCorner >=0 ; iCorner--)
		                //for(iCorner = 0; iCorner <3 ; iCorner++)
		                {
		                        iVertex = a2iTriangleConnectionTable[iFlagIndex][3*iTriangle+iCorner];

		                       if (countV==0) fprintf(f,"draw triangle ");
		                       if (countV==3) {
		                    	   fprintf(f,"\ndraw triangle ");
		                    	   fprintf(f,"{%1.2f %1.2f %1.2f}  ", asEdgeVertex[iVertex].fX, asEdgeVertex[iVertex].fY, asEdgeVertex[iVertex].fZ);
		                    	   countV=1;
		                       } else {
		                       fprintf(f,"{%1.2f %1.2f %1.2f}  ", asEdgeVertex[iVertex].fX, asEdgeVertex[iVertex].fY, asEdgeVertex[iVertex].fZ);
		                       countV++;}


		                 //       sprintf(vertice,"%.2f,%.2f,%.2f,\n", asEdgeNorm[iVertex].fX,   -asEdgeNorm[iVertex].fY,   -asEdgeNorm[iVertex].fZ);
		                 //       strcat(lineaN, vertice);

		                 //       sprintf(vertice,"%.5f,%.5f,%.5f\n", asEdgeVertex[iVertex].fX, asEdgeVertex[iVertex].fY, asEdgeVertex[iVertex].fZ);
		                 //       strcat(lineaV, vertice);
		       //                 glNormal3f(-asEdgeNorm[iVertex].fX,   -asEdgeNorm[iVertex].fY,   -asEdgeNorm[iVertex].fZ);
		       //                 glVertex3f(asEdgeVertex[iVertex].fX, asEdgeVertex[iVertex].fY, asEdgeVertex[iVertex].fZ);
					count++;

		                }
		        }
		        fprintf(f,"\n");
		    //    if (count >0) {
		    //    fprintf(stderr,"\n count %d\n", count);
		    //    fprintf(stderr,"%s\n", lineaN);
		    //    fprintf(stderr,"%s\n", lineaV);
		    //    }

          }
		}



	fclose(f);
	printf("lsms> Number of triangles = %d\n",count/3);

}

