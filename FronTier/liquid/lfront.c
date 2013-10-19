/***************************************************************************************************************************
 *		FRONTIER_EB.c
 ***************************************************************************************************************************/
#include <stdio.h>
#include <math.h>
#include <FronTier.h>
#include "lfront.h"


#define RADIUS		0.5
#define MAX_ITER	50

FRONTIER_EB::FRONTIER_EB(Front &front):m_pFront(&front)
{
}

void FRONTIER_EB::getCrossing(double *pa, double *pb, double *crossing, double *nor)
{
	double a[2], b[2];
	
	a[0] = pa[0];	a[1] = pa[1];
	b[0] = pb[0];	b[1] = pb[1];
	
	int ca, cb, cc;
	ca = getComponent(a);
	cb = getComponent(b);
	
	for(int i=0; i<MAX_ITER; i++)
	{
	    crossing[0] = 0.5*(a[0]+b[0]);
	    crossing[1] = 0.5*(a[1]+b[1]);
	    cc = getComponent(crossing);
	    if(cc==ca)
	    {
	        a[0] = crossing[0];
		a[1] = crossing[1];
	    }
	    else
	    {
		b[0] = crossing[0];
		b[1] = crossing[1];
	    }
	}	

	double dist = sqrt(crossing[0]*crossing[0]+crossing[1]*crossing[1]);
	nor[0] = crossing[0]/dist;
	nor[1] = crossing[1]/dist;

	double center[2] = {0,0};
	
	if(getComponent(center)!=getComponent(pa))
	{
	    nor[0] = -nor[0];
	    nor[1] = -nor[1];
	}
}

// input : p[]
// output: q[]
//
void FRONTIER_EB::getNearestInterfacePoint(double *p, double *q)
{
	INTERFACE *intfc = m_pFront->interf;
	double t;
	HYPER_SURF_ELEMENT *phse;
	HYPER_SURF *phs;
	nearest_interface_point(p,getComponent(p),intfc,NO_BOUNDARIES,
				NULL,q,&t,&phse,&phs);
}

int FRONTIER_EB::getComponent(double *p)
{
	INTERFACE *intfc = m_pFront->interf;
	int comp = component(p,intfc);
	return comp; 
}

void FRONTIER_EB::save(char *filename)
{
	
	RECT_GRID *rect_grid = m_pFront->rect_grid;
	INTERFACE *intfc    = m_pFront->interf;
		
	int i, j;
	int xmax = rect_grid->gmax[0];
	int ymax = rect_grid->gmax[1];
	double x, y;
	
	FILE *hfile = fopen(filename, "w");
	if(hfile==NULL)
	{
		printf("\n can't open %s in SaveAsTecplot_rect_grid_and_interface().", filename);
		exit(0);
	}
	
	// secondly print out the interface
		
	if(exists_interface(intfc))
	{
		CURVE		**cur;
		CURVE		*curve;
		BOND		*bond;
			
		for(cur=intfc->curves; cur && *cur; cur++)	// very strange and important here
		{
			curve = *cur;
			fprintf(hfile, "ZONE I=%d J=%d F=POINT \n", curve->num_points, 1);
			bond=curve->first;
			fprintf(hfile, "%.4f %.4f \n", bond->start->_coords[0], bond->start->_coords[1]);
			for(bond=curve->first; bond!=NULL; bond=bond->next)
				fprintf(hfile, "%.4f %.4f \n", bond->end->_coords[0], bond->end->_coords[1]);					
		}					
	}		
	// close the file handler
	fclose(hfile);
}



