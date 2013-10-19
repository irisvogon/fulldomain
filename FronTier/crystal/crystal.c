/************************************************************************************
FronTier is a set of libraries that implements differnt types of Front Traking algorithms.
Front Tracking is a numerical method for the solution of partial differential equations 
whose solutions have discontinuities.  


Copyright (C) 1999 by The University at Stony Brook. 
 

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

******************************************************************************/


/*
*				example2.c:
*
*		User initialization example for Front Package:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	This is example of three circles all moving the a normal velocity.
*	Bifurcation occurs when they meet each other. FronTier solves
*	the bifurcation automatically.
*
*/

#include <FronTier.h>

/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/


/********************************************************************
 *	Velocity function parameters for the front	 	    *
 ********************************************************************/
#define		MAX_NUM_VERTEX_IN_CELL		20

typedef struct {
        double center[3];
        double radius;
} TEST_SPHERE_PARAMS;

struct _SOLUTE_PARAMS
{
        int dim;
        double coeff;
	double C_0;	/* solute concentration in ambient state */
	double C_eq;	/* solute concentration in equilibrium with solid */
	double rho_s;	/* density of the precipitated solid phase */
	double D;	/* diffusion coefficient of the solute concentration */
	double k;	/* local reaction rate coefficient */
	double *solute;
	COMPONENT comp_solute;
};
typedef struct _SOLUTE_PARAMS SOLUTE_PARAMS;

typedef double STATE;

	/*  Local Application Function Declarations */

static void 	solute_main_driver(Front*,TEST_SPHERE_PARAMS);
static void 	read_domain_infor_from_input(char*,F_BASIC_DATA*);
static void	read_crystal_params_from_input(char*,SOLUTE_PARAMS*);
static void	read_time_control_from_input(char*,Front*);
static void 	init_solute_function(double**,double**,Front*);
static void 	update_solute_function(double*,double*,Front*);
static void 	print_front_and_interior_state(const char*,Front*,double*);
static void 	read_print_front_and_interior_state(const char*,Front*,double**,
				double**);
static void	solute_point_propagate(Front*,POINTER,POINT*,POINT*,	
			HYPER_SURF_ELEMENT*,HYPER_SURF*,double,double*);
static double 	solute_func(double*,COMPONENT,double);
static double   solute_of_state(POINTER);
static double 	crystal_curve(POINTER,double*);
static double 	sine_curve(POINTER,double*);
static double 	sphere_surf(POINTER,double*);
static bool 	fractal_dimension(Front*,TEST_SPHERE_PARAMS,double*,double*);

char *in_name,*restart_state_name,*restart_name,*out_name;
bool RestartRun;
bool ReadFromInput;
int RestartStep;
bool binary = YES;

#define		CRYSTAL_COMP		1
#define		SOLUTE_COMP		2

int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	SOLUTE_PARAMS solute_params; /* velocity function parameters */
	TEST_SPHERE_PARAMS s_params;

	FrontInitStandardIO(argc,argv,&f_basic);

	/* Initialize basic computational data */

	f_basic.size_of_intfc_state = sizeof(STATE);

	in_name      		= f_basic.in_name;
	restart_state_name      = f_basic.restart_state_name;
        out_name     		= f_basic.out_name;
        restart_name 		= f_basic.restart_name;
        RestartRun   		= f_basic.RestartRun;
        ReadFromInput   	= f_basic.ReadFromInput;
	RestartStep 		= f_basic.RestartStep;

	sprintf(restart_state_name,"%s-state.ts%s",restart_name,
			right_flush(RestartStep,7));
	sprintf(restart_name,"%s.ts%s",restart_name,right_flush(RestartStep,7));
#if defined(__MPI__)
        sprintf(restart_name,"%s-nd%s",restart_name,right_flush(pp_mynode(),4));
        sprintf(restart_state_name,"%s-nd%s",restart_state_name,
			right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */
	if (ReadFromInput)
	{
	    read_domain_infor_from_input(in_name,&f_basic);
	}
	else
	{
	    int i;
	    if (f_basic.dim == 1)
            {
                f_basic.L[0] = 0.0;
                f_basic.U[0] = 1.0;
                f_basic.gmax[0] = 80;
                f_basic.boundary[0][0] = REFLECTION_BOUNDARY;
                f_basic.boundary[0][1] = REFLECTION_BOUNDARY;
            }
	    else
	    {
	    	for (i = 0; i < f_basic.dim; ++i)
	    	{
		    f_basic.L[i] = 0.0;	
		    f_basic.U[i] = 1.0;	
		    f_basic.gmax[i] = 80;	
	    	    f_basic.boundary[i][0] = f_basic.boundary[i][1] 
				= PERIODIC_BOUNDARY;
	    	}
	    }
	}

	FrontStartUp(&front,&f_basic);

	if (!RestartRun)
	{
	    /* Initialize interface through level function */

	    level_func_pack.neg_component = CRYSTAL_COMP;
	    level_func_pack.pos_component = SOLUTE_COMP;
	    if (f_basic.dim == 1)
	    {
	    	float **point;
                bi_array(&point,1,1,FLOAT);
                point[0][0] = 0.03;
                level_func_pack.num_points = 1;
                level_func_pack.point_array = point;
	    }
	    else if (f_basic.dim == 2)
	    {
		s_params.center[0] = 0.5;
	        s_params.center[1] = 0.5;
	        s_params.radius = 0.01;
		level_func_pack.func_params = (POINTER)&s_params;
	    	level_func_pack.func = crystal_curve;
	    }
	    else
	    {
		s_params.center[0] = 0.5;
	        s_params.center[1] = 0.5;
	        s_params.center[2] = 0.5;
	        s_params.radius = 0.03;
		level_func_pack.func_params = (POINTER)&s_params;
	    	level_func_pack.func = sphere_surf;
	    }

	    FrontInitIntfc(&front,&level_func_pack);
	}

	/* Initialize velocity field function */

	solute_params.dim = f_basic.dim;
	solute_params.coeff = 0.1;
	solute_params.comp_solute = SOLUTE_COMP;

	velo_func_pack.func_params = (POINTER)&solute_params;
	velo_func_pack.func = NULL;

	FrontInitVelo(&front,&velo_func_pack);

        /* For geometry-dependent velocity, use first
        * order point propagation function, higher order
        * propagation requires surface propagate, currently
        * in writing, not yet in use. The following override
        * the assigned fourth_order_point_propagate.
        */

        front._point_propagate = solute_point_propagate;

	/* Propagate the front */

	solute_main_driver(&front,s_params);

	clean_up(0);
}

static  void solute_main_driver(
        Front *front,
	TEST_SPHERE_PARAMS s_params)
{
        int ip,im,status,count;
        Front *newfront;
        double dt,dt_frac,CFL;
        bool is_print_time, is_movie_time, time_limit_reached;
        char s[10];
        double fcrds[MAXD];
        int  i,dim = front->rect_grid->dim;
	double *solute,*solute_new,*tmp_solute;		/* solutent */
	INTERFACE *grid_intfc;
	SOLUTE_PARAMS *solute_params;
	double h_min;
	double solute_time_limit; 
	double frac_dim,radius;
	FILE *Radius_file,*FracDim_file;
	bool bdry_reached = NO;
	/*TMP*/
	FILE *sfile = fopen("solute1d","w");

	Curve_redistribution_function(front) = expansion_redistribute;

	h_min = front->rect_grid->h[0];
	for (i = 1; i < dim; ++i)
	    h_min = min(h_min,front->rect_grid->h[i]);

	solute_params = (SOLUTE_PARAMS*)front->vparams;
	if (ReadFromInput)
	{
	    read_crystal_params_from_input(in_name,solute_params);
	    read_time_control_from_input(in_name,front);
            CFL = Time_step_factor(front);
	}
	else
	{
	    solute_params->D = 0.5;
	    solute_params->k = 80.0;
	    solute_params->C_eq = 1.0;
	    solute_params->C_0 = 2.0;
	    solute_params->rho_s = 2.0;
	    front->max_time = 0.1;
	    front->max_step = 10000;
	    front->print_time_interval = 0.05;
	    front->movie_frame_interval = 0.005;
	    Frequency_of_redistribution(front,GENERAL_WAVE) = 5;
            CFL = Time_step_factor(front) = 0.40;
	}

	solute_time_limit = 0.5*sqr(h_min)/solute_params->D/(double)dim;

	if (RestartRun)
	{
	    FrontScatIntfcBuffer(front);
	    read_print_front_and_interior_state(restart_state_name,front,
	    			&solute,&solute_new);
	}
	else
	{
            front->time = 0.0;
            front->dt = 0.0;
	    front->step = 0;

	    init_solute_function(&solute,&solute_new,front);
	}

	FrontScatExpDualGridArray(solute,front);
	solute_params->solute = solute;

        if (!RestartRun)
        {
            // Front standard output
            FrontPrintOut(front,out_name);
            FrontMovieFrame(front,out_name,binary);
	    if (dim != 1)
	    	bdry_reached = fractal_dimension(front,s_params,&frac_dim,
					&radius);
	    if (pp_mynode() == 0)
	    {
		char fname[200];
		sprintf(fname,"%s-radius",out_name);
	    	Radius_file = fopen(fname,"w");
		sprintf(fname,"%s-FracDim",out_name);
	    	FracDim_file = fopen(fname,"w");

		fprintf(Radius_file,"\"Crystal radius\"\n\n");
		fprintf(Radius_file,"%f  %f\n",front->time,radius);
		fprintf(FracDim_file,"\"Fractal dimension\"\n\n");
		fprintf(FracDim_file,"%f  %f\n",front->time,frac_dim);
		fflush(Radius_file);
		fflush(FracDim_file);
	    }

	    // Problem specific output
            print_front_and_interior_state(out_name,front,solute);

            status = FrontAdvance(front->dt,&dt_frac,front,&newfront,
                                (POINTER)NULL);
            ip = im = 1;
            front->dt = CFL*FrontHypTimeStep(front);
            front->dt = min(front->dt,CFL*solute_time_limit);
        }
        else
        {
            ip = (int)(front->time/front->print_time_interval + 1.0);
            im = (int)(front->time/front->movie_frame_interval + 1.0);
        }

        front->dt = FrontOutputTimeControl(front,
                        &is_movie_time,&is_print_time,
                        &time_limit_reached,&im,&ip);

        for (;;)
        {
	    solute_params->solute = solute;
            status = FrontAdvance(front->dt,&dt_frac,front,&newfront,
                                (POINTER)NULL);
	    FrontFreeGridIntfc(front);
            assign_interface_and_free_front(front,newfront);
	    update_solute_function(solute,solute_new,front);
	    FrontScatExpDualGridArray(solute_new,front);

	    tmp_solute = solute;
	    solute = solute_new;
	    solute_new = tmp_solute;		/* Switching new to current */

            ++front->step;
            front->time += front->dt;

            front->dt = CFL*FrontHypTimeStep(front);
            front->dt = min(front->dt,CFL*solute_time_limit);

            printf("\ntime = %f   step = %7d   dt = %f\n",
                        front->time,front->step,front->dt);
            fflush(stdout);

            if (is_print_time || time_limit_reached)
	    {
		// Front standard output
                FrontPrintOut(front,out_name);
		// Problem specific output
		print_front_and_interior_state(out_name,front,solute);
	    }
            if (is_movie_time || time_limit_reached)
		// Front standard output
                FrontMovieFrame(front,out_name,binary);
	    if (dim != 1)
	    	bdry_reached = fractal_dimension(front,s_params,&frac_dim,
					&radius);
	    else if (is_movie_time)	/*TMP*/
	    {
	    	fprintf(sfile,"\n\n");
                for (int i = 0; i < 80; ++i)
                    fprintf(sfile," %f  %f\n",i/80.0,solute[i]);
                fflush(sfile);
	    }
	    if (bdry_reached) time_limit_reached = YES;
	    if (pp_mynode() == 0)
	    {
		fprintf(Radius_file,"%f  %f\n",front->time,radius);
		fprintf(FracDim_file,"%f  %f\n",front->time,frac_dim);
		fflush(Radius_file);
		fflush(FracDim_file);
	    }

            if (time_limit_reached)
                    break;
	    /* Output section, next dt may be modified */

            front->dt = FrontOutputTimeControl(front,
                        &is_movie_time,
                        &is_print_time,
                        &time_limit_reached,
                        &im,&ip);
        }
	if (pp_mynode() == 0)
	{
	    fclose(Radius_file);
	    fclose(FracDim_file);
	}
}       /* end solute_main_driver */

LOCAL double crystal_curve(
        POINTER func_params,
        double *coords)
{

        TEST_SPHERE_PARAMS *s_params = (TEST_SPHERE_PARAMS*)func_params;
        double dist, theta;
	double *cen = s_params->center;
	double radius = s_params->radius;

        dist =   sqrt(sqr(coords[0]-cen[0]) + sqr(coords[1]-cen[1]));
        theta = asin(fabs(coords[1]-0.5)/dist);
	if (coords[0]-0.5 < 0 && coords[1]-0.5 > 0)
	    theta = PI - theta;
	else if (coords[0]-0.5 < 0 && coords[1]-0.5 < 0)
	    theta = PI + theta;
	else if (coords[0]-0.5 > 0 && coords[1]-0.5 < 0)
	    theta = 2*PI - theta;
        return dist - radius + .003*sin(6.0*theta);
}       /* end crystal_curve */

static	void init_solute_function(
	double **psolute_old,
	double **psolute_new,
	Front *front)
{
	INTERFACE *intfc = front->interf;
	INTERFACE *grid_intfc;
	RECT_GRID *grid;
	int i,j,k,dim,num_nodes,*gmax;
	double *L,*h,coords[MAXD];
	COMPONENT comp,*gr_comp;
	Table *T;
	double *solute_old,*solute_new;
	STATE *sl,*sr;
	POINT *p;
	HYPER_SURF *hs;
	HYPER_SURF_ELEMENT *hse;
	SOLUTE_PARAMS *solute_params = (SOLUTE_PARAMS*)front->vparams;

	/* Initialize states at the interface */
	next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
	{
	    slsr(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    if (negative_component(hs) == SOLUTE_COMP)
	    {
	    	*sl = (STATE)solute_params->C_0;
	    }
	    else if (positive_component(hs) == SOLUTE_COMP)
	    {
	    	*sr = (STATE)solute_params->C_0;
	    }
	}

	/* Initialize states in the interior regions */

	FrontMakeGridIntfc(front);
        grid_intfc = front->grid_intfc;
	grid = &topological_grid(grid_intfc);
	dim = grid->dim;
	T = table_of_interface(grid_intfc);
	gr_comp = T->components;
	num_nodes = 1;
	gmax = grid->gmax;
	L = grid->L;
	h = grid->h;
	for (i = 0; i < dim; ++i)
	    num_nodes *= (gmax[i] + 1);

	solute_old = (double*)malloc(num_nodes*sizeof(double));
	solute_new = (double*)malloc(num_nodes*sizeof(double));
	switch (dim)
	{
        case 1:
            for (i = 0; i <= gmax[0]; ++i)
            {
                comp = gr_comp[d_index1d(i,gmax)];
                coords[0] = L[0] + i*h[0];
                solute_old[d_index1d(i,gmax)] = solute_func(coords,comp,
                                        (STATE)solute_params->C_0);
                solute_new[d_index1d(i,gmax)] = solute_func(coords,comp,
                                        (STATE)solute_params->C_0);
            }
            break;
	case 2:
	    for (i = 0; i <= gmax[0]; ++i)
	    for (j = 0; j <= gmax[1]; ++j)
	    {
	    	comp = gr_comp[d_index2d(i,j,gmax)];
	    	coords[0] = L[0] + i*h[0];
	    	coords[1] = L[1] + j*h[1];
	    	solute_old[d_index2d(i,j,gmax)] = solute_func(coords,comp,
					(STATE)solute_params->C_0);
	    	solute_new[d_index2d(i,j,gmax)] = solute_func(coords,comp,
					(STATE)solute_params->C_0);
	    }
	    break;
	case 3:
	    for (i = 0; i <= gmax[0]; ++i)
	    for (j = 0; j <= gmax[1]; ++j)
	    for (k = 0; k <= gmax[2]; ++k)
	    {
	    	comp = gr_comp[d_index3d(i,j,k,gmax)];
	    	coords[0] = L[0] + i*h[0];
	    	coords[1] = L[1] + j*h[1];
	    	coords[2] = L[2] + k*h[2];
	    	solute_old[d_index3d(i,j,k,gmax)] = solute_func(coords,comp,
					(STATE)solute_params->C_0);
	    	solute_new[d_index3d(i,j,k,gmax)] = solute_func(coords,comp,
					(STATE)solute_params->C_0);
	    }
	    break;
	}
	*psolute_old = solute_old;
	*psolute_new = solute_new;
}	/* init_solute_function */

static	void update_solute_function(
	double *solute_old,
	double *solute_new,
	Front *front)
{
	INTERFACE *intfc = front->interf;
	INTERFACE *grid_intfc;
	RECT_GRID *grid;
	int i,j,k,l,m,ic,icn,*gmax,icoords[MAXD],nc,dim;
	int gmin[MAXD],ipn[MAXD];
	double *h,*L,coords[MAXD],crx_coords[MAXD];
	double  solute,solute_nb,dgrad[MAXD];
	Table *T;
	SOLUTE_PARAMS *solute_params = (SOLUTE_PARAMS*)front->vparams;
	double coef = solute_params->D*front->dt;
	COMPONENT comp,*gr_comp;
	bool fr_crx_grid_seg;
	const GRID_DIRECTION dir[3][2] = 
		{{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};

	start_clock("update_solute_function");
	FrontMakeGridIntfc(front);
        grid_intfc = front->grid_intfc;
	grid = &topological_grid(grid_intfc);
	T = table_of_interface(grid_intfc);
	dim = grid->dim;
	gmax = grid->gmax;
	h = grid->h;
	L = grid->L;
	for (i = 0; i < dim; ++i) gmin[i] = 0;

	gr_comp = T->components;
	switch (dim)
	{
        case 1:
            for (i = 1; i < gmax[0]; ++i)
            {
                 icoords[0] = i;
                 ic = d_index1d(i,gmax);
                 comp = gr_comp[ic];
                 if (comp != SOLUTE_COMP)
                 {
                     solute_new[ic] = 0.0;
                     continue;
                 }
                 solute =  solute_old[ic];
                 solute_new[ic] = solute;
                 for (l = 0; l < dim; ++l)
                 {
                     dgrad[l] = 0.0;
                     for (m = 0; m < 2; ++m)
                     {
                         fr_crx_grid_seg = FrontStateAtGridCrossing(front,
                                        icoords,dir[l][m],comp,solute_of_state,
                                        &solute_nb,crx_coords);
                         if (fr_crx_grid_seg)
                         {
                             dgrad[l] += (solute_nb - solute)/h[l];
                         }
                         else
                         {
                             next_ip_in_dir(icoords,dir[l][m],ipn,gmin,gmax);
                             icn = d_index1d(ipn[0],gmax);
                             dgrad[l] += (solute_old[icn] - solute)/h[l];
                         }
                     }
                     solute_new[ic] += coef*dgrad[l]/h[l];
                 }
            }
            break;
	case 2:
	    for (i = 1; i < gmax[0]; ++i)
	    for (j = 1; j < gmax[1]; ++j)
	    {
	    	icoords[0] = i;
	    	icoords[1] = j;
	    	ic = d_index2d(i,j,gmax);
	    	comp = gr_comp[ic];
	    	if (comp != SOLUTE_COMP) 
	    	{
		    solute_new[ic] = 0.0;
	    	    continue;
	        }
	        solute = solute_old[ic];
		solute_new[ic] = solute;
		for (l = 0; l < dim; ++l)
		{
	            dgrad[l] = 0.0;
                    for (m = 0; m < 2; ++m)
                    {
                        fr_crx_grid_seg = FrontStateAtGridCrossing(front,
                                        icoords,dir[l][m],comp,solute_of_state,
                                        &solute_nb,crx_coords);
                        if (fr_crx_grid_seg)
                        {
                            dgrad[l] += (solute_nb - solute)/h[l];
                        }
                        else
                        {
                            next_ip_in_dir(icoords,dir[l][m],ipn,gmin,gmax);
                            icn = d_index2d(ipn[0],ipn[1],gmax);
                            dgrad[l] += (solute_old[icn] - solute)/h[l];
                        }
                    }
		    solute_new[ic] += coef*dgrad[l]/h[l];
		}
	    }
	    break;
	case 3:
	    for (i = 1; i < gmax[0]; ++i)
	    for (j = 1; j < gmax[1]; ++j)
	    for (k = 1; k < gmax[2]; ++k)
	    {
	    	icoords[0] = i;
	    	icoords[1] = j;
	    	icoords[2] = k;
	    	ic = d_index3d(i,j,k,gmax);
	    	comp = gr_comp[ic];
	    	if (comp != SOLUTE_COMP) 
	    	{
		    solute_new[ic] = 0.0;
	    	    continue;
	        }
	        solute = solute_old[ic];
		solute_new[ic] = solute;
		for (l = 0; l < dim; ++l)
		{
	            dgrad[l] = 0.0;
                    for (m = 0; m < 2; ++m)
                    {
                        fr_crx_grid_seg = FrontStateAtGridCrossing(front,
                                        icoords,dir[l][m],comp,solute_of_state,
                                        &solute_nb,crx_coords);
                        if (fr_crx_grid_seg)
                        {
                            dgrad[l] += (solute_nb - solute)/h[l];
                        }
                        else
                        {
                            next_ip_in_dir(icoords,dir[l][m],ipn,gmin,gmax);
                            icn = d_index3d(ipn[0],ipn[1],ipn[2],gmax);
                            dgrad[l] += (solute_old[icn] - solute)/h[l];
                        }
                    }
		    solute_new[ic] += coef*dgrad[l]/h[l];
		}
	    }
	    break;
	}
	stop_clock("update_solute_function");
}	/* update_solute_function */

static  double solute_of_state(
        POINTER state)
{
        STATE *solute_state = (STATE*)state;
        return (double)(*solute_state);
}       /* end solute_of_state */

static	double solute_func(
	double *coords,
	COMPONENT comp,
	double C_0)
{
	if (comp != SOLUTE_COMP) return 0.0;
	return C_0;
}	/* end solute_init_func */

static  void solute_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        double              dt,
        double              *V)
{
        double vel[MAXD];
        int i, dim = front->rect_grid->dim;
	double nor[MAXD];
        double p1[MAXD];
        double *p0 = Coords(oldp);
        double dn,*h = front->rect_grid->h;
        double s0,s1,grad_s;
        STATE *sl,*sr,*state;
	SOLUTE_PARAMS *solute_params = (SOLUTE_PARAMS*)front->vparams;
        COMPONENT comp = solute_params->comp_solute;
        double *solute = solute_params->solute;
        double D = solute_params->D;
        double k = solute_params->k;
        double C_0 = solute_params->C_0;
        double C_eq = solute_params->C_eq;
        double rho_s = solute_params->rho_s;

        if (wave_type(oldhs) < FIRST_PHYSICS_WAVE_TYPE)
        {
            for (i = 0; i < dim; ++i)
                Coords(newp)[i] = Coords(oldp)[i];
            slsr(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
            state = (negative_component(oldhs) == comp) ?
                        (STATE*)left_state(newp) :
                        (STATE*)right_state(newp);
            *state = (negative_component(oldhs) == comp) ? *sl : *sr;
            return;
        }

        GetFrontNormal(oldp,oldhse,oldhs,nor,front);
        if (negative_component(oldhs) == comp)
            for (i = 0; i < dim; ++i)
                nor[i] *= -1.0;
        dn = grid_size_in_direction(nor,h,dim);
        for (i = 0; i < dim; ++i)
            p1[i] = p0[i] + nor[i]*dn;
        slsr(oldp,oldhse,oldhs,(POINTER*)&sl,(POINTER*)&sr);
        state = (negative_component(oldhs) == comp) ? sl : sr;
        s0 = (double)(*state);

	if (!FrontIntrpStateAtPoint(front,comp,p1,solute,solute_of_state,&s1))
        {
            s1 = (double)(*state);
        }

        grad_s = (s1 - s0)/dn;
        for (i = 0; i < dim; ++i)
        {
            //vel[i] = nor[i]*k*(s0 - C_eq)/rho_s;
            vel[i] = nor[i]*D*grad_s/rho_s;
        }
        for (i = 0; i < dim; ++i)
        {
            Coords(newp)[i] = Coords(oldp)[i] + dt*vel[i];
            set_max_front_speed(i,fabs(vel[i]),NULL,Coords(newp),front);
        }

	/* Update the state of the new interface point */
	state = (negative_component(oldhs) == comp) ? (STATE*)left_state(newp) :
			(STATE*)right_state(newp);
	s0 += dt*(D*grad_s - k*(s0 - C_eq))/dn;
        *state = (STATE)s0;
}       /* solute_point_propagate */


LOCAL double sine_curve(
        POINTER func_params,
        double *coords)
{

        double dist, theta;

	double y_intfc = 0.5 + 0.2*sin(6.0*PI*coords[0]);

        dist =   y_intfc - coords[1];
        return dist;
}       /* end sine_curve */


static	void print_front_and_interior_state(
	const char *out_name,
	Front *front,
	double *solute)
{
	INTERFACE *intfc = front->interf;
	INTERFACE *grid_intfc;
	RECT_GRID *grid;
	int i,j,k,dim,*gmax;
	double *L,*h,coords[MAXD];
	Table *T;
	STATE *sl,*sr;
	POINT *p;
	HYPER_SURF *hs;
	HYPER_SURF_ELEMENT *hse;
	char filename[100];
	FILE *outfile;
	bool grid_intfc_made = NO;

	sprintf(filename,"%s-state.ts%s",out_name,right_flush(front->step,7));
#if defined(__MPI__)
        sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */
	outfile = fopen(filename,"w");
	
	/* Initialize states at the interface */
	fprintf(outfile,"Interface solute state:\n");
	int count = 0;
	next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
	{
	    count++;
	    slsr(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    fprintf(outfile,"%20.16f %20.16f\n",(double)*sl,(double)*sr);
	}

	/* Initialize states in the interior regions */

	if (front->grid_intfc == NULL)
	{
	    FrontMakeGridIntfc(front);
	    grid_intfc = front->grid_intfc;
	    grid_intfc_made = YES;
	}
	else
	    grid_intfc = front->grid_intfc;
	grid = &topological_grid(grid_intfc);
	dim = grid->dim;
	gmax = grid->gmax;
	L = grid->L;
	h = grid->h;

	fprintf(outfile,"\nInterior solute state:\n");
	switch (dim)
	{
	case 2:
	    for (i = 0; i <= gmax[0]; ++i)
	    for (j = 0; j <= gmax[1]; ++j)
	    {
	        fprintf(outfile,"%20.16f\n",solute[d_index2d(i,j,gmax)]);
	    }
	    break;
	case 3:
	    for (i = 0; i <= gmax[0]; ++i)
	    for (j = 0; j <= gmax[1]; ++j)
	    for (k = 0; k <= gmax[2]; ++k)
	    {
	        fprintf(outfile,"%20.16f\n",solute[d_index3d(i,j,k,gmax)]);
	    }
	}
	if (grid_intfc_made)
	    FrontFreeGridIntfc(front);
}	/* print_front_and_interior_state */

static	void read_print_front_and_interior_state(
	const char *in_name,
	Front *front,
	double **psolute_old,
	double **psolute_new)
{
	INTERFACE *intfc = front->interf;
	INTERFACE *grid_intfc;
	RECT_GRID *grid;
	int i,j,k,dim,*gmax,num_nodes;
	double *L,*h,coords[MAXD];
	Table *T;
	STATE *sl,*sr;
	POINT *p;
	HYPER_SURF *hs;
	HYPER_SURF_ELEMENT *hse;
	FILE *infile;
	double x,*solute,*solute_new;

	infile = fopen(in_name,"r");

	/* Initialize states at the interface */
	next_output_line_containing_string(infile,"Interface solute state:");
	next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
	{
	    slsr(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    fscanf(infile,"%lf",&x);
	    *sl = (STATE)x;
	    fscanf(infile,"%lf",&x);
	    *sr = (STATE)x;
	}

	/* Initialize states in the interior regions */

	FrontMakeGridIntfc(front);
	grid_intfc = front->grid_intfc;
	grid = &topological_grid(grid_intfc);
	dim = grid->dim;
	gmax = grid->gmax;
	L = grid->L;
	h = grid->h;
	num_nodes = 1;
	for (i = 0; i < dim; ++i)
	    num_nodes *= (gmax[i] + 1);

	solute = (double*)malloc(num_nodes*sizeof(double));
	solute_new = (double*)malloc(num_nodes*sizeof(double));
	next_output_line_containing_string(infile,"Interior solute state:");
	switch (dim)
	{
	case 2:
	    for (i = 0; i <= gmax[0]; ++i)
	    for (j = 0; j <= gmax[1]; ++j)
	    {
	    	fscanf(infile,"%lf",&solute[d_index2d(i,j,gmax)]);
	    }
	    break;
	case 3:
	    for (i = 0; i <= gmax[0]; ++i)
	    for (j = 0; j <= gmax[1]; ++j)
	    for (k = 0; k <= gmax[2]; ++k)
	    {
	    	fscanf(infile,"%lf",&solute[d_index3d(i,j,k,gmax)]);
	    }
	}
	*psolute_old = solute;
	*psolute_new = solute_new;
	fclose(infile);
}	/* end read_print_front_and_interior_state */


static double sphere_surf(
        POINTER func_params,
        double *coords)
{
        TEST_SPHERE_PARAMS *s_params = (TEST_SPHERE_PARAMS*)func_params;
	double x0,y0,z0,R;
	double distance;

        x0 = s_params->center[0];
        y0 = s_params->center[1];
        z0 = s_params->center[2];
	R = s_params->radius;

	distance = sqrt(sqr(coords[0] - x0) + sqr(coords[1] - y0) +
			sqr(coords[2] - z0)) - R;

        return distance;

}       /* end sphere_surf */



static	void read_domain_infor_from_input(
	char *in_name,
	F_BASIC_DATA *f_basic)
{
	FILE *infile;
	char input_string[200],sbdry[200];
	int i;

	infile = fopen(in_name,"r");
	for (i = 0; i < f_basic->dim; ++i)
	{
	    sprintf(input_string,"Domain limit in %d-th dimension:",i);
	    fgetstring(infile,input_string);
	    fscanf(infile,"%lf %lf",&f_basic->L[i],&f_basic->U[i]);
	}
	fgetstring(infile,"Computational grid:");
	for (i = 0; i < f_basic->dim; ++i)
	{
	    fscanf(infile,"%d",&f_basic->gmax[i]);
	}
	for (i = 0; i < f_basic->dim; ++i)
	{
	    sprintf(input_string,"Lower boundary in %d-th dimension:",i);
	    fgetstring(infile,input_string);
	    fscanf(infile,"%s",sbdry);
	    switch (sbdry[0])
	    {
	    case 'P':
	    case 'p':
	    	f_basic->boundary[i][0] = PERIODIC_BOUNDARY;
		break;
	    case 'D':
	    case 'd':
	    	f_basic->boundary[i][0] = DIRICHLET_BOUNDARY;
	    	break;
	    case 'N':
	    case 'n':
	    	f_basic->boundary[i][0] = NEUMANN_BOUNDARY;
	    	break;
	    default:
	    	printf("Unknown boundary!\n");
		clean_up(ERROR);
	    }
	    sprintf(input_string,"Upper boundary in %d-th dimension:",i);
	    fgetstring(infile,input_string);
	    fscanf(infile,"%s",sbdry);
	    switch (sbdry[0])
	    {
	    case 'P':
	    case 'p':
	    	f_basic->boundary[i][1] = PERIODIC_BOUNDARY;
		break;
	    case 'D':
	    case 'd':
	    	f_basic->boundary[i][1] = DIRICHLET_BOUNDARY;
	    	break;
	    case 'N':
	    case 'n':
	    	f_basic->boundary[i][1] = NEUMANN_BOUNDARY;
	    	break;
	    default:
	    	printf("Unknown boundary!\n");
		clean_up(ERROR);
	    }
	}
	fclose(infile);
}

static	void	read_crystal_params_from_input(
	char *in_name,
	SOLUTE_PARAMS *solute_params)
{
	FILE *infile;

	infile = fopen(in_name,"r");
	fgetstring(infile,"Diffusion coefficient:");
	fscanf(infile,"%lf",&solute_params->D);
	fgetstring(infile,"Growth rate:");
	fscanf(infile,"%lf",&solute_params->k);
	fgetstring(infile,"Equilibrium concentration:");
	fscanf(infile,"%lf",&solute_params->C_eq);
	fgetstring(infile,"Ambient concentration:");
	fscanf(infile,"%lf",&solute_params->C_0);
	fgetstring(infile,"Crystal density:");
	fscanf(infile,"%lf",&solute_params->rho_s);
	fclose(infile);
}

static	void	read_time_control_from_input(
	char *in_name,
	Front *front)
{
	FILE *infile;

	infile = fopen(in_name,"r");
	fgetstring(infile,"Max time:");
	fscanf(infile,"%lf",&front->max_time);
	fgetstring(infile,"Max step:");
	fscanf(infile,"%d",&front->max_step);
	fgetstring(infile,"Print interval:");
	fscanf(infile,"%lf",&front->print_time_interval);
	fgetstring(infile,"Movie frame interval:");
	fscanf(infile,"%lf",&front->movie_frame_interval);
	fgetstring(infile,"CFL factor:");
	fscanf(infile,"%lf",&(Time_step_factor(front)));
	fgetstring(infile,"Redistribution interval:");
	fscanf(infile,"%d",&(Frequency_of_redistribution(front,GENERAL_WAVE)));
	fclose(infile);
}

static bool fractal_dimension(
	Front *front,
	TEST_SPHERE_PARAMS s_params,
	double *frac_dim,
	double *radius)
{
	double coords[MAXD],*center = s_params.center;
	double dist,r_sqr,r_max,r_min = s_params.radius;
	INTERFACE *grid_intfc,*intfc = front->interf;
	int i,j,k,*gmax,dim = intfc->dim;
	POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	struct Table *T;
	RECT_GRID *grid;
	double *L,*U,*h;
	COMPONENT *gr_comp,comp;
	int N,Nc;
	bool grid_intfc_made = NO;
	double ratio;
	bool bdry_reached;
	double margin[MAXD];

	/* Find maximum radius of crystal growth */
	r_max = 0.0;
	grid = computational_grid(intfc);
	L = grid->GL;
	U = grid->GU;
	for (i = 0; i < dim; ++i)
	    margin[i] = 0.01*(U[i] - L[i]);
	bdry_reached = NO;
	next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    if (wave_type(hs) < FIRST_PHYSICS_WAVE_TYPE)
	    	continue;
	    r_sqr = 0.0;
	    for (i = 0; i < dim; ++i)
	    {
		if (Coords(p)[i] >= U[i] - margin[i] ||
		    Coords(p)[i] <= L[i] + margin[i])
		    bdry_reached = YES;
	    	r_sqr += sqr(Coords(p)[i] - center[i]);
	    }
	    if (r_max < r_sqr) r_max = r_sqr;
	}
	r_max = sqrt(r_max);
	*radius = r_max;
#if defined (__MPI__)
	pp_global_max(radius,1);
#endif /* defined (__MPI__)

	/* Preparation for counting */

        if (front->grid_intfc == NULL)
        {
            FrontMakeGridIntfc(front);
            grid_intfc = front->grid_intfc;
            grid_intfc_made = YES;
        }
        else
            grid_intfc = front->grid_intfc;
        grid = &topological_grid(grid_intfc);
        gmax = grid->gmax;
        L = grid->L;
        h = grid->h;
	T = table_of_interface(grid_intfc);
        gr_comp = T->components;

	/* Start counting */
	N = Nc = 0;
	switch (dim)
	{
	case 2:
	    for (i = 0; i <= gmax[0]; ++i)
            for (j = 0; j <= gmax[1]; ++j)
            {
	    	comp = gr_comp[d_index2d(i,j,gmax)];
		coords[0] = L[0] + i*h[0];
                coords[1] = L[1] + j*h[1];
		dist = sqrt(sqr(coords[0] - center[0]) +
			    sqr(coords[1] - center[1]));
	    	if (dist > r_min && dist < r_max)
		{
		    ++N;
		    if (comp == CRYSTAL_COMP)
		    	++Nc;
		}
	    }
	    break;
	case 3:
	    for (i = 0; i <= gmax[0]; ++i)
            for (j = 0; j <= gmax[1]; ++j)
	    for (k = 0; k <= gmax[2]; ++k)
            {
	    	comp = gr_comp[d_index3d(i,j,k,gmax)];
		coords[0] = L[0] + i*h[0];
                coords[1] = L[1] + j*h[1];
		coords[2] = L[2] + k*h[2];
		dist = sqrt(sqr(coords[0] - center[0]) +
			    sqr(coords[1] - center[1]) +
			    sqr(coords[2] - center[2]));
	    	if (dist > r_min && dist < r_max)
		{
		    ++N;
		    if (comp == CRYSTAL_COMP)
		    	++Nc;
		}
	    }
	}
#if defined (__MPI__)
	pp_global_isum(&N,1);
	pp_global_isum(&Nc,1);
#endif /* defined (__MPI__) */
	if (grid_intfc_made)
	    FrontFreeGridIntfc(front);
	ratio = ((double)N)/((double)Nc);
	*frac_dim = (double)dim + log(ratio)/log(h[0]);
	return bdry_reached;
}	/* end fractal_dimension */
	
