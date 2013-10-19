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
        float center[3];
        float radius;
} TEST_SPHERE_PARAMS;

struct _SOLUTE_PARAMS
{
        int dim;
        float coeff;
	float C_eq;	/* solute concentration in equilibrium with solid */
	float rho_s;	/* density of the precipitated solid phase */
	float D;	/* diffusion coefficient of the solute concentration */
	float k;	/* local reaction rate coefficient */
	float *solute;
	COMPONENT comp_solute;
};
typedef struct _SOLUTE_PARAMS SOLUTE_PARAMS;

typedef float STATE;

	/*  Local Application Function Declarations */

static void 	solute_main_driver(Front*);
static void 	init_solute_function(float**,float**,Front*);
static void 	update_solute_function(float*,float*,Front*);
static void 	collect_cell_ptst(INTRP_CELL*,int*,COMPONENT,INTERFACE*,float*);
static float 	solute_func(float*,COMPONENT);
static float 	crystal_curve(POINTER,float*);
static float 	sine_curve(POINTER,float*);
static float 	solute_from_crxing(int*,GRID_DIRECTION,COMPONENT,INTERFACE*);
static float 	crxing_state(CRXING*,COMPONENT);
static void 	show_intfc_solute_states(INTERFACE*);
static void 	s_hyp_solution(float*,COMPONENT,HYPER_SURF*,SIDE,Front*,
				POINTER,POINTER,POINTER);
static void	solute_point_propagate(Front*,POINTER,POINT*,POINT*,	
			HYPER_SURF_ELEMENT*,HYPER_SURF*,float,float*);
static bool 	test_point_in_tri(float*,float**);
static bool 	test_point_in_tetra(float*,float**);
static bool 	build_linear_element(INTRP_CELL*,float*);
static float 	lin_intrp_var(float*,INTRP_CELL*);
static void 	print_front_and_interior_state(const char*,Front*,float*);
static void 	read_print_front_and_interior_state(const char*,Front*,float**,
				float**);
static float 	sphere_func(POINTER,float*);
static	void 	read_domain_infor_from_input(char*,F_BASIC_DATA*);
static	void	read_crystal_params_from_input(char*,SOLUTE_PARAMS*);
static	void	read_time_control_from_input(char*,Front*);

#define Cross_prod(x0,y0,x1,y1,xm,ym)                                   \
        ((x0 - xm)*(y1 - ym) - (x1 - xm)*(y0 - ym))

char *in_name,*restart_state_name,*restart_name,*out_name;
bool RestartRun;
bool ReadFromInput;
int RestartStep;
bool binary;

static float lin_cell_tol;

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
	    for (i = 0; i < f_basic.dim; ++i)
	    {
		f_basic.L[i] = 0.0;	
		f_basic.U[i] = 1.0;	
		f_basic.gmax[i] = 80;	
	    	f_basic.boundary[i][0] = f_basic.boundary[i][1] 
				= PERIODIC_BOUNDARY;
	    }
	}

	FrontStartUp(&front,&f_basic);

	/* Start testing */

	if (!RestartRun)
	    front.step = 0;

	if (!RestartRun)
	{
	    /* Initialize interface through level function */

	    level_func_pack.neg_component = 1;
	    level_func_pack.pos_component = 2;
	    if (f_basic.dim == 2)
	    	level_func_pack.func = crystal_curve;
	    else
	    {
		s_params.center[0] = 0.5;
	        s_params.center[1] = 0.5;
	        s_params.center[2] = 0.5;
	        s_params.radius = 0.1;
		level_func_pack.func_params = (POINTER)&s_params;
	    	level_func_pack.func = sphere_func;
	    }

	    FrontInitIntfc(&front,&level_func_pack);
	}

	/* Initialize velocity field function */

	solute_params.dim = f_basic.dim;
	solute_params.coeff = 0.1;
	solute_params.comp_solute = 2;

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

	solute_main_driver(&front);

	clean_up(0);
}

static  void solute_main_driver(
        Front *front)
{
        int ip,im,status,count;
        Front *newfront;
        float dt,dt_frac,CFL;
        bool is_print_time, is_movie_time, time_limit_reached;
        char s[10];
        float fcrds[MAXD];
        int  dim = front->rect_grid->dim;
	float *solute,*solute_new,*tmp_solute;		/* solutent */
	INTERFACE *grid_intfc;
	SOLUTE_PARAMS *solute_params;
	float h_min = min(front->rect_grid->h[0],front->rect_grid->h[1]);
	float solute_time_limit; 

	Curve_redistribution_function(front) = expansion_redistribute;

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
	    solute_params->rho_s = 2.0;
	    front->max_time = 0.5;
	    front->max_step = 2000;
	    front->print_time_interval = 0.05;
	    front->movie_frame_interval = 0.005;
	    Frequency_of_redistribution(front,GENERAL_WAVE) = 5;
            CFL = Time_step_factor(front) = 0.40;
	}

	solute_time_limit = 0.5*sqr(h_min)/solute_params->D/(float)dim;

	if (RestartRun)
	{
	    scatter_front(front);
	    read_print_front_and_interior_state(restart_state_name,front,
	    			&solute,&solute_new);
	}
	else
	{
            redistribute(front,YES,NO);
            front->time = 0.0;
            front->dt = 0.0;

	    init_solute_function(&solute,&solute_new,front);
	}

	scatter_top_grid_float_array(solute,front);
	solute_params->solute = solute;

        if (!RestartRun)
        {
            // Front standard output
            FrontPrintOut(front,out_name);
            FrontMovieFrame(front,out_name,binary);

	    // Problem specific output
            print_front_and_interior_state(out_name,front,solute);

            status = FrontAdvance(front->dt,&dt_frac,front,&newfront,
                                (POINTER)NULL);
            ip = im = 1;
            front->dt = CFL*FrontHypTimeStep(front);
            front->dt = min(front->dt,CFL*solute_time_limit);
            front->dt = FrontOutputTimeControl(front,
                        &is_movie_time,
                        &is_print_time,
                        &time_limit_reached,
                        &im,&ip);
        }
        else
        {
            ip = (int)(front->time/front->print_time_interval + 1.0);
            im = (int)(front->time/front->movie_frame_interval + 1.0);
        }

        for (;;)
        {
	    solute_params->solute = solute;
            status = FrontAdvance(front->dt,&dt_frac,front,&newfront,
                                (POINTER)NULL);
	    set_copy_intfc_states(YES);
	    free_grid_intfc(front->grid_intfc);
            assign_interface_and_free_front(front,newfront);
	    update_solute_function(solute,solute_new,front);
	    scatter_top_grid_float_array(solute_new,front);

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

            if (time_limit_reached)
                    break;
	    /* Output section, next dt may be modified */

            front->dt = FrontOutputTimeControl(front,
                        &is_movie_time,
                        &is_print_time,
                        &time_limit_reached,
                        &im,&ip);
        }
}       /* end solute_main_driver */

LOCAL float crystal_curve(
        POINTER func_params,
        float *coords)
{

        float dist, theta;

        dist =   sqrt(sqr(coords[0]-0.5) + sqr(coords[1]-0.5));
        theta = asin(fabs(coords[1]-0.5)/dist);
	if (coords[0]-0.5 < 0 && coords[1]-0.5 > 0)
	    theta = PI - theta;
	else if (coords[0]-0.5 < 0 && coords[1]-0.5 < 0)
	    theta = PI + theta;
	else if (coords[0]-0.5 > 0 && coords[1]-0.5 < 0)
	    theta = 2*PI - theta;
        return dist - .03 + .005*sin(6.0*theta);
}       /* end crystal_curve */

static	void init_solute_function(
	float **psolute_old,
	float **psolute_new,
	Front *front)
{
	INTERFACE *intfc = front->interf;
	INTERFACE *grid_intfc;
	RECT_GRID *grid;
	int i,j,k,dim,num_nodes,*gmax;
	float *L,*h,coords[MAXD];
	COMPONENT comp,*gr_comp;
	Table *T;
	float *solute_old,*solute_new;
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
	    if (negative_component(hs) == 2)
	    {
	    	*sl = (STATE)solute_params->C_eq;
	    }
	    else if (positive_component(hs) == 2)
	    {
	    	*sr = (STATE)solute_params->C_eq;
	    }
	}

	/* Initialize states in the interior regions */

	front->grid_intfc = grid_intfc = 
		make_grid_intfc(intfc,EXPANDED_DUAL_GRID,NULL);
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

	solute_old = (float*)malloc(num_nodes*sizeof(float));
	solute_new = (float*)malloc(num_nodes*sizeof(float));
	switch (dim)
	{
	case 2:
	    for (i = 0; i <= gmax[0]; ++i)
	    for (j = 0; j <= gmax[1]; ++j)
	    {
	    	comp = gr_comp[d_index2d(i,j,gmax)];
	    	coords[0] = L[0] + i*h[0];
	    	coords[1] = L[1] + j*h[1];
	    	solute_old[d_index2d(i,j,gmax)] = solute_func(coords,comp);
	    	solute_new[d_index2d(i,j,gmax)] = solute_func(coords,comp);
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
	    	solute_old[d_index3d(i,j,k,gmax)] = solute_func(coords,comp);
	    	solute_new[d_index3d(i,j,k,gmax)] = solute_func(coords,comp);
	    }
	    break;
	}
	*psolute_old = solute_old;
	*psolute_new = solute_new;
}	/* init_solute_function */

static	void update_solute_function(
	float *solute_old,
	float *solute_new,
	Front *front)
{
	INTERFACE *intfc = front->interf;
	INTERFACE *grid_intfc;
	RECT_GRID *grid;
	int i,j,k,l,m,ic,icn,*gmax,icoords[MAXD],nc,dim;
	int gmin[MAXD],ipn[MAXD];
	float *h,*L,coords[MAXD];
	float  solute,solute_nb,dgrad[MAXD];
	Table *T;
	SOLUTE_PARAMS *solute_params = (SOLUTE_PARAMS*)front->vparams;
	float coef = solute_params->D*front->dt;
	static CRXING *crxs[MAX_NUM_CRX];
	COMPONENT comp,*gr_comp;
	const GRID_DIRECTION dir[3][2] = 
		{{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};

	start_clock("update_solute_function");
	set_copy_intfc_states(YES);
	front->grid_intfc = grid_intfc = 
		make_grid_intfc(intfc,EXPANDED_DUAL_GRID,NULL);
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
	case 2:
	    for (i = 1; i < gmax[0]; ++i)
	    for (j = 1; j < gmax[1]; ++j)
	    {
	    	icoords[0] = i;
	    	icoords[1] = j;
	    	ic = d_index2d(i,j,gmax);
	    	comp = gr_comp[ic];
	    	if (comp != 2) 
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
		    	nc = GridSegCrossing(crxs,icoords,dir[l][m],grid_intfc);
			if (nc == 0)
			{
			    next_ip_in_dir(icoords,dir[l][m],ipn,gmin,gmax);
	    		    icn = d_index2d(ipn[0],ipn[1],gmax);
	    		    dgrad[l] += (solute_old[icn] - solute)/h[l];
			}
			else
			{
	    		    solute_nb = crxing_state(crxs[0],comp);
	    		    dgrad[l] += (solute_nb - solute)/h[l];
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
	    	if (comp != 2) 
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
		    	nc = GridSegCrossing(crxs,icoords,dir[l][m],grid_intfc);
			if (nc == 0)
			{
			    next_ip_in_dir(icoords,dir[l][m],ipn,gmin,gmax);
	    		    icn = d_index3d(ipn[0],ipn[1],ipn[2],gmax);
	    		    dgrad[l] += (solute_old[icn] - solute)/h[l];
			}
			else
			{
	    		    solute_nb = crxing_state(crxs[0],comp);
	    		    dgrad[l] += (solute_nb - solute)/h[l];
			}
		    }
		    solute_new[ic] += coef*dgrad[l]/h[l];
		}
	    }
	    break;
	}
	stop_clock("update_solute_function");
}	/* update_solute_function */

static	float solute_func(
	float *coords,
	COMPONENT comp)
{
	if (comp != 2) return 0.0;
	return 2.0;
}	/* end solute_init_func */

static float solute_from_crxing(
	int *icoords,
	GRID_DIRECTION dir,
	COMPONENT comp,
	INTERFACE *grid_intfc)
{
	CRXING *crx,**crxings;
	int nc;
	POINT *p;
	HYPER_SURF *hs;
	POINTER sl,sr;

	nc = GridSegCrossing(crxings,icoords,dir,grid_intfc);
	switch (dir)
	{
	case WEST:
	    crx = crxings[0];
	    break;
	case SOUTH:
	    crx = crxings[0];
	    break;
	case LOWER:
	    crx = crxings[0];
	    break;
	case EAST:
	    crx = crxings[nc-1];
	    break;
	case NORTH:
	    crx = crxings[nc-1];
	    break;
	case UPPER:
	    crx = crxings[nc-1];
	    break;
	}
	return crxing_state(crx,comp);
}	/* end solute_from_crxing */

static void s_hyp_solution(
	float *coords,
	COMPONENT comp,
	HYPER_SURF *hs,
	SIDE side,
	Front *front,
	POINTER solute,
	POINTER state,
	POINTER dflt_state)
{
	int icoords[MAXD];
	INTERFACE *grid_intfc = front->grid_intfc;
	float *solute_mesh = (float*)solute;
	static INTRP_CELL *blk_cell;
	RECT_GRID *gr = &topological_grid(grid_intfc);
	STATE *st;
	int i,dim = gr->dim;

	if (blk_cell == NULL)
	{
	    scalar(&blk_cell,sizeof(INTRP_CELL));
	    uni_array(&blk_cell->var,MAX_NUM_VERTEX_IN_CELL,sizeof(float));
	    bi_array(&blk_cell->coords,MAX_NUM_VERTEX_IN_CELL,MAXD,
	    		sizeof(float));
	    bi_array(&blk_cell->p_lin,MAXD+1,MAXD,sizeof(float));
	    uni_array(&blk_cell->var_lin,MAXD+1,sizeof(float));
	    lin_cell_tol = 1.0;
	    for (i = 0; i < dim; ++i)
	    	lin_cell_tol *= 0.00001*gr->h[i];
	}
	if (!rect_in_which(coords,icoords,gr))
	{
	    ft_assign(state,dflt_state,sizeof(STATE));
	    return;
	}
	collect_cell_ptst(blk_cell,icoords,comp,grid_intfc,solute_mesh);
	if (build_linear_element(blk_cell,coords))
	{
	    st = (STATE*)state;
	    *st = (STATE)lin_intrp_var(coords,blk_cell);
	}
	else
	{
	    if (!nearest_intfc_state(coords,comp,front->interf,state,NULL,NULL))
	    {
	    	STATE *dflt_st;
	    	st = (STATE*)state;
	    	dflt_st = (STATE*)dflt_state;
	    	*st = (STATE)(*dflt_st);
	    }
	}
}	/* end s_hyp_solution */


static bool build_linear_element(
	INTRP_CELL *blk_cell,
	float *coords)
{
	int dim = blk_cell->dim;
	float **ps = blk_cell->coords;
	float *vars = blk_cell->var;
	float **p = blk_cell->p_lin;
	float *var = blk_cell->var_lin;
	int i,j,k,l,nv = blk_cell->nv;
	float dist[MAX_NUM_VERTEX_IN_CELL];

	for (i = 0; i < nv; i++)
	{
	    dist[i] = 0.0;
	    for (j = 0; j < dim; ++j)
	    {
	   	dist[i] += sqr(coords[j] - ps[i][j]);
	    }
	}
	for (i = 0; i < nv; ++i)
	{
	    for (j = i+1; j < nv; ++j)
	    {
		if (dist[i] > dist[j])
		{
		    float tmp;
		    tmp = dist[j];
		    dist[j] = dist[i];
		    dist[i] = tmp;
		    tmp = vars[j];
		    vars[j] = vars[i];
		    vars[i] = tmp;
		    for (k = 0; k < dim; ++k)
		    {
		    	tmp = ps[j][k];
		    	ps[j][k] = ps[i][k];
		    	ps[i][k] = tmp;
		    }
		}
	    }
	}
	switch(dim)
	{
	case 2:
	    for (i = 0; i < nv; i++)
	    {
	    	p[0] = blk_cell->coords[i];
	    	var[0] =  blk_cell->var[i];
	    	for (j = i+1; j < nv; j++)
		{
	    	    p[1] = blk_cell->coords[j];
	    	    var[1] = blk_cell->var[j];
		    for (k = j+1; k < nv; k++)
		    {
	    	    	p[2] = blk_cell->coords[k];
	    	    	var[2] = blk_cell->var[k];
			if (test_point_in_tri(coords,p) == YES)
	    		    return FUNCTION_SUCCEEDED;
		    }
		}
	    }
	    break;
	case 3:
	    for (i = 0; i < nv; i++)
	    {
	    	p[0] = blk_cell->coords[i];
	    	var[0] = blk_cell->var[i];
	    	for (j = i+1; j < nv; j++)
		{
	    	    p[1] = blk_cell->coords[j];
	    	    var[1] = blk_cell->var[j];
		    for (k = j+1; k < nv; k++)
		    {
	    	    	p[2] = blk_cell->coords[k];
	    	    	var[2] = blk_cell->var[k];
			for (l = k+1; l < nv; l++)
			{
	    	    	    p[3] = blk_cell->coords[l];
	    	    	    var[3] = blk_cell->var[l];
			    if (test_point_in_tetra(coords,p) == YES)
			    {
	    		        return FUNCTION_SUCCEEDED;
			    }
			}
		    }
		}
	    }
	}
	return FUNCTION_FAILED;
}	/* end build_linear_element */

static void collect_cell_ptst(
	INTRP_CELL *blk_cell,
	int *icoords,
	COMPONENT comp,
	INTERFACE *grid_intfc,
	float *solute)
{
	Table *T = table_of_interface(grid_intfc);
	RECT_GRID *gr = &topological_grid(grid_intfc);
	int dim = gr->dim;
	int *gmax = gr->gmax;
	float *L = gr->L;
	float *h = gr->h;
	COMPONENT *gr_comp = T->components;
	static COMPONENT cell_comp2d[2][2];
	static COMPONENT cell_comp3d[2][2][2];
	int i,j,k,index,nv,nc;
	CRXING *crx,*crxs[MAX_NUM_CRX];
	GRID_DIRECTION dir;
	int ic[MAXD];
	
	blk_cell->is_bilinear = YES;
	blk_cell->dim = dim;
	nv = 0;
	switch (dim)
	{
	case 2:
	    for (i = 0; i < 2; ++i)
	    for (j = 0; j < 2; ++j)
	    {
	    	ic[0] = icoords[0] + i;
	    	ic[1] = icoords[1] + j;
	    	index = d_index2d(ic[0],ic[1],gmax);
	    	cell_comp2d[i][j] = gr_comp[index];
	    	if (gr_comp[index] == comp)
	    	{
	    	    blk_cell->coords[nv][0] = L[0] + ic[0]*h[0];
	    	    blk_cell->coords[nv][1] = L[1] + ic[1]*h[1];
		    blk_cell->var[nv] = solute[index];
		    nv++;
	    	}
	    	else
	    	    blk_cell->is_bilinear = NO;
	    }
	    break;
	case 3:
	    for (i = 0; i < 2; ++i)
	    for (j = 0; j < 2; ++j)
	    for (k = 0; k < 2; ++k)
	    {
	    	ic[0] = icoords[0] + i;
	    	ic[1] = icoords[1] + j;
	    	ic[2] = icoords[2] + k;
	    	index = d_index3d(ic[0],ic[1],ic[2],gmax);
	    	cell_comp3d[i][j][k] = gr_comp[index];
	    	if (gr_comp[index] == comp)
	    	{
	    	    blk_cell->coords[nv][0] = L[0] + ic[0]*h[0];
	    	    blk_cell->coords[nv][1] = L[1] + ic[1]*h[1];
	    	    blk_cell->coords[nv][2] = L[2] + ic[2]*h[2];
		    blk_cell->var[nv] = solute[index];
		    nv++;
	    	}
	    	else
	    	    blk_cell->is_bilinear = NO;
	    }
	    break;
	}
	if (blk_cell->is_bilinear == YES) 
	{
	    blk_cell->nv = nv;
	    return;
	}
	switch (dim)
	{
	case 2:
	    for (i = 0; i < 2; ++i)
	    for (j = 0; j < 2; ++j)
	    {
	    	ic[0] = icoords[0] + i;
	    	ic[1] = icoords[1] + j;
	    	if (cell_comp2d[i][j] == comp)
	    	{
	    	    if (cell_comp2d[(i+1)%2][j] != comp)
	    	    {
		    	dir = (i < (i+1)%2) ? EAST : WEST;
		    	nc = GridSegCrossing(crxs,ic,dir,grid_intfc);
		    	crx = crxs[0];
		    	if (nc == 0)
		    	{
		    	    screen("ERROR: no crxing between (%d %d) "
			           "and (%d %d)\n",icoords[0]+i,icoords[1]+j,
					icoords[0]+(i+1)%2,icoords[1]+j);
		    	}
		    	blk_cell->var[nv] = crxing_state(crx,comp);
		    	blk_cell->coords[nv][0] = Coords(crx->pt)[0];
		    	blk_cell->coords[nv][1] = Coords(crx->pt)[1];
		    	nv++;
	    	    }
	    	    if (cell_comp2d[i][(j+1)%2] != comp)
	    	    {
		    	dir = (j < (j+1)%2) ? NORTH : SOUTH;
		    	nc = GridSegCrossing(crxs,ic,dir,grid_intfc);
		    	crx = crxs[0];
		    	if (nc == 0)
		    	{
		    	    screen("ERROR: no crxing between (%d %d) "
			           "and (%d %d)\n",icoords[0]+i,icoords[1]+j,
					icoords[0]+i,icoords[1]+(j+1)%2);
		    	}
		    	blk_cell->var[nv] = crxing_state(crx,comp);
		    	blk_cell->coords[nv][0] = Coords(crx->pt)[0];
		    	blk_cell->coords[nv][1] = Coords(crx->pt)[1];
		    	nv++;
	    	    }
	    	}
	    }
	    break;
	case 3:
	    for (i = 0; i < 2; ++i)
	    for (j = 0; j < 2; ++j)
	    for (k = 0; k < 2; ++k)
	    {
	    	ic[0] = icoords[0] + i;
	    	ic[1] = icoords[1] + j;
	    	ic[2] = icoords[2] + k;
	    	if (cell_comp3d[i][j][k] == comp)
	    	{
	    	    if (cell_comp3d[(i+1)%2][j][k] != comp)
	    	    {
		    	dir = (i < (i+1)%2) ? EAST : WEST;
		    	nc = GridSegCrossing(crxs,ic,dir,grid_intfc);
		    	crx = crxs[0];
		    	if (nc == 0)
		    	{
		    	    screen("ERROR: no crxing between (%d %d %d) "
			           "and (%d %d %d)\n",icoords[0]+i,icoords[1]+j,
				   icoords[2]+k,icoords[0]+(i+1)%2,icoords[1]+j,
				   icoords[2]+k);
		    	}
		    	blk_cell->var[nv] = crxing_state(crx,comp);
		    	blk_cell->coords[nv][0] = Coords(crx->pt)[0];
		    	blk_cell->coords[nv][1] = Coords(crx->pt)[1];
		    	blk_cell->coords[nv][2] = Coords(crx->pt)[2];
		    	nv++;
	    	    }
	    	    if (cell_comp3d[i][(j+1)%2][k] != comp)
	    	    {
		    	dir = (j < (j+1)%2) ? NORTH : SOUTH;
		    	nc = GridSegCrossing(crxs,ic,dir,grid_intfc);
		    	crx = crxs[0];
		    	if (nc == 0)
		    	{
		    	    screen("ERROR: no crxing between (%d %d %d) "
			           "and (%d %d %d)\n",icoords[0]+i,icoords[1]+j,
				   icoords[2]+k,icoords[0]+i,icoords[1]+(j+1)%2,
				   icoords[2]+k);
		    	}
		    	blk_cell->var[nv] = crxing_state(crx,comp);
		    	blk_cell->coords[nv][0] = Coords(crx->pt)[0];
		    	blk_cell->coords[nv][1] = Coords(crx->pt)[1];
		    	blk_cell->coords[nv][2] = Coords(crx->pt)[2];
		    	nv++;
	    	    }
	    	    if (cell_comp3d[i][j][(k+1)%2] != comp)
	    	    {
		    	dir = (k < (k+1)%2) ? UPPER : LOWER;
		    	nc = GridSegCrossing(crxs,ic,dir,grid_intfc);
		    	crx = crxs[0];
		    	if (nc == 0)
		    	{
		    	    screen("ERROR: no crxing between (%d %d %d) "
			           "and (%d %d %d)\n",icoords[0]+i,icoords[1]+j,
				   icoords[2]+k,icoords[0]+i,icoords[1]+j,
				   icoords[2]+(k+1)%2);
		    	}
		    	blk_cell->var[nv] = crxing_state(crx,comp);
		    	blk_cell->coords[nv][0] = Coords(crx->pt)[0];
		    	blk_cell->coords[nv][1] = Coords(crx->pt)[1];
		    	blk_cell->coords[nv][2] = Coords(crx->pt)[2];
		    	nv++;
	    	    }
	    	}
	    }
	    break;
	}
	blk_cell->nv = nv;
}	/* end collect_cell_ptst */
	

static float crxing_state(
	CRXING *crx,
	COMPONENT comp)
{
	HYPER_SURF *hs = crx->hs;
	HYPER_SURF_ELEMENT *hse = Hyper_surf_element(crx->bond);
	POINT *p = crx->pt;
	POINTER sl,sr;
	STATE *lstate,*rstate;

	slsr(p,hse,hs,&sl,&sr);
	lstate = (STATE*)sl;
	rstate = (STATE*)sr;
	if (comp == negative_component(hs))
	    return (float)(*lstate);
	else if (comp == positive_component(hs))
	    return (float)(*rstate);
	else
	{
	    screen("ERROR: In crxing_state() component does not match\n");
	}
}	/* end crxing_state */

static void show_intfc_solute_states(
	INTERFACE *intfc)
{
	STATE *sl,*sr;
	POINT *p;
	HYPER_SURF *hs;
	HYPER_SURF_ELEMENT *hse;
	STATE *state;
	static int count = 0;
	char oname[100];
	FILE *ofile;

	sprintf(oname,"intfc-state.%d",count++);
	ofile = fopen(oname,"w");

	/* Initialize states at the interface */
	next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
	{
            if (wave_type(hs) < FIRST_PHYSICS_WAVE_TYPE)
		continue;
	    slsr(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
	    if (negative_component(hs) == 2)
	    {
		state = (STATE*)sl;
	    }
	    else if (positive_component(hs) == 2)
	    {
		state = (STATE*)sr;
	    }
	    fprintf(ofile,"%f %f\n",Coords(p)[0],(float)(*state));
	}
	fprintf(ofile,"\n\n");
	fclose(ofile);
}	/* end show_intfc_states */


static  void solute_point_propagate(
        Front *front,
        POINTER wave,
        POINT *oldp,
        POINT *newp,
        HYPER_SURF_ELEMENT *oldhse,
        HYPER_SURF         *oldhs,
        float              dt,
        float              *V)
{
        float vel[MAXD];
        int i, dim = front->rect_grid->dim;
	float nor[MAXD];
        float p1[MAXD];
        float *p0 = Coords(oldp);
        float dn,*h = front->rect_grid->h;
        float s0,s1,grad_s;
        STATE *sl,*sr,*state;
	SOLUTE_PARAMS *solute_params = (SOLUTE_PARAMS*)front->vparams;
        COMPONENT comp = solute_params->comp_solute;
        float *solute = solute_params->solute;
        float D = solute_params->D;
        float k = solute_params->k;
        float C_eq = solute_params->C_eq;
        float rho_s = solute_params->rho_s;
	STATE ans;

        if (wave_type(oldhs) < FIRST_PHYSICS_WAVE_TYPE)
        {
            for (i = 0; i < dim; ++i)
                Coords(newp)[i] = Coords(oldp)[i];
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
        s0 = (float)(*state);
        s_hyp_solution(p1,comp,oldhs,UNKNOWN_SIDE,front,(POINTER)solute,
                                (POINTER)&ans,state);
        s1 = (float)(ans);

	s0 = (k*C_eq + D*s1/dn)/(k + D/dn);
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
	s0 += dt*(D*grad_s - k*(s0 - C_eq));
        *state = (STATE)s0;
}       /* solute_point_propagate */

static 	bool test_point_in_tetra(
	float	*coords,
	float	**p)
{
	float		*p0, *p1, *p2, *p3;
	float		a[3],b[3],c[3],v[3];
	float		D,Dp;
	int		i,j;

	for (i = 0; i < 4; i++)
	{
	    p0 = p[i];
	    p1 = p[(i+1)%4];
	    p2 = p[(i+2)%4];
	    p3 = p[(i+3)%4];

	    for (j = 0; j < 3; j++)
	    {
	    	a[j] = p1[j] - p0[j];
	    	b[j] = p2[j] - p0[j];
	    	c[j] = p3[j] - p0[j];
	    	v[j] = coords[j] - p0[j];
	    }
	    D = Det3d(a,b,c);
	    if (fabs(D) < lin_cell_tol) return NO;
	    Dp = Det3d(a,b,v);
	    if ((D > 0.0 && Dp < 0.0) ||
	        (D < 0.0 && Dp > 0.0))
		return NO;
	}
	return YES;
}	/* end test_point_in_tetra */


static	bool	test_point_in_tri(
	float	*coords,
	float	**p)
{
	float	x0, y0, x1, y1, x2, y2, cp;

	x0 = p[0][0];	y0 = p[0][1];
	x1 = p[1][0];	y1 = p[1][1];
	x2 = p[2][0];	y2 = p[2][1];

	cp = Cross_prod(x2,y2,x0,y0,x1,y1);

	if (fabs(cp) < lin_cell_tol) 
	    return NO;                  /* degenerated triangle */
	else if (cp > 0.0)		/* counterclockwise triangle */
	{
	    cp = Cross_prod(x1,y1,coords[0],coords[1],x0,y0);
	    if (cp < 0.0)
	    	return NO;
	    cp = Cross_prod(x2,y2,coords[0],coords[1],x1,y1);
	    if (cp < 0.0)
	    	return NO;
	    cp = Cross_prod(x0,y0,coords[0],coords[1],x2,y2);
	    if (cp < 0.0)
	    	return NO;
	}
	else				      /* clockwise triangle */
	{
	    cp = Cross_prod(x1,y1,coords[0],coords[1],x0,y0);
	    if (cp > 0.0)
	    	return NO;
	    cp = Cross_prod(x2,y2,coords[0],coords[1],x1,y1);
	    if (cp > 0.0)
	    	return NO;
	    cp = Cross_prod(x0,y0,coords[0],coords[1],x2,y2);
	    if (cp > 0.0)
	    	return NO;
	}
	return YES;
}	/*end test_point_in_tri*/

static	float lin_intrp_var(
	float *crds,
	INTRP_CELL *blk_cell)
{
	float **p = blk_cell->p_lin;
	float *var = blk_cell->var_lin;
	float ans;
	int dim = blk_cell->dim;
	float f[MAXD+1];

	switch (dim)
	{
	case 2:
	{
	    float x0, y0, x1, y1, x2, y2;
	    float xx, yy;
	    float den;

	    x0 = p[0][0];    y0 = p[0][1];
	    x1 = p[1][0] - x0;    y1 = p[1][1] - y0;
	    x2 = p[2][0] - x0;    y2 = p[2][1] - y0;
	    xx = crds[0] - x0;        yy = crds[1] - y0;
	    den = x1*y2 - y1*x2;
	    f[1] = (xx*y2 - yy*x2) / den;
	    f[2] = (x1*yy - y1*xx) / den;
	    f[0] = 1.0 - f[1] - f[2];
	    f[0] = max(0.0,f[0]);	f[0] = min(1.0,f[0]);
	    f[1] = max(0.0,f[1]);	f[1] = min(1.0,f[1]);
	    f[2] = max(0.0,f[2]);	f[2] = min(1.0,f[2]);
	    ans = f[0]*var[0] + f[1]*var[1] + f[2]*var[2];
	}
	break;
	case 3:
	{
	    float v10, v11, v12;
	    float v20, v21, v22;
	    float v30, v31, v32;
	    float q0, q1, q2;
	    float *p0, *p1, *p2, *p3;
	    float den;

	    p0 = p[0];    p2 = p[2];
	    p1 = p[1];    p3 = p[3];
	    q0 = crds[0] - p0[0]; q1 = crds[1] - p0[1]; q2 = crds[2] - p0[2];
	    v10 = p1[0] - p0[0]; v11 = p1[1] - p0[1]; v12 = p1[2] - p0[2];
	    v20 = p2[0] - p0[0]; v21 = p2[1] - p0[1]; v22 = p2[2] - p0[2];
	    v30 = p3[0] - p0[0]; v31 = p3[1] - p0[1]; v32 = p3[2] - p0[2];
	    den = QDet3d(v1,v2,v3);
	    if (fabs(den) < MACH_EPS)
	    {
	        f[0] = 0.0;
	        f[1] = 0.0;
	        f[2] = 0.0;
	        f[3] = 1.0;
	    }

	    f[1] = QDet3d(q,v2,v3)/den;
	    f[2] = QDet3d(v1,q,v3)/den;
	    f[3] = QDet3d(v1,v2,q)/den;
	    f[0] = 1.0 - f[1] - f[2] - f[3];
	    f[0] = max(0.0,f[0]);	f[0] = min(1.0,f[0]);
	    f[1] = max(0.0,f[1]);	f[1] = min(1.0,f[1]);
	    f[2] = max(0.0,f[2]);	f[2] = min(1.0,f[2]);
	    f[3] = max(0.0,f[3]);	f[3] = min(1.0,f[3]);
	    ans = f[0]*var[0] + f[1]*var[1] + f[2]*var[2] + f[3]*var[3];
	}
	}
	return ans;
}	/* end lin_intrp_var */

LOCAL float sine_curve(
        POINTER func_params,
        float *coords)
{

        float dist, theta;

	float y_intfc = 0.5 + 0.2*sin(6.0*PI*coords[0]);

        dist =   y_intfc - coords[1];
        return dist;
}       /* end sine_curve */

static	void print_front_and_interior_state(
	const char *out_name,
	Front *front,
	float *solute)
{
	INTERFACE *intfc = front->interf;
	INTERFACE *grid_intfc;
	RECT_GRID *grid;
	int i,j,k,dim,*gmax;
	float *L,*h,coords[MAXD];
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
	    fprintf(outfile,"%20.16f %20.16f\n",(float)*sl,(float)*sr);
	}

	/* Initialize states in the interior regions */

	if (front->grid_intfc == NULL)
	{
	    grid_intfc_made = YES;
	    grid_intfc = front->grid_intfc =
		make_grid_intfc(intfc,EXPANDED_DUAL_GRID,NULL);
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
	    free_grid_intfc(grid_intfc);
}	/* print_front_and_interior_state */

static	void read_print_front_and_interior_state(
	const char *in_name,
	Front *front,
	float **psolute_old,
	float **psolute_new)
{
	INTERFACE *intfc = front->interf;
	INTERFACE *grid_intfc;
	RECT_GRID *grid;
	int i,j,k,dim,*gmax,num_nodes;
	float *L,*h,coords[MAXD];
	Table *T;
	STATE *sl,*sr;
	POINT *p;
	HYPER_SURF *hs;
	HYPER_SURF_ELEMENT *hse;
	FILE *infile;
	float x,*solute,*solute_new;

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

	grid_intfc = front->grid_intfc =
		make_grid_intfc(intfc,EXPANDED_DUAL_GRID,NULL);
	grid = &topological_grid(grid_intfc);
	dim = grid->dim;
	gmax = grid->gmax;
	L = grid->L;
	h = grid->h;
	num_nodes = 1;
	for (i = 0; i < dim; ++i)
	    num_nodes *= (gmax[i] + 1);

	solute = (float*)malloc(num_nodes*sizeof(float));
	solute_new = (float*)malloc(num_nodes*sizeof(float));
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


static float sphere_func(
        POINTER func_params,
        float *coords)
{
        TEST_SPHERE_PARAMS *s_params = (TEST_SPHERE_PARAMS*)func_params;
	float x0,y0,z0,R;
	float distance;

        x0 = s_params->center[0];
        y0 = s_params->center[1];
        z0 = s_params->center[2];
	R = s_params->radius;

	distance = sqrt(sqr(coords[0] - x0) + sqr(coords[1] - y0) +
			sqr(coords[2] - z0)) - R;

        return distance;

}       /* end sphere_func */



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
