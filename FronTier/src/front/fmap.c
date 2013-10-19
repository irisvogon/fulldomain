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
*				fmap.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#include <front/fdecs.h>		/* includes int.h, table.h */

LOCAL 	bool build_linear_element(INTRP_CELL*,float*);
LOCAL 	void collect_cell_ptst(INTRP_CELL*,int*,COMPONENT,Front*,float*,
			float (*func)(Locstate));
LOCAL	bool test_point_in_seg(float*,float**);
LOCAL	bool test_point_in_tri(float*,float**);
LOCAL 	bool test_point_in_tetra(float*,float**);

EXPORT	int FrontAdvance(
        float    dt,
        float    *dt_frac,
        Front    *front,
        Front    **newfront,
        POINTER  wave)
{
	int status,count;
	float start_dt = dt;

	*dt_frac = 1.0;
	front->dt_frac = dt_frac;
        status = advance_front(front->dt,dt_frac,front,newfront,wave);

        count = 0;
        while (status == MODIFY_TIME_STEP || status == REPEAT_TIME_STEP)
        {
	    if (status == MODIFY_TIME_STEP)
            	front->dt = (*dt_frac)*start_dt;
            status = advance_front(front->dt,dt_frac,front,newfront,wave);
            count++;
            if (count > 10) 
	    {
		screen("ERROR: in FrontAdvance() modified step 10 times\n");
	    	clean_up(ERROR);
	    }
        }
}	/* end FrontAdvance */

EXPORT	float FrontHypTimeStep(
	Front *front)
{
	float fcrds[MAXD];
	float max_dt;

	/* f_max_front_time_step */
	max_dt = (*front->max_front_time_step)(front,fcrds);
#if defined(__MPI__)
	pp_global_min(&max_dt,1);
#endif /* defined(__MPI__) */
	return max_dt;
}	/* end FrontTimeStep */

EXPORT	float FrontOutputTimeControl(
	Front *front,
	bool *is_movie_time,
	bool *is_print_time,
	bool *time_limit_reached,
	int *im,		/* Printing number */
	int *ip)		/* Movie Frame number */
{
	float time = front->time;
	float dt = front->dt;
	float new_dt;
	float dt1,dt2,dt3;

	*is_movie_time = *is_print_time = *time_limit_reached = NO;

	dt1 = (*im)*front->movie_frame_interval - time;
	dt2 = (*ip)*front->print_time_interval - time;
	dt3 = front->max_time - time;
	new_dt = min3(dt1,dt2,dt3);

	if (front->step >= front->max_step)
            *time_limit_reached = YES;

	if (new_dt > dt) 
	    return dt;

	if (dt1 == new_dt)
	{
	    *is_movie_time = YES;
	    (*im)++;
	}
        if (dt2 == new_dt)
        {
            *is_print_time = YES;
	    (*ip)++;
        }
        if (dt3 == new_dt)
        {
            *time_limit_reached = YES;
        }
	return new_dt;
}	/* end FrontOutputTimeControl */

EXPORT  void FrontRedistMesh(
        Front *fr)
{
        INTERFACE *intfc = fr->interf;
        CURVE **c;
        SURFACE **s;
        int dim = intfc->dim;
        bool force_redistribute = YES;

        switch (dim)
        {
        case 2:
            Curve_redistribute(fr,&force_redistribute);
            break;
        case 3:
	    reset_normal_on_intfc(fr->interf);
            Surface_redistribute(fr,&force_redistribute);
        }
}       /* end FrontRedistMesh */

EXPORT  void FrontSetSpacing(
        Front *fr,
        float spacing)
{
        Front_spacing(fr,GENERAL_WAVE) = spacing;
}       /* end FrontSetSpacing */

EXPORT	void FrontSetTriParams(
	Front *fr,
	float max_tri_area_fac,
	float min_tri_area_fac,
	float min_angle_at_vertex,	/* In unit of degree */
	float max_scaled_tri_side)
{
	float sslen,std_area,hmin,a;

	hmin = fr->rect_grid->h[0];
	if (hmin > fr->rect_grid->h[1])
	    hmin = fr->rect_grid->h[1];
	if (hmin > fr->rect_grid->h[2])
	    hmin = fr->rect_grid->h[2];

	sslen = Front_spacing(fr,GENERAL_WAVE);
	std_area = 0.25*sqrt(3.0)*sqr(sslen)*sqr(hmin);
	Max_tri_sqr_area(fr,GENERAL_WAVE) = sqr(max_tri_area_fac*std_area);
	Min_tri_sqr_area(fr,GENERAL_WAVE) = sqr(min_tri_area_fac*std_area);
	Max_bond_len(fr,GENERAL_WAVE) = sslen*sqrt(max_tri_area_fac)*hmin;
	Min_bond_len(fr,GENERAL_WAVE) = sslen*sqrt(min_tri_area_fac)*hmin;
	a = radians(min_angle_at_vertex);
	Aspect_ratio_tolerance(fr,GENERAL_WAVE) = 0.5*sin(a)/(2.0 - cos(a));
	Max_scaled_tri_side_sqr_length(fr) = sqr(max_scaled_tri_side);
}	/* end FrontSetTriParams */

EXPORT 	void FrontResetTriParams(
	Front *fr)
{
	float max_tri_area_fac = 2.0;
	float min_tri_area_fac = 0.5;
	float min_angle_at_vertex = radians(15.0);
	float max_scaled_tri_side = 1.3;

	FrontSetTriParams(fr,max_tri_area_fac,
			     min_tri_area_fac,
			     min_angle_at_vertex,
			     max_scaled_tri_side);
}	/* end FrontResetTriParams */


EXPORT	void	FrontInitStandardIO(
	int		argc,
	char		**argv,
	F_BASIC_DATA 	*f_basic)
{
	char *in_name      = f_basic->in_name;
	char *out_name     = f_basic->out_name;
	char *restart_name = f_basic->restart_name;
	int  *subdomains   = f_basic->subdomains;
	f_basic->ReadFromInput = NO;
	f_basic->RestartRun = NO;
#if defined(__MPI__)
        pp_init(&argc,&argv);
#endif /* defined(__MPI__) */
	argc--;
	argv++;
	strcpy(out_name,"intfc");
	while (argc >= 1)
	{
	    if (argv[0][0] != '-')
	    {
		printf("Usage: example -o output\n");
		exit(1);
	    }
	    switch(argv[0][1]) {
	    case 'i':
	    case 'I':
	    	f_basic->ReadFromInput = YES;
	    	zero_scalar(in_name,200);
                strcpy(in_name,argv[1]);
                argc -= 2;
		argv += 2;
		break;
	    case 'r':
	    case 'R':
	    	f_basic->RestartRun = YES;
	    	zero_scalar(restart_name,200);
                strcpy(restart_name,argv[1]);
                argc -= 2;
		argv += 2;
		break;
	    case 't':
	    case 'T':
	    	f_basic->RestartStep = atoi(argv[1]);
                argc -= 2;
                argv += 2;
                break;
	    case 'd':
	    case 'D':
	    	f_basic->dim = atoi(argv[1]);
                argc -= 2;
                argv += 2;
                break;
	    case 'o':
	    case 'O':
		zero_scalar(out_name,200);
		strcpy(out_name,argv[1]);
#if defined(__MPI__)
                sprintf(out_name,"%s.%d",out_name,pp_mynode());
#endif /* defined(__MPI__) */
		freopen(out_name,"w",stdout);
		zero_scalar(out_name,200);
		strcpy(out_name,argv[1]);
		argc -= 2;
		argv += 2;
		break;
#if defined(__MPI__)
            case 'p':
            case 'P':
	    {
                int i,total_num_proc = 1;
                for (i = 0; i < MAXD; ++i)
                {
                    if (argc < 2 || argv[1][0] == '-') break;
                    argc -= 1;
                    argv += 1;
                    subdomains[i] = atoi(argv[0]);
                    total_num_proc *= subdomains[i];
                }
                argc -= 1;
                argv += 1;
                if (total_num_proc != pp_numnodes())
                {
                    printf("total number of processors for the partition %d "
                           "does not equal to requested np %d\n",
                           total_num_proc,pp_numnodes());
                    clean_up(ERROR);
                }
	    }
#endif /* defined(__MPI__) */
	    }
	}
}	/* end FrontInitStatndardIO */

EXPORT	void FrontMovieFrame(
	Front *front,
	char *out_name,
	bool print_in_binary)
{
	show_front_output(front,out_name,print_in_binary);
}	/* end FrontMovieFrame */

EXPORT	void FrontPrintOut(
	Front *front,
	char *out_name)
{
	print_front_output(front,out_name);
}	/* end FrontPrintOut */

EXPORT	void FrontFreeAll(
	Front *front)
{
	free_front(front);
}	/* end FrontFreeAll */

EXPORT	void FrontMakeGridIntfc(
	Front *front)
{
	front->grid_intfc = make_grid_intfc(front->interf,
			EXPANDED_DUAL_GRID,NULL);
}	/* end FrontMakeGridIntfc */

EXPORT	void FrontFreeGridIntfc(
	Front *front)
{
	free_grid_intfc(front->grid_intfc);
}	/* end FrontFreeGridIntfc */

EXPORT	void FrontScatIntfcBuffer(
	Front *front)
{
	scatter_front(front);
}	/* end FrontScatIntfcBuffer */

EXPORT	void FrontScatExpDualGridArray(
	float *grid_array,
	Front *front)
{
	scatter_top_grid_float_array(grid_array,front);
}	/* end FrontScatExpDualGridArray */

EXPORT	bool FrontStateAtGridCrossing(
	Front *front,
	int *icoords,
	GRID_DIRECTION dir,
	COMPONENT comp,
	float (*state_func)(Locstate),
	float *ans,
	float *crx_coords)
{
	Locstate state;
	HYPER_SURF *hs;
	INTERFACE *grid_intfc = front->grid_intfc;
	static CRXING *crxs[MAX_NUM_CRX];
	int i,nc,dim = grid_intfc->dim;

	nc = GridSegCrossing(crxs,icoords,dir,grid_intfc);
	if (nc == 0) return NO;
	hs = crxs[0]->hs;
	if (comp == negative_component(hs))
	{
	    state = left_state(crxs[0]->pt);
	    for (i = 0; i < dim; ++i)
	    	crx_coords[i] = Coords(crxs[0]->pt)[i];
	}
	else if (comp == positive_component(hs))
	{
	    state = right_state(crxs[0]->pt);
	    for (i = 0; i < dim; ++i)
	    	crx_coords[i] = Coords(crxs[0]->pt)[i];
	}
	else
	{
	    screen("ERROR: In FrontStateAtGridCrossing(),"
	           "component does not match\n");
	    return NO;
	}
	*ans = (*state_func)(state);
	return YES;
}	/* end FrontStateAtGridCrossing */

#define         MAX_NUM_VERTEX_IN_CELL          20
LOCAL float lin_cell_tol;

EXPORT	bool FrontIntrpStateAtPoint(
	Front *front,
	COMPONENT comp,
	float *coords,
	float *grid_array,
	float (*get_state)(Locstate),
	float *ans)
{
	int icoords[MAXD];
	INTERFACE *grid_intfc = front->grid_intfc;
	static INTRP_CELL *blk_cell;
	RECT_GRID *gr = &topological_grid(grid_intfc);
	int i,dim = gr->dim;

	if (blk_cell == NULL)
	{
	    scalar(&blk_cell,sizeof(INTRP_CELL));
	    uni_array(&blk_cell->var,MAX_NUM_VERTEX_IN_CELL,sizeof(float));
	    bi_array(&blk_cell->coords,MAX_NUM_VERTEX_IN_CELL,MAXD,sizeof(float));
	    bi_array(&blk_cell->p_lin,MAXD+1,MAXD,sizeof(float));
	    uni_array(&blk_cell->var_lin,MAXD+1,sizeof(float));
	    lin_cell_tol = 1.0;
	    for (i = 0; i < dim; ++i)
	    	lin_cell_tol *= 0.00001*gr->h[i];
	}

	if (!rect_in_which(coords,icoords,gr))
	{
	    return NO;
	}
	collect_cell_ptst(blk_cell,icoords,comp,front,grid_array,get_state);
	if (build_linear_element(blk_cell,coords))
	{
	    *ans = FrontLinIntrp(coords,blk_cell,NO);
	    return YES;
	}
	else
	    return NO;
}	/* end FrontIntrpStateAtPoint */

LOCAL bool build_linear_element(
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
	case 1:
	    for (i = 0; i < nv; i++)
	    {
	    	p[0] = blk_cell->coords[i];
	    	var[0] =  blk_cell->var[i];
	    	for (j = i+1; j < nv; j++)
		{
	    	    p[1] = blk_cell->coords[j];
	    	    var[1] = blk_cell->var[j];
		    if (test_point_in_seg(coords,p) == YES)
	    		return FUNCTION_SUCCEEDED;
		}
	    }
	    break;
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

LOCAL void collect_cell_ptst(
	INTRP_CELL *blk_cell,
	int *icoords,
	COMPONENT comp,
	Front *front,
	float *grid_array,
	float (*get_state)(Locstate))
{
	INTERFACE *grid_intfc = front->grid_intfc;
	Table *T = table_of_interface(grid_intfc);
	RECT_GRID *gr = &topological_grid(grid_intfc);
	int dim = gr->dim;
	int *gmax = gr->gmax;
	float *L = gr->L;
	float *h = gr->h;
	COMPONENT *gr_comp = T->components;
	static COMPONENT cell_comp1d[2];
	static COMPONENT cell_comp2d[2][2];
	static COMPONENT cell_comp3d[2][2][2];
	int i,j,k,index,nv,nc;
	CRXING *crx,*crxs[MAX_NUM_CRX];
	GRID_DIRECTION dir;
	int ic[MAXD];
	bool fr_crx_grid_seg;
	float state_at_crx;
	float crx_coords[MAXD];
	
	blk_cell->is_bilinear = YES;
	blk_cell->dim = dim;
	nv = 0;
	switch (dim)
	{
	case 1:
	    for (i = 0; i < 2; ++i)
	    {
	    	ic[0] = icoords[0] + i;
	    	index = d_index1d(ic[0],gmax);
	    	cell_comp1d[i] = gr_comp[index];
	    	if (gr_comp[index] == comp)
	    	{
	    	    blk_cell->coords[nv][0] = L[0] + ic[0]*h[0];
		    blk_cell->var[nv] = grid_array[index];
		    nv++;
	    	}
	    	else
	    	    blk_cell->is_bilinear = NO;
	    }
	    break;
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
		    blk_cell->var[nv] = grid_array[index];
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
		    blk_cell->var[nv] = grid_array[index];
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
	case 1:
	    for (i = 0; i < 2; ++i)
	    {
	    	ic[0] = icoords[0] + i;
	    	if (cell_comp1d[i] == comp)
		{
	    	    if (cell_comp1d[(i+1)%2] != comp)
		    {
		    	dir = (i < (i+1)%2) ? EAST : WEST;
			fr_crx_grid_seg = FrontStateAtGridCrossing(front,ic,dir,
				comp,get_state,&state_at_crx,crx_coords);
		    	if (!fr_crx_grid_seg)
		    	{
		    	    screen("ERROR: no crxing between (%d) and (%d)\n",
			    		icoords[0]+i,icoords[0]+(i+1)%2);
		    	}
		    	blk_cell->var[nv] = state_at_crx;
		    	blk_cell->coords[nv][0] = crx_coords[0];
		    	nv++;
		    }
		}
	    }
	    break;
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
			fr_crx_grid_seg = FrontStateAtGridCrossing(front,ic,dir,
				comp,get_state,&state_at_crx,crx_coords);
		    	if (!fr_crx_grid_seg)
		    	{
		    	    screen("ERROR: no crxing between (%d %d) "
			           "and (%d %d)\n",icoords[0]+i,icoords[1]+j,
					icoords[0]+(i+1)%2,icoords[1]+j);
		    	}
		    	blk_cell->var[nv] = state_at_crx;
		    	blk_cell->coords[nv][0] = crx_coords[0];
		    	blk_cell->coords[nv][1] = crx_coords[1];
		    	nv++;
	    	    }
	    	    if (cell_comp2d[i][(j+1)%2] != comp)
	    	    {
		    	dir = (j < (j+1)%2) ? NORTH : SOUTH;
			fr_crx_grid_seg = FrontStateAtGridCrossing(front,ic,dir,
				comp,get_state,&state_at_crx,crx_coords);
		    	if (!fr_crx_grid_seg)
		    	{
		    	    screen("ERROR: no crxing between (%d %d) "
			           "and (%d %d)\n",icoords[0]+i,icoords[1]+j,
					icoords[0]+i,icoords[1]+(j+1)%2);
		    	}
		    	blk_cell->var[nv] = state_at_crx;
		    	blk_cell->coords[nv][0] = crx_coords[0];
		    	blk_cell->coords[nv][1] = crx_coords[1];
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
			fr_crx_grid_seg = FrontStateAtGridCrossing(front,ic,dir,
				comp,get_state,&state_at_crx,crx_coords);
		    	if (!fr_crx_grid_seg)
		    	{
		    	    screen("ERROR: no crxing between (%d %d %d) "
			           "and (%d %d %d)\n",icoords[0]+i,icoords[1]+j,
				   icoords[2]+k,icoords[0]+(i+1)%2,icoords[1]+j,
				   icoords[2]+k);
		    	}
		    	blk_cell->var[nv] = state_at_crx;
		    	blk_cell->coords[nv][0] = crx_coords[0];
		    	blk_cell->coords[nv][1] = crx_coords[1];
		    	blk_cell->coords[nv][2] = crx_coords[2];
		    	nv++;
	    	    }
	    	    if (cell_comp3d[i][(j+1)%2][k] != comp)
	    	    {
		    	dir = (j < (j+1)%2) ? NORTH : SOUTH;
			fr_crx_grid_seg = FrontStateAtGridCrossing(front,ic,dir,
				comp,get_state,&state_at_crx,crx_coords);
		    	if (!fr_crx_grid_seg)
		    	{
		    	    screen("ERROR: no crxing between (%d %d %d) "
			           "and (%d %d %d)\n",icoords[0]+i,icoords[1]+j,
				   icoords[2]+k,icoords[0]+i,icoords[1]+(j+1)%2,
				   icoords[2]+k);
		    	}
		    	blk_cell->var[nv] = state_at_crx;
		    	blk_cell->coords[nv][0] = crx_coords[0];
		    	blk_cell->coords[nv][1] = crx_coords[1];
		    	blk_cell->coords[nv][2] = crx_coords[2];
		    	nv++;
	    	    }
	    	    if (cell_comp3d[i][j][(k+1)%2] != comp)
	    	    {
		    	dir = (k < (k+1)%2) ? UPPER : LOWER;
			fr_crx_grid_seg = FrontStateAtGridCrossing(front,ic,dir,
				comp,get_state,&state_at_crx,crx_coords);
		    	if (!fr_crx_grid_seg)
		    	{
		    	    screen("ERROR: no crxing between (%d %d %d) "
			           "and (%d %d %d)\n",icoords[0]+i,icoords[1]+j,
				   icoords[2]+k,icoords[0]+i,icoords[1]+j,
				   icoords[2]+(k+1)%2);
		    	}
		    	blk_cell->var[nv] = state_at_crx;
		    	blk_cell->coords[nv][0] = crx_coords[0];
		    	blk_cell->coords[nv][1] = crx_coords[1];
		    	blk_cell->coords[nv][2] = crx_coords[2];
		    	nv++;
	    	    }
	    	}
	    }
	    break;
	}
	blk_cell->nv = nv;
}	/* end collect_cell_ptst */


LOCAL 	bool test_point_in_tetra(
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


LOCAL	bool	test_point_in_seg(
	float	*coords,
	float	**p)
{
	if ((p[0][0]-lin_cell_tol < coords[0] &&
	     p[1][0]+lin_cell_tol > coords[0]) ||
	    (p[1][0]-lin_cell_tol < coords[0] &&
	     p[0][0]+lin_cell_tol > coords[0])) 
	    return YES;
	else
	    return NO;
}
	/* end test_point_in_seg */

LOCAL	bool	test_point_in_tri(
	float	*coords,
	float	**p)
{
	float	x0, y0, x1, y1, x2, y2, cp;
	float	v1[2],v2[2];

	x0 = p[0][0];	y0 = p[0][1];
	x1 = p[1][0];	y1 = p[1][1];
	x2 = p[2][0];	y2 = p[2][1];

	v1[0] = x2 - x1;	v1[1] = y2 - y1;
	v2[0] = x0 - x1;	v2[1] = y0 - y1;
	Cross2d(v1,v2,cp);

	if (fabs(cp) < lin_cell_tol) 
	    return NO;                  /* degenerated triangle */
	else if (cp > 0.0)		/* counterclockwise triangle */
	{
	    v1[0] = x1 - x0;		v1[1] = y1 - y0;
	    v2[0] = coords[0] - x0;	v2[1] = coords[1] - y0;
	    Cross2d(v1,v2,cp);
	    if (cp < 0.0)
	    	return NO;
	    v1[0] = x2 - x1;		v1[1] = y2 - y1;
	    v2[0] = coords[0] - x1;	v2[1] = coords[1] - y1;
	    Cross2d(v1,v2,cp);
	    if (cp < 0.0)
	    	return NO;
	    v1[0] = x0 - x2;		v1[1] = y0 - y2;
	    v2[0] = coords[0] - x2;	v2[1] = coords[1] - y2;
	    Cross2d(v1,v2,cp);
	    if (cp < 0.0)
	    	return NO;
	}
	else				      /* clockwise triangle */
	{
	    v1[0] = x1 - x0;		v1[1] = y1 - y0;
	    v2[0] = coords[0] - x0;	v2[1] = coords[1] - y0;
	    Cross2d(v1,v2,cp);
	    if (cp > 0.0)
	    	return NO;
	    v1[0] = x2 - x1;		v1[1] = y2 - y1;
	    v2[0] = coords[0] - x1;	v2[1] = coords[1] - y1;
	    Cross2d(v1,v2,cp);
	    if (cp > 0.0)
	    	return NO;
	    v1[0] = x0 - x2;		v1[1] = y0 - y2;
	    v2[0] = coords[0] - x2;	v2[1] = coords[1] - y2;
	    Cross2d(v1,v2,cp);
	    if (cp > 0.0)
	    	return NO;
	}
	return YES;
}	/*end test_point_in_tri*/

EXPORT	float FrontLinIntrp(
	float *crds,
	INTRP_CELL *blk_cell,
	bool reuse_coeffs)
{
	float **p = blk_cell->p_lin;
	float *var = blk_cell->var_lin;
	float ans;
	int dim = blk_cell->dim;
	static float f[MAXD+1];

	switch (dim)
	{
	case 1:
	{
	    float den = fabs(p[0][0] - p[1][0]);
	    f[0] = fabs(crds[0] - p[1][0])/den;
	    f[1] = fabs(crds[0] - p[0][0])/den;
	    if (f[0] < 0.0)
	    {
	    	f[0] = 0.0;	f[1] = 1.0;
	    }
	    else if (f[0] > 1.0)
	    {
	    	f[0] = 1.0;	f[1] = 0.0;
	    }
	    ans = f[0]*var[0] + f[1]*var[1];
	}
	break;
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
}	/* end FrontLinIntrp */

EXPORT  float FrontBilinIntrp(
        float *crds,
        INTRP_CELL *blk_cell,
        bool reuse_coeffs)
{
	static float f[3][2];
	int i,j,k,il,iu,index;
	int dim = blk_cell->dim;
	float ans;
	float l[MAXD],u[MAXD];

	if (!reuse_coeffs)
	{
	    il = 0;
	    iu = (dim == 2) ? 3 : 7;
	    for (i = 0; i < dim; ++i)
	    {
	    	l[i] = blk_cell->coords[il][i];
		u[i] = blk_cell->coords[iu][i];
	    	f[i][0] = (u[i] - crds[i])/(u[i] - l[i]);
		f[i][0] = min(f[i][0],1.0);
		f[i][0] = max(f[i][0],0.0);
	    	f[i][1] = 1.0 - f[i][0];
	    }
	}
	index = 0;
	ans = 0.0;
	switch (dim)
	{
	case 2:
	    for (i = 0; i < 2; ++i)
	    for (j = 0; j < 2; ++j)
	    {
	    	ans += blk_cell->var[index++]*f[0][i]*f[1][j];
	    }
	    break;
	case 3:
	    for (i = 0; i < 2; ++i)
	    for (j = 0; j < 2; ++j)
	    for (k = 0; k < 2; ++k)
	    {
	    	ans += blk_cell->var[index++]*f[0][i]*f[1][j]*f[2][k];
	    }
	}
	return ans;

}	/* end FrontBilinIntrp */

EXPORT	bool FrontCpuAdaptSubdomain(
	Front *front,
	float cpu_time,
	int *lexpand,
	int *uexpand)
{
	return cpu_adapt_front(front,cpu_time,lexpand,uexpand);
}	/* end FrontCpuAdaptSubdomain */
