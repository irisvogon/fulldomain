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
*				hnpt.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains n-point stencil code for finite difference solution
*	of the hyperbolic equations within a single component of a tracked
*	interface problem.
*	Point sources in the region are found and passed appropriately.
*/

#define DEBUG_STRING	"hnpt"

#include <hyp/hdecs.h>

enum _CRX_SIDE {
	PREV = -1,
	NEXT =  1
};
typedef enum _CRX_SIDE CRX_SIDE;


LOCAL	int	gmax[MAXD], vsten_rad, is_vec_method;
LOCAL	int	dim;
LOCAL   int	lbuf[MAXD], ubuf[MAXD];
LOCAL	float	h[MAXD];
LOCAL	float	L[MAXD], VL[MAXD], VU[MAXD];

	/* Redefine rect grid macros for local efficiency */
	/* these will have to be modified to support variable mesh grids */
#undef cell_center
#define cell_center(indx,i,gr)	(L[i] + ((indx) + 0.5)*h[i])

	/* LOCAL Function Declarations */
LOCAL	COMPONENT	set_stencil_comps_and_crosses(int,Stencil*);
LOCAL	bool	interior_sten_reg(int,int,Stencil*);
LOCAL	void	copy_states(int,CRX_SIDE,int,int,Stencil*,CRXING*);
LOCAL	void	find_outer_states(int,int,Stencil*,CRX_SIDE);
LOCAL	void	set_dirichlet_bdry_sten_states(int,HYPER_SURF*,Stencil*,
                                               CRX_SIDE);
LOCAL	void	set_neumann_bdry_sten_states(int,int,int,COMPONENT,CRXING*,
					     Stencil*,CRX_SIDE);
LOCAL	void	update_reg_grid_state(float,float,int,int*,int,
				      Stencil*,POINT*,COMPONENT);
LOCAL   bool Debug;

#if defined(TWOD) || defined(THREED)
LOCAL	void	set_static_coords(int**,POINT*,int,int,RECT_GRID*);
#endif /* defined(TWOD) || defined(THREED) */

LOCAL	void	interp_states(int,CRX_SIDE,int,int,int,Stencil*,CRXING*);
	void	interp_dirichlet_states(int,int,int,int,Stencil*,CRXING*);

EXPORT	void set_hyp_npt_globals(
	Wave		*wave)
{
	RECT_GRID	*r_grid;
	int		i;

	r_grid = wave->rect_grid;
	dim = r_grid->dim;
	for (i = 0; i < dim; ++i)
	{
	    gmax[i] = r_grid->gmax[i];
	    h[i] = r_grid->h[i];
	    L[i] = r_grid->L[i];
	    VL[i] = r_grid->VL[i];
	    VU[i] = r_grid->VU[i];
	    lbuf[i] = r_grid->lbuf[i];
	    ubuf[i] = r_grid->ubuf[i];
	}
	vsten_rad = vsten_radius(wave);
	is_vec_method = is_vector_method(wave);
}		/*end set_hyp_npt_globals*/

/*
*	The following globals are shared between hyp_npt() and
*	the subroutine update_reg_grid_state().  They are used
*	to pass quantities that do not change inside the inner loop
*	of the update step.  The are passed as globals for efficiency.
*/

LOCAL	COMPONENT	ext_comp;
LOCAL	float		hdir, dir[MAXD];
LOCAL	GRID_DIRECTION	prev_side, next_side;
LOCAL	INTERFACE	*intfc = NULL;
LOCAL	int		endpt;


/*
*			hyp_npt():
*
*	Single direction sweep for a n-point scheme.
*/

EXPORT void hyp_npt(
	int		swp_num,
	int		*iperm,	    /* sweep based permutation of coord dir */
	float		ds,
	float		dt,
	Wave		*wave,
	Wave		*newwave,
	Front		*fr,
	Front		*newfr,	    /* newfr needed if hlw is to support */
		    		    /*	changing topologies	     */
	COMPONENT	max_comp)
{
	POINT		Basep;
	int 		i, idirs[MAXD];
	int		imin[3], imax[3];
	static	Stencil	*sten = NULL;
	static	int	num_pts;
#if defined(TWOD) || defined(THREED)
	static	int	**icoords = NULL;
	RECT_GRID	*gr = wave->rect_grid;
#endif /* defined(TWOD) || defined(THREED) */

	for (i = 0; i < dim; ++i)
	{
	    idirs[i] = iperm[(i+swp_num)%dim];
	    dir[i] = 0.0;
	}
	dir[idirs[0]] = 1.0;

	intfc = fr->interf;
	ext_comp = exterior_component(intfc);

	if (sten == NULL)
	{
	    /* This assumes that the stencil size is constant. */
	    num_pts = wave->npts_sten;
	    endpt = stencil_radius(wave);

	    sten = alloc_stencil(num_pts,fr);
#if defined(TWOD) || defined(THREED)
	    icoords = sten->icoords;
#endif /* defined(TWOD) || defined(THREED) */
	}

	hdir = ds;		/* grid spacing in sweep direction */
	switch (idirs[0])
	{
	case 0:	/*X SWEEP */
	    prev_side = WEST;
	    next_side = EAST;
	    if (debugging("hyp_npt"))
	    	(void) printf("prev_side = WEST, next_side = EAST\n");
	    break;
	case 1: /*Y SWEEP */
	    prev_side = SOUTH;
	    next_side = NORTH;
	    if (debugging("hyp_npt"))
	    	(void) printf("prev_side = SOUTH, next_side = NORTH\n");
	    break;
	case 2: /* Z SWEEP */
	    prev_side = LOWER;
	    next_side = UPPER;
	    if (debugging("hyp_npt"))
	    	(void) printf("prev_side = LOWER, next_side = UPPER\n");
	    break;
	}

	sten->fr   = fr;	sten->newfr   = newfr;
	sten->wave = wave;	sten->newwave = newwave;

	set_sweep_limits(wave,swp_num,idirs,imin,imax);
	//#bjet2
	set_limits_for_open_bdry(wave,fr,swp_num,idirs,imin,imax);

	switch (dim)
	{
#if defined(ONED)
	case 1:
	{
	    int	i0;
	    int	i0max;
	    int	i0min;

	    i0min = imin[0];	i0max = imax[0]; 
	    sten->reg_stencil = NO;
	    for (i0 = i0min; i0 < i0max; ++i0) 
	    {
	    	update_reg_grid_state(ds,dt,swp_num,iperm,
				      i0,sten,&Basep,max_comp);
	    }
	    break;
	}
#endif /* defined(ONED) */
#if defined(TWOD)
	case 2:
	{
	    int	i0, i1;
	    int	i0max, i1max;
	    int	i0min, i1min;

	    i0min = imin[0];	i0max = imax[0]; 
	    i1min = imin[1];	i1max = imax[1]; 
	    for (i1 = i1min; i1 < i1max; ++i1)
	    {
	    	set_static_coords(icoords,&Basep,idirs[1],i1,gr);
	    	sten->reg_stencil = NO;
	    	for (i0 = i0min; i0 < i0max; ++i0) 
	    	{
	    	    update_reg_grid_state(ds,dt,swp_num,iperm,
	    	    		          i0,sten,&Basep,max_comp);
	    	}
	    }
	    break;
	}
#endif /* defined(TWOD) */
#if defined(THREED)
	case 3:
	{
	    int	i0, i1, i2;
	    int	i0max, i1max, i2max;
	    int	i0min, i1min, i2min;

	    i0min = imin[0];	i0max = imax[0]; 
	    i1min = imin[1];	i1max = imax[1]; 
	    i2min = imin[2];	i2max = imax[2];

	    for (i2 = i2min; i2 < i2max; ++i2)
	    {
	    	set_static_coords(icoords,&Basep,idirs[2],i2,gr);
	    	for (i1 = i1min; i1 < i1max; ++i1)
	    	{
	    	    set_static_coords(icoords,&Basep,idirs[1],i1,gr);
		    sten->reg_stencil = NO;
		    for (i0 = i0min; i0 < i0max; ++i0) 
		    {
                        for(i = 0; i < num_pts; i++)
                          sten->boundarystate_store[i] = 0;
			update_reg_grid_state(ds,dt,swp_num,iperm,i0,
					      sten,&Basep,max_comp);
		    }
		}
	    }
	    break;
	}
#endif /* defined(THREED) */
	}
	debug_print("hyp_npt","Leaving hyp_npt(), dir = %d\n",idirs[0]);
}		/*end hyp_npt*/


#if defined(TWOD) || defined(THREED)
/*
*			set_static_coords():
*/

/*ARGSUSED*/
LOCAL	void set_static_coords(
	int		**icoords,
	POINT		*basep,
	int		idir,
	int		is,
	RECT_GRID	*gr)
{
	int		i;

	for (i = -endpt; i <= endpt; ++i)
	    icoords[i][idir] = is;
	Coords(basep)[idir] = cell_center(is,idir,gr);
}		/*end set_static_coords*/
#endif /* defined(TWOD) || defined(THREED) */


/*
*			interior_sten_reg():
*
*	If called after a vector solver, decides whether the stencil at a given
*	point was regular for the vector solver.  In other words, decides
*	whether solution needs to be recomputed at a given point.
*/

LOCAL	bool interior_sten_reg(
	int		is,
	int		idir,
	Stencil		*sten)
{
	Wave		*wave = sten->wave;
	Front		*newfr = sten->newfr;
	int 		i;
	int 		**icoords = sten->icoords, icrds[MAXD];
	COMPONENT 	cmp;
	COMPONENT	new_comp = sten->newcomp, *comp = sten->comp;
	CRXING		*cross;

	if (ComponentIsFlowSpecified(new_comp,newfr))
	    return NO;
	for (i = -endpt; i <= endpt; ++i)
	    if ((equivalent_comps(comp[i],new_comp,newfr->interf) == NO) ||
		(sten->nc[i] != 0))
	    return NO;
	if (endpt >= vsten_rad)
	    return YES;

	for (i = 0; i < dim; ++i)
	    icrds[i] = icoords[0][i];
	for (i = endpt+1; i <= vsten_rad; ++i)
	{
	    icrds[idir] = is - i;
	    cmp = Find_rect_comp(idir,icrds,wave);
	    cross = Rect_crossing(icrds,next_side,wave);
	    if ((equivalent_comps(cmp,new_comp,newfr->interf) == NO) ||
		(cross != NULL))
	    	return NO;
		
	    icrds[idir] = is + i;
	    cmp = Find_rect_comp(idir,icrds,wave);
	    cross = Rect_crossing(icrds,prev_side,wave);
	    if ((equivalent_comps(cmp,new_comp,newfr->interf) == NO) ||
		(cross != NULL))
		return NO;
	}
	return YES;
}		/*end interior_sten_reg*/


void    print_states_gas(Locstate);
int	get_st_from_interp(Locstate, int*, Stencil*);
int	get_st_from_ghost(Locstate, int*, Stencil*);
bool	extrp_st_along_normal(Locstate, float*, Stencil*, Front*);

/*
*			update_reg_grid_state():
*
*	Computes the solution at a given point based on at stencil of
*	surrounding points.
*/
void printf_state(int , Locstate);

LOCAL	void update_reg_grid_state(
	float		ds,
	float		dt,
	int		swp_num,
	int		*iperm,
	int		is,		//icoords index in sweep direction
	Stencil		*sten,
	POINT		*basep,
	COMPONENT	max_comp)
{
        Wave		*wave = sten->wave;
	Wave		*newwave = sten->newwave;
	Front		*fr = sten->fr;
	Front		*newfr = sten->newfr;
	COMPONENT	*comp;
	COMPONENT	new_comp;		/* mid comp wrt newfr */
	HYPER_SURF	*hs;
	int		i, j, index, idir = iperm[swp_num];
	int		**icoords = sten->icoords;

	for (i = -endpt; i <= endpt; ++i)
	    icoords[i][idir] = is + i;
	
        sten->prev_reg_stencil = sten->reg_stencil;

		/* Find components and crosses*/

	new_comp = set_stencil_comps_and_crosses(idir,sten);
        //xiaoxue
        new_comp = 3;
	comp = sten->comp;

        if (is_vec_method)
	{
	    if (interior_sten_reg(is,idir,sten))
	    	return;		/*Don't overwrite vector solution*/
	    else
	    	sten->prev_reg_stencil = NO;
	}

//TMP
jump_is_vec_method:

	for (j = 0; j < dim; ++j)
	{
	    if (j == idir)
	    {
	    	Coords(sten->p[0])[idir] = cell_center(is,idir,wave->rect_grid);
		for (i = 1; i <= endpt; ++i)
		{
		    Coords(sten->p[-i])[idir] = Coords(sten->p[0])[idir]-i*hdir;
		    Coords(sten->p[i])[idir]  = Coords(sten->p[0])[idir]+i*hdir;
		}
	    }
	    else
	    {
	    	for (i = -endpt; i <= endpt; ++i)
		    Coords(sten->p[i])[j] = Coords(basep)[j];
	    }
	}
	
        sten->reg_stencil = YES;
	
        if (RegionIsFlowSpecified(Rect_state(icoords[0],newwave),
			          Rect_state(icoords[0],wave),
			          Coords(sten->p[0]),new_comp,comp[0],fr))
	{
	    sten->reg_stencil = NO;
	    return;
	}

	if(the_point(sten->p[0]))
	{
	    printf("#debug interp_npt\n");
	    add_to_debug("interp_npt");
	}

		/* Find mid state */
	if (equivalent_comps(comp[0],new_comp,newfr->interf) == YES)
        {
	    sten->st[0] = Rect_state(icoords[0],wave);
        }
	else 
	{
            detect_and_load_mix_state(iperm[swp_num],sten,0);
	    sten->reg_stencil = NO;
	    return;            
	}
	ft_assign(left_state(sten->p[0]),sten->st[0],fr->sizest);
	ft_assign(right_state(sten->p[0]),sten->st[0],fr->sizest);


	find_outer_states(is,idir,sten,PREV);
	find_outer_states(is,idir,sten,NEXT);
      	

	if(the_point(sten->p[0]))
	{
	    (void) printf("STENCIL: "); print_Stencil(sten);
	    for (i = -endpt; i <= endpt; ++i)
	    {
	    	(void) printf("state[%d]:	",i);
		print_states_gas(sten->st[i]);
	    	//(*fr->print_state)(sten->st[i]);
	    }
	    fflush(NULL);
	    remove_from_debug("interp_npt");
	    //if(sten->fr->step == 1)
		//clean_up(ERROR);
	}
//#endif /* defined(DEBUG_HYP_NPT) */

 		/* Check for point source/sink in icoords */

	index = (wave->num_point_sources) ?
	    is_source_block(wave,newfr->interf,new_comp,icoords[0]) : -1;

 		/* Update state */

	//point_FD
	npt_solver(ds,dt,Rect_state(icoords[0],newwave),dir,swp_num,iperm,
		   &index,sten,wave);
	if(the_point(sten->p[0]))
	{
	    printf("ans=\n");
	    print_states_gas(Rect_state(icoords[0],newwave));
	}

#if defined(DEBUG_HYP_NPT)
	if (deb_hyp_npt)
	{
	    (void) printf("ANSWER: ");
	    (*fr->print_state)(Rect_state(icoords[0],newwave));
	    remove_from_debug("mymuscl");
	}
#endif /* defined(DEBUG_HYP_NPT) */
//        Locstate ppnew = NULL;
//        for (i = -endpt; i <= endpt; ++i)
//        printf("boundarystate[%d] = %d\n", i, sten->boundarystate[i]);

        for (i = -endpt; i <= endpt; ++i)
        {
            if(sten->boundarystate[i]) 
            {
                free(sten->ststore[i + endpt]);
            }
        }
        
}		/*end update_reg_grid_state*/

/*
*			set_stencil_comps_and_crosses():
*
*	This function loads the components of each point and any crosses
*	that occur in the interior of the stencil.  For the crosses,
*	each stencil entry contains the list of crosses lying on its interior
*	side, and these lists are ordered towards the outer points.
*	Thus sten->crx[0] is meaningless, and sten->crx[1] contains the list
*	of crosses (if any) in order FROM point 0 TO point 1.  Periodic
*	boundaries are essentially invisible (cf Find_rect_comp).
*/

LOCAL	COMPONENT set_stencil_comps_and_crosses(
	int		idir,
	Stencil		*sten)
{
	Wave		*wave = sten->wave;
	Wave		*newwave = sten->newwave;
	int		**icoords = sten->icoords;
	TRI_GRID	*tg = wave_tri_soln(wave)->tri_grid;
	int		i, *nc = sten->nc;
	CRXING		***crx = sten->crx;
	COMPONENT	*comp = sten->comp;

	sten->newcomp = Rect_comp(icoords[0],newwave);
	comp[0] = Rect_comp(icoords[0],wave);
	crx[0][0] = NULL;
	nc[0] = 0;
	sten->hs[0] = NULL;/*TODO REMOVE*/

	for (i = 1; i <= endpt; ++i)
	{
	    sten->hs[-i] = sten->hs[i] = NULL;/*TODO REMOVE*/
	    comp[-i] = Find_rect_comp(idir,icoords[-i],wave);
	    nc[-i] = crossings_in_direction(crx[-i],icoords[-i+1],prev_side,tg);
	    comp[i] = Find_rect_comp(idir,icoords[i],wave);
	    nc[i] = crossings_in_direction(crx[i],icoords[i-1],next_side,tg);
	}
	return sten->newcomp;
}		/*set_stencil_comps_and_crosses*/

/*
*			Find_rect_comp():
*
*	Finds the component of a given point using the macro Rect_comp.
*/

EXPORT	COMPONENT	Find_rect_comp(
	int		idir,
	int		*icoords,
	Wave		*wave)
{
	int		is, imax;
	
	is = icoords[idir];
	imax = gmax[idir] + ubuf[idir];

	if ((is < -lbuf[idir]) || (is >= imax))
	    return ext_comp;
	else
        {
	    return Rect_comp(icoords,wave);
        }	
}		/*Find_rect_comp*/

/*
*			find_outer_states():
*
*	Loads appropriate states for the outer indices of a stencil on either
*	the prev or next sides.  If an interface is crossed, an appropriate
*	state is picked off the interface, and copied out to the end of the
*	stencil (cf copy_states).  Neumann boundaries are also handled
*	appropriately (cf set_neumann_bdry_sten_states).
*/
//#define RK_PRINT(...) {printf("%s:%d: ",__FILE__,__LINE__);printf(__VA_ARGS__);}
//#define RK_PRINT(...) {}
void set_flow_bdry_sten_states(int i, HYPER_SURF *hs, Stencil *sten, 
				CRX_SIDE side, int endpt);
void interp_flow_bdry_states(int start, HYPER_SURF *hs, int side, int endpt, int new_comp,
				Stencil *sten, CRXING *cross);
void assign_reflective_state(Stencil*, int, int, int);

void printf_state(int , Locstate);

LOCAL	void find_outer_states(
	int	 is,
	int	 idir,
	Stencil	 *sten,
	CRX_SIDE side)
{
        Wave	   *wave = sten->wave;
	Front	   *fr = sten->fr;
	Front	   *newfr = sten->newfr;
	COMPONENT  *comp = sten->comp, new_comp = sten->newcomp;
	int	   **icoords = sten->icoords;
	int	   i, indx;
	CRXING	   *crx;
	HYPER_SURF *hs;
	CRXING	   *new_crx[MAX_NUM_CRX];

	crx = NULL;

        //ywall treatment 
        float* coordsx = Rect_coords(icoords[0],wave);
        float GLx = fr->rect_grid->GL[0];
        float GUx = fr->rect_grid->GU[0];
        float dh[3];
        for(i = 0; i < 3; i++)
            dh[i] = fr->rect_grid->h[i];
        
        if((idir == 0 || idir == 1 || idir == 2) && coordsx[2] > 3.9)
        {
            for  (i = 1; i <= endpt; ++i)
	    {
	      indx = i * side;
              sten->st[indx] = sten->st[0];
            }
            return;
        }
        else if(idir == 1 && coordsx[2] < 3.9) //ywall reflection
        {
            //xiaoxue
            bool interside = 1;
            for  (i = 1; i <= endpt; ++i)
	    {
	      indx = i * side;
              if( comp[indx] != 3) interside = 0;
	      if (equivalent_comps(comp[indx],new_comp,newfr->interf) == YES)
	        sten->st[indx] = Rect_state(icoords[indx],wave);
                continue;
            }

            if(interside == 1) return;

            float ywall[2];
            ywall[0] = 0.0;
            ywall[1] = 3.75;
            float* coordsywall;
            int firstGhostCell = 0;
            for(i = 1; i <= endpt; i++)
            {
                indx = i*side;
                if(comp[indx] != 3) {
                    firstGhostCell = indx;
                    break;
                }
            }
            int intCell = firstGhostCell - side;
            int ghostCell = firstGhostCell;
            while(ghostCell <= endpt && ghostCell >= - endpt && intCell <= endpt && ghostCell >= - endpt)
            {
                
                assign_reflective_state(sten, ghostCell, intCell, 1);
                intCell += -side;
                ghostCell += side;
            }
            return;
        }
        else if( idir == 2 || (idir == 0 && side == 1 && coordsx[0] > -7.06 && coordsx[0] < 0.0)) //zwall and ramp reflection
        {
            bool interside = 1;
            for  (i = 1; i <= endpt; ++i)
	    {
	      indx = i * side;
              if( comp[indx] != 3) interside = 0;
	      if (equivalent_comps(comp[indx],new_comp,newfr->interf) == YES)
	        sten->st[indx] = Rect_state(icoords[indx],wave);
                continue;
            }

            if(interside == 1) return;
            
            float zwall[3];
            float zwallx[2];
            zwall[0] = 1.6;
            zwall[1] = 2.4;
            zwall[3] = 3.9;
            zwallx[0] = -7.0;
            zwallx[1] = 0.0;
            int firstGhostCell = 0;
            for(i = 1; i <= endpt; i++)
            {
                indx = i*side;
                if(comp[indx] != 3) {
                    firstGhostCell = indx;
                    break;
                }
            }
            int intCell = firstGhostCell - side;
            int ghostCell = firstGhostCell;
            
            int xdir;
            if(idir == 2)
            {
                if( coordsx[2] > 3.0) xdir = 2;
                else if(coordsx[0] > 0) xdir = 2;
                else  xdir = 3;
            }
            else
            {
                xdir = 3;
            }
            while(ghostCell <= endpt && ghostCell >= - endpt && intCell <= endpt && ghostCell >= - endpt)
            {
                assign_reflective_state(sten, ghostCell, intCell, xdir);
                intCell += -side;
                ghostCell += side;
            }
            return;
        }
        else if( idir == 0 && side == -1)// inflow  
        {
              bool interside = 1;
            for  (i = 1; i <= endpt; ++i)
	    {
	      indx = i * side;
              if( comp[indx] != 3) interside = 0;
	      if ( comp[indx] == 3 )
	        sten->st[indx] = Rect_state(icoords[indx],wave);
                continue;
            }
            if(interside == 1) return;
            
            if(debugging("inflow"))
                printf("coordsx = [%f %f %f] comp = [%d %d %d]\n",coordsx[0], coordsx[1], coordsx[2], comp[-2], comp[-1], comp[0]);
            int inflowindx;
            for(i = 1; i <= endpt; i++)
            {
                indx = i*side;
                if(comp[indx] != 3) {
                    inflowindx = indx + 1;
                    break;
                }
            }
            if(debugging("inflow"))
                printf("inflowindx = %d\n", inflowindx);

            int ghostCell = inflowindx - 1;
            while(ghostCell <= endpt && ghostCell >= - endpt) 
            {
                sten->st[ghostCell] = Rect_state(icoords[inflowindx],wave);
                 if(debugging("inflow"))
                    printf_state(ghostCell,sten->st[ghostCell]);
                ghostCell += side;
            }
            return;
        }
        else if( idir == 0 && coordsx[0] >= GUx - 2*dh[0]) //outflow
        {
            bool interside = 1;
            for  (i = 1; i <= endpt; ++i)
	    {
	      indx = i * side;
              if( comp[indx] != 3) interside = 0;
	      if ( comp[indx] == 2 )
	        sten->st[indx] = Rect_state(icoords[indx],wave);
                continue;
            }
            if(interside == 1) return;
            
            int outflowindx;
            for(i = 1; i <= endpt; i++)
            {
                indx = i*side;
                if(comp[indx] != 3) {
                    outflowindx = indx - 1;
                    break;
                }
            }
            
            int ghostCell = outflowindx + 1;
            while(ghostCell <= endpt && ghostCell >= - endpt)
            {
                sten->st[indx] = Rect_state(icoords[outflowindx],wave);
                ghostCell += side;
            } 
            return;
        }
        else if(coordsx[2] > 3.9)
        {
            //xiaoxue
            bool interside = 1;
            for  (i = 1; i <= endpt; ++i)
	    {
	      indx = i * side;
              if( comp[indx] != 3) interside = 0;
	      if (equivalent_comps(comp[indx],new_comp,newfr->interf) == YES)
	        sten->st[indx] = Rect_state(icoords[indx],wave);
                continue;
            }
            if(interside == 1) return;

            int firstGhostCell = 0;
            for(i = 1; i <= endpt; i++)
            {
                indx = i*side;
                if(comp[indx] != 3) {
                    firstGhostCell = indx;
                    break;
                }
            }
            int ghostCell = firstGhostCell;
            while(ghostCell <= endpt && ghostCell >= - endpt)
            {
                if(idir != 2)
                  assign_reflective_state(sten, ghostCell, 0, idir);
                else
                  sten->st[ghostCell] = sten->st[0];
                ghostCell += side;
            }
            return;
        }
        printf("error! not included in the if-else cased in find_outer_state\n");
        for  (i = 1; i <= endpt; ++i)
	{
	    indx = i * side;

	    if (equivalent_comps(comp[indx],new_comp,newfr->interf) == YES)
	    {
	        sten->st[indx] = Rect_state(icoords[indx],wave);
	        continue;
	    }
	    if (detect_and_load_mix_state(idir,sten,indx))
	        continue;
	    sten->reg_stencil = NO;
	    crx = sten->crx[indx][0];

		if(crx == NULL)
		{
			//copy_states if it is outside the domain.
			if (is+indx < -lbuf[idir] || is+indx >= gmax[idir]+ubuf[idir])
			{
				double *crds = Rect_coords(icoords[indx],wave);
				sten->st[indx] = sten->st[indx-side];
				copy_states(i,side,endpt,new_comp,sten,NULL);
				return;
			}
		}
	    
	    if (crx != NULL)
	    {
			double *crds = Rect_coords(icoords[indx],wave);
	        /*there is a cross between indx and in_one*/

	        (void) crossings_in_direction(new_crx,icoords[(i-1)*side],
	                               (side == PREV) ? prev_side : next_side,
	                               wave_tri_soln(sten->newwave)->tri_grid);

	        hs = crx->hs;

			if ((new_crx[0] == NULL) ||
				(wave_type(hs) != wave_type(new_crx[0]->hs)))
			{
			   /* crx has moved out of stencil OR
				* another curve has moved between indx and the
				* old crx over the course of the time step. */

				if(debugging("sample_speed"))
				{
					interp_states(i,side,endpt,new_comp,idir,sten,crx);
				}
				else
				{
					copy_states(i,side,endpt,new_comp,sten,crx);
				}
	            return;
			}
	        if ((wave_type(hs) == NEUMANN_BOUNDARY) &&
	            (fr->neumann_bdry_state))
	        {
	            set_neumann_bdry_sten_states(i,is,idir,new_comp,
	                                         crx,sten,side);
	            return;
	        }
                /*
	        else if (wave_type(hs) == DIRICHLET_BOUNDARY)
//	                 || wave_type(hs) == FLOW_BOUNDARY)
	        {
                    //fprintf(stdout, "before DIRICHLET_BOUNDARY\n");
		    if(debugging("sample_speed"))
		    {
			    interp_dirichlet_states(i,side,endpt,new_comp,sten,crx);
		    }
		    else
		    { 
			    set_dirichlet_bdry_sten_states(i,hs,sten,side);
		    }
                    interp_flow_bdry_states(i,hs,side,endpt,new_comp,sten,crx);
                    //fprintf(stdout, "after DIRICHLET_BOUNDARY\n");

		    if(NO)
		    {
			    int	i;

			    (void) printf("STENCIL: "); print_Stencil(sten);
			    for (i = -endpt; i <= endpt; ++i)
			    {
				(void) printf("state[%d]:	",i);
				print_states_gas(sten->st[i]);
				    //(*fr->print_state)(sten->st[i]);
			    }
			    fflush(NULL);
		    }

		    return;
	        }
                */
			// NEW BDRY CONDITION DEBUGGING RK
                else if (wave_type(hs) == FLOW_BOUNDARY 
                         || wave_type(hs) == DIRICHLET_BOUNDARY)
		{
		    double *crds = Rect_coords(icoords[indx],wave);

			#if 0 // dirichlet condition
		    if(debugging("sample_speed"))
		    {
			    interp_dirichlet_states(i,side,endpt,new_comp,sten,crx);
		    }
		    else
		    {
			    set_dirichlet_bdry_sten_states(i,hs,sten,side);
		    }

		    return;
#else		
		    set_dirichlet_bdry_sten_states(i,hs,sten,side);
		    interp_flow_bdry_states(i,hs,side,endpt,new_comp,sten,crx);
		    //set_flow_bdry_sten_states(i,hs,sten,side, endpt);

#endif

//		    set_dirichlet_bdry_sten_states(i,hs,sten,side);
//		    return;
		}
	        else
	        {  
	            /* crx is interior */
				double *crds = Rect_coords(icoords[indx],wave);

				if(debugging("sample_speed"))
				{
					interp_states(i,side,endpt,new_comp,idir,sten,crx);
				}
				else
				{
					copy_states(i,side,endpt,new_comp,sten,crx);
				}
					return;
	        }
	    }
	    else
	    {
				double *crds = Rect_coords(icoords[indx],wave);
			int		   gr_side = ERROR;
	        float		   *coords;
	        static const float OFFSET = 0.01; // TOLERANCE

	        coords = Coords(sten->p[indx]);
	        if (is+indx < -lbuf[idir])
			{
				gr_side = 0;
	            coords[idir] = VL[idir] + hdir*OFFSET;
	        }
			else if (is+indx >= gmax[idir]+ubuf[idir])
			{
				gr_side = 1;
	            coords[idir] = VU[idir] - hdir*OFFSET;
	        }
	
			if(gr_side != ERROR &&
			   rect_boundary_type(fr->interf, idir, gr_side) == 
				OPEN_BOUNDARY)
			{
				sten->st[indx] = sten->st[indx-side];
				copy_states(i,side,endpt,new_comp,sten,crx);
				return;
			}
		
			if(debugging("sample_speed"))
			{
				interp_states(i,side,endpt,new_comp,idir,sten,crx);
				return;
			}

			sten->st[indx] = sten->worksp[indx];
			if (!nearest_crossing_state_with_comp(icoords[indx],
	    		Coords(sten->p[indx]),new_comp,
				wave_tri_soln(sten->newwave)->tri_grid,
				&sten->st[indx]))
			{
	            nearest_intfc_state_and_pt(coords,new_comp,
								   sten->fr,sten->newfr,
								   sten->st[indx],
								   Coords(sten->p[indx]),
								   &sten->hs[indx]);
			}
	        
			copy_states(i,side,endpt,new_comp,sten,crx);
	        return;
	    }
	}
}		/*end find_outer_states*/


/*
*			set_neumann_bdry_sten_states():
*
*	Sets stencil states using reflection through a Neumann boundary.
*/

LOCAL	void	set_neumann_bdry_sten_states(
	int	  i,
	int	  is,
	int	  idir,
	COMPONENT new_comp,
	CRXING	  *nbdry,
	Stencil	  *sten,
	CRX_SIDE  side)
{
    	HYPER_SURF *hs = nbdry->hs;
	Wave	   *wave = sten->wave;
	Front	   *fr = sten->fr;
	CRXING	   *rcrx;
	float 	   coords[MAXD], coordsref[MAXD], n[MAXD];
	int	   j, indx, ri;
	size_t	   sizest = fr->sizest;

	for (; i <= endpt; ++i)
	{
	    indx = i * side;
	    if ((*fr->neumann_bdry_state)(Coords(sten->p[indx]),new_comp,
					  nbdry->pt,hs,fr,
					  (POINTER)wave,sten->worksp[indx]))
	    {
	        sten->st[indx] = sten->worksp[indx];
	        for (j = 0; j < dim; ++j)
	    	    Coords(sten->p[indx])[j] = Coords(nbdry->pt)[j];
	        ft_assign(left_state(sten->p[indx]),sten->st[indx],sizest);
	        ft_assign(right_state(sten->p[indx]),sten->st[indx],sizest);
	    }
	    else 
	    {
	        if (sten->crx[indx][0] != nbdry)
	    	    copy_states(i,side,endpt,new_comp,sten,NULL);
	        else
	    	    copy_states(i,side,endpt,new_comp,sten,nbdry);
	        return;
	    }

	    if (is_bdry_hs(hs))
	    {
	        /* TODO:
		 * Currently this section of set_neumann_bdry_sten_states() is
	         * only implemented for rectangular Neumann boundaries
	         * as indicated by the above test.  Interior
	         * or oblique boundaries do not use this boundary
	         * state function. */

	        /* Have set the first state, now need to check for another
	         * cross on the other side of the Neumann boundary. Some
	         * technical issues arise due to the way neumann_bdry_state()
	         * reflects coords.
	         */

	        if (i == endpt)
		    break;
	        ri = (side == PREV) ? -indx - 2*is - 1 :
			              -indx + 2*(gmax[idir] - is) - 1;
	        if (side == PREV)
	        {
	    	    rcrx = (ri >= 0) ? sten->crx[ri+1][0] :
				       sten->crx[ri][sten->nc[ri]-1];
	        }
	        else if (side == NEXT)
	        {
	    	    rcrx = (ri <= 0) ? sten->crx[ri-1][0] :
	    			       sten->crx[ri][sten->nc[ri]-1];
	        }
	        if (rcrx != NULL) 
	        {
	    	    indx = (++i) * side;
	    	    for (j = 0; j < dim; ++j) 
	    	        coords[j] = Coords(rcrx->pt)[j];
	    	    if ((reflect_pt_about_Nbdry(coords,coordsref,n,new_comp,
					        nbdry->hs,fr))
		        && 
		        (*fr->neumann_bdry_state)(coordsref,new_comp,nbdry->pt,
					          hs,fr,(POINTER)wave,
					          sten->worksp[indx]))
		    {
		        sten->st[indx] = sten->worksp[indx];
		        for (j = 0; j < dim; ++j)
		    	    Coords(sten->p[indx])[j] = coordsref[i];
		        ft_assign(left_state(sten->p[indx]),sten->st[indx],sizest);
		        ft_assign(right_state(sten->p[indx]),sten->st[indx],sizest);
		    }
		    copy_states(i,side,endpt,new_comp,sten,NULL);
		    return;
	        }
	    }
	}
}		/*end set_neumann_bdry_sten_states*/

/*
*			set_dirichlet_bdry_sten_states():
*
*	Sets stencil states through a Dirichlet boundary.
*/

//#bjet2
LOCAL	void	set_dirichlet_bdry_sten_states(
	int	   i,
	HYPER_SURF *hs,
	Stencil	   *sten,
	CRX_SIDE   side)
{
	Wave		*wave = sten->wave;
	Front		*fr = sten->fr;
	int		indx; 
	size_t		sizest = fr->sizest;

	for (; i <= endpt; ++i)
	{
	    indx = i * side;
	    sten->st[indx] = sten->worksp[indx];
	    evaluate_dirichlet_boundary_state(Coords(sten->p[indx]),hs,fr,wave,
	    				     sten->st[indx]);
	    ft_assign(left_state(sten->p[indx]),sten->st[indx],sizest);
	    ft_assign(right_state(sten->p[indx]),sten->st[indx],sizest);
	}
}		/*end set_dirichlet_bdry_sten_states*/

/*
*			copy_states():
*
*	We have crossed a non-reflecting curve and are ready to simply copy an
*       appropriate state out to the end of the stencil.  If the cross is not
*	NULL, then it is assumed that the start state has not been set, and
*	an appropriate state is loaded and then copied.
*/

LOCAL	void copy_states(
	int	 start,
	CRX_SIDE side,
	int	 endpt,
	int	 new_comp,
	Stencil	 *sten,
	CRXING	 *cross)
{
	int		i, j, indx, in_one;	
	size_t		sizest = sten->fr->sizest;

	if (cross != NULL)
	{
	    indx = start*side;
	    sten->st[indx] = state_with_comp(cross->pt,cross->hs,new_comp);
	    if (sten->st[indx] != NULL)
	    {
	        for (j = 0; j < dim; ++j)
	    	    Coords(sten->p[indx])[j] = Coords(cross->pt)[j];
	        ft_assign(left_state(sten->p[indx]),sten->st[indx],sizest);
	        ft_assign(right_state(sten->p[indx]),sten->st[indx],sizest);
	    }
	    else
	    {
	        sten->st[indx] = sten->worksp[indx];
	        nearest_intfc_state_and_pt(Coords(sten->p[indx]),new_comp,
					   sten->fr,sten->newfr,sten->st[indx],
					   Coords(sten->p[indx]),
					   &sten->hs[indx]);
	    }
	}

	for (i = start+1; i <= endpt; ++i)
	{
	    indx = i*side;
	    in_one = indx - side;
	    sten->st[indx] = sten->st[in_one];
	    for (j = 0; j < dim; ++j)
	    	Coords(sten->p[indx])[j] = Coords(sten->p[in_one])[j];
	    ft_assign(left_state(sten->p[indx]),sten->st[indx],sizest);
	    ft_assign(right_state(sten->p[indx]),sten->st[indx],sizest);
	}
	return;
}		/*end copy_states*/

float	interp_kappa(
	POINT	*p,
	TRI	*t)
{
	float	kappa, f[3], *p0,*p1,*p2;
	int	i;

	p0 = Coords(Point_of_tri(t)[0]);
	p1 = Coords(Point_of_tri(t)[1]);
	p2 = Coords(Point_of_tri(t)[2]);

	linear_interp_coefs_three_pts(f,Coords(p),p0,p1,p2);
	
	kappa = 0.0;
	for(i=0; i<3; i++)
	{
	    if(f[i] < 0.0)
		f[i] = 0.0;
	    if(f[i] > 1.0)
		f[i] = 1.0;

	    kappa += Point_of_tri(t)[i]->curvature*f[i];
	    
	    //printf("%15.8e  ", Point_of_tri(t)[i]->curvature);
	}

	//printf("  %d %d %d  ", Boundary_point(Point_of_tri(t)[0]), 
	//		       Boundary_point(Point_of_tri(t)[1]),
	//		       Boundary_point(Point_of_tri(t)[2]));
	
	return kappa;
}


void	interp_geom(
	float		*kappa,
	float		*nor,
	POINT		*p,
	TRI		*t,
	HYPER_SURF	*hs,
	Front		*fr)
{
	float	len, f[3], v[3], *p0,*p1,*p2;
	POINT	*pt;
	int	i, j;

	if(t == NULL)
	{
	    *kappa = 0.0;
	    return;
	}

	p0 = Coords(Point_of_tri(t)[0]);
	p1 = Coords(Point_of_tri(t)[1]);
	p2 = Coords(Point_of_tri(t)[2]);

	linear_interp_coefs_three_pts(f,Coords(p),p0,p1,p2);
	
	for(j=0; j<3; j++)
	    nor[j] = 0.0;
	*kappa = 0.0;
	
	for(i=0; i<3; i++)
	{
	    if(f[i] < 0.0)
		f[i] = 0.0;
	    if(f[i] > 1.0)
		f[i] = 1.0;

	    pt = Point_of_tri(t)[i];
	    normal(pt, Hyper_surf_element(t), hs, v, fr);
	    
	    *kappa += pt->curvature*f[i];
	    for(j=0; j<3; j++)
		nor[j] += v[j]*f[i];
	    
	    //printf("%15.8e  ", Point_of_tri(t)[i]->curvature);
	}

	len = Mag3d(nor);
	for(j=0; j<3; j++)
	    nor[j] /= len;

	//printf("  %d %d %d  ", Boundary_point(Point_of_tri(t)[0]), 
	//		       Boundary_point(Point_of_tri(t)[1]),
	//		       Boundary_point(Point_of_tri(t)[2]));
}


int	interp_st_in_direction(Locstate,int*,int,int,Stencil*);

bool	extrp_st_normal_edge(Locstate,float*,float*,int,int,Stencil*,Front*);

void	get_st_from_exterp(Locstate,Locstate,Locstate,
			Stencil*,int, HYPER_SURF*,float,float*);

// For point indx comp != new_comp
LOCAL	void interp_states(
	int	 start,
	CRX_SIDE side,
	int	 endpt,
	int	 new_comp,
	int	 idir,
	Stencil	 *sten,
	CRXING	 *cross)
{
	int		i, j, indx, in_one, newcomp;	
	size_t		sizest = sten->fr->sizest;
	INTERFACE	*intfc = sten->fr->interf;
	Locstate	s1, s2;
	int		**icoords = sten->icoords, inone;
	CRXING		*new_crx[MAX_NUM_CRX];
	HYPER_SURF	*hs;
	Front		*fr;
	float		kappa, len, nor[3], fnor[3], fnor1[3], *crds1, *crds2, *pt;
	static Locstate	sll = NULL, srr = NULL, sttmp;

	DEBUG_ENTER(interp_states)
	
	if(sll == NULL)
	{
            alloc_state(intfc,&sll,sizest);
            alloc_state(intfc,&srr,sizest);
            alloc_state(intfc,&sttmp,sizest);
	}
	//TMP limiter

	indx = start*side;
	newcomp = Rect_comp(sten->icoords[indx], sten->newwave);

	if(cross != NULL)
	{
	    for (j = 0; j < 3; ++j)
		Coords(sten->p[indx])[j] = Coords(cross->pt)[j];
	}
	else
	{
	    for (j = 0; j < 3; ++j)
		Coords(sten->p[indx])[j] = Coords(sten->p[indx-side])[j];
	}

	//printf("#interp_state %d %d\n", newcomp, new_comp);
	
	if(newcomp == new_comp)
	{
	    //get_st_from_interp(sten->worksp[indx], sten->icoords[indx], sten);

if(debugging("norexp"))
{
	    if(!extrp_st_along_normal(sten->worksp[indx],Coords(sten->p[indx]),sten,sten->fr))
		get_st_from_ghost(sten->worksp[indx], sten->icoords[indx], sten);
}
else
{
	    get_st_from_ghost(sten->worksp[indx], sten->icoords[indx], sten);
}

	}
	else
	{
	    inone = (start-1)*side;
	    
	    s1 = sten->st[inone];
	    s2 = Rect_state(icoords[indx], sten->wave);
	    crds1 = Rect_coords(icoords[inone],sten->wave);
	    crds2 = Rect_coords(icoords[indx],sten->wave);
	    
	    difference(crds2, crds1, nor, 3);

	    len = Mag3d(nor);
	    for(j=0; j<3; j++)
		nor[j] /= len;

	    if(debugging("interp_npt"))
		print_general_vector("s2crds", 
		    Rect_coords(sten->icoords[indx],sten->wave), 3, "\n");

	    ft_assign(fnor, nor, 3*FLOAT);

	    if(cross != NULL)
	    {
		hs = cross->hs;
		//kappa = interp_kappa(cross->pt, cross->tri);
		
		interp_geom(&kappa, fnor, cross->pt, cross->tri, hs, sten->fr);
		pt = Coords(cross->pt);
		fr = sten->fr;
	    }
	    else
	    {
		//newcomp != new_comp, there must be a crx between two nodes.
		(void) crossings_in_direction(new_crx,icoords[(start-1)*side],
	                               (side == PREV) ? prev_side : next_side,
	                               wave_tri_soln(sten->newwave)->tri_grid);
		if(new_crx[0] == NULL)
		{
		    printf("ERROR interp_states, no crx in newwave.\n");
		    (void) printf("STENCIL: "); 
		    print_Stencil(sten);
		    clean_up(ERROR);
		}
		hs = new_crx[0]->hs;
		//kappa = interp_kappa(new_crx[0]->pt, new_crx[0]->tri);

		interp_geom(&kappa, fnor, new_crx[0]->pt, new_crx[0]->tri, hs, sten->newfr);
		pt = Coords(new_crx[0]->pt);
		fr = sten->newfr;
	    }

	    if(debugging("norexp"))
	    {
		bool	status1, status2;

		status1 = extrp_st_normal_edge(sll, pt, fnor, 2, idir, sten, fr);
		for(i=0; i<3; i++)
		    fnor1[i] = -fnor[i];
		status2 = extrp_st_normal_edge(srr, pt, fnor1, 3, idir, sten, fr);

		if(NO && status1 && status2)
		{
		    printf("#sll");
		    print_states_gas(sll);
		    if(new_comp == 2)
			print_states_gas(s1);
		    else
			print_states_gas(s2);
		    
		    printf("#srr");
		    print_states_gas(srr);
		    if(new_comp == 2)
			print_states_gas(s2);
		    else
			print_states_gas(s1);
		    
		    //clean_up(ERROR);
		}
		//else
		//{
		//    print_general_vector("#fpt=", pt, 3, "\n");
		//}
		
		//assume fnor pointing to comp 2

if(debugging("revjump"))
{
		if(new_comp == 3)
		    for(i=0; i<3; i++)
			nor[i] = -fnor[i];
		else
		    ft_assign(nor, fnor, 3*FLOAT);
}
else
{
		if(new_comp == 2)
		    for(i=0; i<3; i++)
			nor[i] = -fnor[i];
		else
		    ft_assign(nor, fnor, 3*FLOAT);
}

		/*
		if(status1 && status2)
		{
		    if(new_comp == 2)
		    {
			s1 = sll;
			s2 = srr;
		    }
		    else
		    {
			s1 = srr;
			s2 = sll;
		    }
		}
		*/
	    }
	    
	    if(debugging("gfm2ndr"))
	    {
		float		v[3], fac1, fac2;
		int		fg, comp2;

		difference(pt, crds1, v, 3);
		fac1 = Dot3d(v, nor);
		
		difference(crds2, pt, v, 3);
		fac2 = Dot3d(v, nor);
		
		fac1 = fac1/(fac1+fac2);
		fac2 = 1.0 - fac1;
	    
		//printf("orig s12 %d  %f  %f\n", new_comp, fac1, fac2);
		//print_states_gas(s1);
		//print_states_gas(s2);
		
		fg = interp_st_in_direction(sttmp,icoords[indx],idir,new_comp,sten);
		if(fg <= 1)
		{
		    ft_assign(sll, s1, sizest);
		    bi_interpolate_intfc_states(intfc, fac2, fac1,
			crds1, s1, crds2, sttmp, sll);

		    s1 = sll;
		}
		
		//printf("fg %d\n", fg);
		
		comp2 = new_comp == 2 ? 3 : 2;
		fg = interp_st_in_direction(sttmp,icoords[inone],idir,comp2,sten);
		if(fg <= 1)
		{
		    ft_assign(srr, s2, sizest);
		    bi_interpolate_intfc_states(intfc, fac2, fac1,
			crds1, sttmp, crds2, s2, srr);
	
		    s2 = srr;
		}
	
		//printf("fg %d\n", fg);
		//printf("slr  \n");
		//print_states_gas(s1);
		//print_states_gas(s2);
	
		//clean_up(ERROR);
	    }

	    //printf("#kappa %15.8e\n", kappa);

	    get_st_from_exterp(sten->worksp[indx], s1, s2, sten, 
	        new_comp, hs, kappa, nor);

	    //printf("#extr st\n");
	    //add_to_debug("interp_npt");

	    //top fixed layer 
	    //if((idir == 0 || idir == 1) && 
	    //    crds2[2] > 4.0-0.012 && crds2[2] < 4.0+0.012)
	    //{
		//ft_assign(sten->worksp[indx],s1,sizest);
	    //}
	}

	sten->st[indx] = sten->worksp[indx];
	
	//add_to_debug("interp_npt");
	//print_int_vector("icrds=", sten->icoords[indx], 3, "\n");
	//printf("#ghost fi\n");
	
	ft_assign(left_state(sten->p[indx]),sten->st[indx],sizest);
	ft_assign(right_state(sten->p[indx]),sten->st[indx],sizest);

	for (i = start+1; i <= endpt; ++i)
	{
	    indx = i*side;
	    in_one = indx - side;
	    sten->st[indx] = sten->st[in_one];
	    for (j = 0; j < dim; ++j)
	    	Coords(sten->p[indx])[j] = Coords(sten->p[in_one])[j];
	    ft_assign(left_state(sten->p[indx]),sten->st[indx],sizest);
	    ft_assign(right_state(sten->p[indx]),sten->st[indx],sizest);
	}
	
	DEBUG_LEAVE(interp_states)
	return;
}


#if defined(TWOD)

/*                  h_symmetry_states_printout
*
*       This function print out the front states of two symmetric
*       points w.r.t. the center line on the interface. It's for
*       the symmetric testing of the single mode interface.
*       First we find the two bonds which cross these two symmetric
*       vertical lines, then locate the coordinates of the two points
*       of the intersections.  The left and right states are printed
*       on these two points, in order to check the symmetry on the front.
*/

LOCAL   int	 h_symmetry_states_printout(
        Front           *front,
        Wave		*wave)
{
        int             i, k, dim = wave->rect_grid->dim;
        int             wave_type;
        int             imin[MAXD], imax[MAXD];
        HYPER_SURF      *hs;
        CURVE           *c;
        BOND            *b;
        float           lp[2], rp[2], s0, s1, s2, s3;
        float           dh[MAXD], sgn_x_l, sgn_x_r;
        Locstate        sl_l, sl_r, sr_l, sr_r;
        INTERFACE       *interf = front->interf;
        COMPONENT       pcomp = positive_component(hs);
        COMPONENT       ncomp = negative_component(hs);

        alloc_state(intfc,&sl_l,front->sizest);
        alloc_state(intfc,&sl_r,front->sizest);
        alloc_state(intfc,&sr_l,front->sizest);
        alloc_state(intfc,&sr_r,front->sizest);

        for (i = 0; i < dim; i++)
        {
            dh[i] = wave->rect_grid->h[i];
            imin[i] = 0;        imax[i] = wave->rect_grid->gmax[i];
            if (wave->rect_grid->lbuf[i] == 0) imin[i] += 1;
            if (wave->rect_grid->ubuf[i] == 0) imax[i] -= 1;
        }

        printf("dh[0] = %lf, dh[1] = %lf\n",dh[0],dh[1]);

        if (make_bond_comp_lists(interf) == FUNCTION_FAILED)
        {
                (void) printf("WARNING in make_bond_lists(), ");
                (void) printf("make_bond_listsd() failed\n");
                return FUNCTION_FAILED;
        }

        k = 0;      

        lp[0] = ((imin[0] + imax[0]) / 4) * (dh[0]);

        rp[0] = (((imin[0] + imax[0]) * 3 / 4) - 1) * (dh[0]);


        (void) next_bond(interf,NULL,NULL);

        while (next_bond(interf,&b,&c))
        {
                if(wave_type(c) < FIRST_SCALAR_PHYSICS_WAVE_TYPE) continue;

	        hs = Hyper_surf(c);
	       	pcomp = positive_component(hs);
         	ncomp = negative_component(hs);

                s0 = (Coords(b->end)[0] - lp[0]);
                s1 = (Coords(b->start)[0] - lp[0]);
                s2 = (Coords(b->end)[0] - rp[0]);
                s3 = (Coords(b->start)[0] - rp[0]);

                sgn_x_l = s0 * s1;
                sgn_x_r = s2 * s3;

                if(sgn_x_l > 0.0 && sgn_x_r > 0.0) continue;

                if(s0 == 0.0)
                {
                  lp[1] = Coords(b->end)[1];

                  k = k + 1;
                }

                if(sgn_x_l < 0.0)
                {
                  lp[1] = ((Coords(b->end)[1] - Coords(b->start)[1]) *
                  (Coords(b->start)[0] - lp[0]) /
                  (Coords(b->end)[0] - Coords(b->start)[0])) +
                   Coords(b->start)[1];

                  k = k + 1;
                }

                if(s2 == 0.0)
                {
                  rp[1] = Coords(b->end)[1];

                  k = k + 1;
                }

                if(sgn_x_r < 0.0)
                {
                  rp[1] = ((Coords(b->end)[1] - Coords(b->start)[1]) *
                  (Coords(b->start)[0] - rp[0]) /
                  (Coords(b->end)[0] - Coords(b->start)[0])) +
                   Coords(b->start)[1];

                  k = k + 1;
                }

        }

        printf("left point = (%lf,%lf)\n",lp[0],lp[1]);

        printf("right point = (%lf,%lf)\n",rp[0],rp[1]);

        printf("k = %d\n",k);

        hyp_solution(lp,ncomp,hs,NEGATIVE_SIDE,front,wave,sl_l,NULL);

        hyp_solution(lp,pcomp,hs,POSITIVE_SIDE,front,wave,sl_r,NULL);

        hyp_solution(rp,ncomp,hs,NEGATIVE_SIDE,front,wave,sr_l,NULL);

        hyp_solution(rp,pcomp,hs,POSITIVE_SIDE,front,wave,sr_r,NULL);

        (*front->print_state)(sl_l);
        (*front->print_state)(sl_r);
        (*front->print_state)(sr_l);
        (*front->print_state)(sr_r);

}

#endif /* defined(TWOD) */



