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
*				dscatter.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains the parallel driver initializing routines, as called in main,
*	and the principle or default initializing subroutines.
*/


#include <driver/ddecs.h>


/*
*               clip_front_to_rect_boundary_type():
*/

EXPORT	void	clip_front_to_rect_boundary_type(
	Front	*front)
{
	INTERFACE *intfc = front->interf;
	PP_GRID	  *pp_grid = front->pp_grid;
	RECT_GRID *t_gr, *c_gr, *zoom_gr, Dual_grid;
	int	  i, tgmax[MAXD];

	zoom_gr = &pp_grid->Zoom_grid;

	/* clip virtual boundaries */
	for (i = 0; i < intfc->dim; i++)
	{
	    if (!buffered_boundary_type(rect_boundary_type(intfc,i,0)))
		zoom_gr->lbuf[i] = 0;
	    if (!buffered_boundary_type(rect_boundary_type(intfc,i,1)))
		zoom_gr->ubuf[i] = 0;
	}

	/* set topological and computational grids */
	t_gr = &topological_grid(intfc);
	c_gr = computational_grid(intfc);
	set_rect_grid(zoom_gr->L,zoom_gr->U,zoom_gr->GL,
		      zoom_gr->GU,zoom_gr->lbuf,
		      zoom_gr->ubuf,zoom_gr->gmax,
		      intfc->dim,&t_gr->Remap,zoom_gr);

	for (i = 0; i < intfc->dim; i++)
	{
	    float h = t_gr->h[i];
	    tgmax[i] = irint((zoom_gr->VU[i]-zoom_gr->VL[i])/h);
	}

	if(intfc->dim == 3)
	{
	    //#bjet2  WARN  assume the top grid is square.
	    copy_rect_grid(c_gr,zoom_gr);
	    copy_rect_grid(front->rect_grid,zoom_gr);
	    set_dual_grid(&Dual_grid,c_gr);
	    set_expanded_grid(&Dual_grid,t_gr);
	    (void) adjust_top_grid_for_square(t_gr,zoom_gr);
	}
	else
	{
	    set_rect_grid(zoom_gr->VL,zoom_gr->VU,zoom_gr->GL,
	    	      zoom_gr->GU,NOBUF,NOBUF,tgmax,
	    	      intfc->dim,&t_gr->Remap,t_gr);
	    (void) adjust_top_grid_for_square(t_gr,zoom_gr);
	    copy_rect_grid(c_gr,zoom_gr);
	    copy_rect_grid(front->rect_grid,zoom_gr);
	}

	if (intfc->dim == 3 && pp_numnodes() > 1)
	{
	    delete_outside_surface(intfc);	
	    communicate_default_comp(front);
	}
}		/*end clip_front_to_rect_boundary_type*/

