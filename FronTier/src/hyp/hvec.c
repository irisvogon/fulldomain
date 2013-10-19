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
*				hvec.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains a uni_arrayized row-by-row (column-by-column) driver for a
*	finite difference stencil solution of the hyperbolic equations
*	within a single component of a tracked interface problem.
*	Point sources in the region are found and passed appropriately.
*/


#include <hyp/hdecs.h>


LOCAL void sweep_comp_segments(int,int*,float*,Wave*,Wave*,Front*,Front*,
                         int*,int,int,float,float,int*);

/*
*			hyp_reg_vec():
*
*		    single directional sweep.
*
*	Perform a regular grid sweep first followed by a correction
*	sweep at the irregular points.
*/
void printf_state(int , Locstate);
/*ARGSUSED*/
EXPORT	void hyp_reg_vec(
	int		swp_num,
	int		*iperm,	/* sweep determined coord permutation	*/
	float		dh,	/* space increment */
	float		dt,	/* time increment */
	Wave		*wv,
	Wave		*newwv,
	Front		*fr,
	Front		*newfr,	/* newfr needed if hlw is to support
				   changing top. */
	COMPONENT	max_comp)
{
	RECT_GRID	*gr = fr->rect_grid;
	float		dir[MAXD];
	int		imin[3], imax[3];
	int		i0min, i0max;
	int		vsize;
	int		dim = gr->dim;
	int		icoords[MAXD];
	int		i, idirs[MAXD];
	int		nrad = vsten_radius(wv);

	for (i = 0; i < dim; ++i)
	{
	    idirs[i] = iperm[(i+swp_num)%dim];
	    dir[i] = 0.0;
	}
	dir[idirs[0]] = 1.0;
	debug_print("hyp_vec","Entered hyp_reg_vec(), dir = %d, swp_num = %d\n",
	      idirs[0],swp_num);

        set_sweep_limits(wv,swp_num,idirs,imin,imax);
	//#bjet2
	set_limits_for_open_bdry(wv,fr,swp_num,idirs,imin,imax);
	
	i0min = imin[0];
	i0max = imax[0];
	if (gr->lbuf[idirs[0]] > 0)
	    i0min = max(i0min-nrad, -gr->lbuf[idirs[0]]);
	if (gr->ubuf[idirs[0]] > 0)
	    i0max = min(i0max+nrad, gr->gmax[idirs[0]]+gr->ubuf[idirs[0]]);
        

        /*     if(idirs[0] == 1) //ysweep
        {
              i0min = -nrad;
              i0max = gr->gmax[idirs[0]] + nrad;
        }*/
        
	vsize = i0max - i0min;

	alloc_phys_vecs(wv,vsize);

	icoords[idirs[0]] = 0;	/* value not used */

#if defined(TIME_HVEC)
	start_clock("Finite difference step");
#endif /* defined(TIME_HVEC) */

	switch (dim)
	{
	case 1:
	    sweep_comp_segments(swp_num,iperm,dir,wv,newwv,fr,newfr,
	    		icoords,i0min,i0max,dt,dh,idirs);
	    break;
#if defined(TWOD)
	case 2:
	{
	    int	i1;
	    int	i1min, i1max;

	    i1min = imin[1];	i1max = imax[1];

	    for( i1 = i1min;  i1 < i1max;  ++i1 )
	    {
	        icoords[idirs[1]] = i1;
	        sweep_comp_segments(swp_num,iperm,dir,wv,newwv,fr,
			newfr,icoords,i0min,i0max,dt,dh,idirs);		
	    }
	    break;
	}
#endif /* defined(TWOD) */
#if defined(THREED)
	case 3:
	{
	    int	i1, i2;
	    int	i1min, i1max, i2min, i2max;

	    i1min = imin[1];	i1max = imax[1];
	    i2min = imin[2];	i2max = imax[2];
	    for( i2 = i2min;  i2 < i2max;  ++i2 )
	    {
	        icoords[idirs[2]] = i2;
	   	for( i1 = i1min;  i1 < i1max;  ++i1 )
	   	{
	            icoords[idirs[1]] = i1;
	            sweep_comp_segments(swp_num,iperm,dir,wv,newwv,
		    	    fr,newfr,icoords,i0min,i0max,dt,dh,idirs);
	   	}
	    }
	    break;
	}
#endif /* defined(THREED) */
	}
#if defined(TIME_HVEC)
	stop_clock("Finite difference step");
#endif /* defined(TIME_HVEC) */

	free_phys_vecs(wv);
	if (debugging("hyp_vec"))
	{
	    (void) printf("New wave states after %d regular sweep:\n",idirs[0]);
	    (*wv->show_wave_states)(newwv);
	}
	
	/*  UPDATE THE IRREGULAR GRID */
#if defined(TIME_HVEC)
	start_clock("sweep on irregular grid");
#endif /* defined(TIME_HVEC) */
  /*      printf("after vector_FD sweep in dir %d!\n", idirs[0]);

        for(int kk = imin[2]; kk < imax[2]; kk++)
        {
            for( int jj = imin[1] ; jj < imax[1]; jj++)
            {
                for(int ii = i0min; ii < i0max; ii++)
                {
                    icoordsx[idirs[0]] = ii;
                    icoordsx[idirs[1]] = jj;
                    icoordsx[idirs[2]] = kk;
                    if(ii == 0 && jj == 7 && kk == 11){
                    xstate = Rect_state(icoordsx,wv);
                    coordsx = Rect_coords(icoordsx,wv);
                    xstate2 = Rect_state(icoordsx,newwv);
                    printf(" icoords = [%d %d %d] coords = [%f %f %f]\n", icoordsx[0], icoordsx[1],icoordsx[2], coordsx[0], coordsx[1], coordsx[2]);
                    if(xstate!= NULL)
                    {
                        printf("wv:    ");
                        printf_state(icoordsx[idirs[0]],xstate);
                    }
                    if(xstate2 != NULL)
                    {
                        printf("newwv: ");
                        printf_state(icoordsx[idirs[0]],xstate2);
                    }
                    }
                }
            }
        }
*/
        hyp_npt(swp_num,iperm,dh,dt,wv,newwv,fr,newfr,max_comp);
        
#if defined(TIME_HVEC)
	stop_clock("sweep on irregular grid");
#endif /* defined(TIME_HVEC) */

	if (debugging("hyp_vec"))
	{
	    (void) printf("New wave states after %d irregular sweep:\n",
			  idirs[0]);
	    (*wv->show_wave_states)(newwv);
	}

	debug_print("hyp_vec","Leaving hyp_reg_vec(), dir = %d\n",idirs[0]);
}			/*end hyp_reg_vec*/


LOCAL void sweep_comp_segments(
	int             swp_num,
        int             *iperm,
	float           *dir,
        Wave            *wv,
        Wave            *newwv, 
        Front           *fr,
        Front           *newfr,
	int		*icoords,
	int		i0min,
	int		i0max,
        float           dt,
        float           dh, 
	int		*idirs)
{
	COMPONENT	comp;
	int		seg_min,seg_max;
	seg_min = i0min;
	while (seg_min != i0max)
	{
	    icoords[idirs[0]] = seg_min;
	    comp = Rect_comp(icoords,wv);
	    for (seg_max = seg_min+1; seg_max < i0max; seg_max++)
	    {
	    	icoords[idirs[0]] = seg_max;
		if (comp != Rect_comp(icoords,wv))
		    break;
	    }

	    vec_solver(swp_num,iperm,dir,wv,newwv,fr,newfr,
	    		icoords,seg_min,seg_max,dt,dh);
            
	    seg_min = seg_max;
	}
        
}	/* end sweep_comp_segment */
