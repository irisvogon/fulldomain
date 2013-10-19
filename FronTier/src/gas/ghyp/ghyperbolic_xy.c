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
*				ghyperbolic.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains g_tri_integral() and g_quad_integral, set to funciton
*		pointers in the wave structure.
*
*	Contains routines used in the hyperbolic step.
*
*/

#include <ghyp/ghyp.h>

#define VAPOR_COMP    4
#define AMBIENT_COMP  3
#define   tan10                 0.176327
#define   cos10                 0.984807
#define   sin10                 0.173648 
#define   ramp_x0                -7.0
#define   ramp_x1               ( -7 + 0.8/tan10)
#define   ramp_z0               1.6
#define   ramp_z1               2.4

	/* LOCAL Function Declarations */
LOCAL	bool	too_many_bad_gas_states(const char*,int,Front*);
LOCAL	void	g_init_obstacle_states(Vec_Gas*,int,int);
LOCAL	int	detect_ambient_vapor_mix(int,Stencil*);
LOCAL   int     local_LF_npt_tang_solver_switch(float,Tan_stencil*,Front*);
LOCAL   int     on_what = 0;

// TMP RK for jet
void    get_constant_state(
        Locstate        sl,
        int             comp,
        float           *coords,
        float           t);

bool	check_change_inflow_st(Locstate,float*);

/*ARGSUSED*/
EXPORT	void point_FD(
	float		dh,
	float		dt,
	Locstate	ans,
	const float	*dir,
	int		swp_num,
	int		*iperm,
	int		*index,
	Stencil		*sten)
{
        printf("enter point_FD\n");
        Front		*fr = sten->fr;
	Front		*newfr = sten->newfr;
	RECT_GRID	*gr = fr->rect_grid;
	Locstate	st, *state;
	Wave 		*wave = sten->wave;
	Wave 		*newwave = sten->newwave;
	float		*rho, *en_den;
	float		*m[MAXD];
	float		*vel_field[MAXD];
	float		**coords;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	float		*vacuum_dens;
	float		*min_pressure;
	float		*min_energy;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	static float	alpha;
	int		i, j;
	int		dim = gr->dim;
	int		idirs[MAXD];
	static Vec_Gas	*vst = NULL;
	static Vec_Src *src = NULL;
	static int	endpt, npts;
	static float    **Q = NULL;
	static Locstate	inst = NULL;
	static Locstate instx = NULL;
	bool		smooth_vel;
	int		imin[3],imax[3];
	int		*gmax = gr->gmax;
	int leng = 0;
	bool upboundary;
	float *coordsx = Coords(sten->pstore[0]);
	int **icoords = sten->icoords;
	int icoordsp[3];
	int comp;
	int  *lbuf, *ubuf;
	lbuf = gr->lbuf;
	ubuf = gr->ubuf;

	if (is_obstacle_state(sten->st[0])) 
	{
	    g_obstacle_state(ans,sten->fr->sizest);
	    return;
	}
	
	for( i = 0; i < dim; i++)    idirs[i] = iperm[(i+swp_num)%3];

	
	for (i = 0; i < dim; ++i)
	{
	    imin[i] = 0;
	    imax[i] = gmax[i];
	    if(lbuf[i] > 0) imin[i] -= lbuf[i];
	    if(ubuf[i] > 0) imax[i] += ubuf[i];	
	}
	
	if (vst == NULL) 
	{
	    /* Assumes stencil size never changes */

	    npts = sten->npts;
	    endpt = npts/2;
//	    alloc_phys_vecs(wave,npts);
	    alloc_phys_vecs(wave,12);
	    vst = g_wave_vgas(wave);
	    src = g_wave_vsrc(wave);
	    g_wave_vgas(wave) = NULL;
	    g_wave_vsrc(wave) = NULL;
	    bi_array(&Q,3,3,FLOAT);
	    vst->Q = (const float* const*)Q;
	    if (is_rotational_symmetry())
		alpha = rotational_symmetry();
	}

	if(inst == NULL)
	    scalar(&inst, fr->sizest);
	if(instx == NULL) scalar(&instx,fr->sizest);

	if (RegionIsFlowSpecified(ans,sten->st[0],Coords(sten->p[0]),
				  sten->newcomp,sten->newcomp,fr))
	    return;

	if((idirs[0] == 0 || idirs[0] == 1) && Coords(sten->pstore[2])[2] > 3.90)
	{
	      	ft_assign(instx,sten->ststore[2],fr->sizest);
		get_constant_state(instx, 3, Coords(sten->pstore[2]),fr->time);
		float* coordsxx = Coords(sten->pstore[2]);
		Dens(ans) = Dens(instx);
		if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
		{
		        for(j = 0; j < Params(sten->ststore[2])->n_comps; j++)
			 pdens(ans)[j] = pdens(instx)[j];
		}
		Energy(ans) = Energy(instx);
		for( i = 0; i < dim; i++)
			Mom(ans)[i] = Mom(instx)[i];
		Set_params(ans,sten->ststore[2]);
		set_type_of_state(ans,GAS_STATE);
		return;
	}

        /* Hard wired code */
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(YES == detect_ambient_vapor_mix(iperm[swp_num],sten))
            {
                Gas_param   **prms = NULL, *vaporprms = NULL;
                uint64_t    igas;

                return_params_list(&prms);
                for (igas = 1; prms && *prms; ++igas, ++prms)
                {
                    if (igas == VAPOR_COMP-1)
                    {
                        vaporprms = *prms;
                        break;
                    }
                }
                return_params_list(&prms);
                for (igas = 1; prms && *prms; ++igas, ++prms)
                {
                    if (igas == AMBIENT_COMP-1)
                        break;
                }
                for(i = -endpt; i <= endpt; ++i)
                {
                    {
                        if(Params(sten->st[i]) == vaporprms)
                        {
                            pdens(sten->st[i])[0] = 0.0;
                            pdens(sten->st[i])[1] = Dens(sten->st[i]);
                            Params(sten->st[i]) = *prms;
                        }
                    }
                }
            }
        }

	clear_Vec_Gas_set_flags(vst);
	for (i = 0; i < 3; ++i)
	    for (j = 0; j < 3; ++j)
	        Q[i][j] = 0.0;
	for (i = 0; i < dim; ++i)
	{
	    idirs[i] = iperm[(i+swp_num)%dim];
	    m[i] = vst->m[i];
	    vel_field[i] = vst->vel_field[i];	
	    Q[i][idirs[i]] = 1.0;
	}
	for (; i < 3; ++i)
	    Q[i][i] = 1.0;
	rho = vst->rho;
	state = vst->state;
	coords = vst->coords;
	en_den = vst->en_den;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	vacuum_dens = vst->vacuum_dens;
	min_pressure = vst->min_pressure;
	min_energy = vst->min_energy;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	//hard wird for inflow
	for(i = 0; i < npts; ++i)
	{
	    float	jetradius=0.1, jetcenter[] = {0.,0.,4.};
//	    float	jetradius=0.1, jetcenter[] = {0.,0.,3.9};	    
	    coords[i] = Coords(sten->pstore[i]);
	    if(sten->newcomp == 2 && 
	       coords[i][2] > jetcenter[2] &&
	       ( sqr(coords[i][0]-jetcenter[0])
	        +sqr(coords[i][1]-jetcenter[1] ) < sqr(jetradius)))
	    {
		ft_assign(inst, sten->ststore[i], fr->sizest);
		get_constant_state(inst,sten->newcomp,coords[i],fr->time);
		for(j=i; j<npts; j++)
		    sten->ststore[j] = inst;
		break;
	    }
	}
	

	//WENO_xiaoxue
	/*The old point FD used with MUSCL scheme need 5 point stencil, The new point_FD used by weno scheme need 9 extra point*/
	    //Hard coding code to assign the state for point_FD.
	    icoordsp[0] = icoords[0][0];
	    icoordsp[1] = icoords[0][1];
	    icoordsp[2] = icoords[0][2];
	    coordsx = Rect_coords(icoordsp ,wave); //coodinate of the point in the middle of the stencil
	    int sten_comp[5]; 
	    comp = Rect_comp(icoordsp,wave);
	    for( i = 0; i < 5; i++)
	    {
		icoordsp[idirs[0]] = icoords[i-2][idirs[0]]; // i - 2 = -2,-1, 0, 1, 2
		sten_comp[i] = Rect_comp(icoordsp,wave);
	    }

	    if( sten_comp[3] == sten_comp[4] && sten_comp[3] + sten_comp[4] == 2*3)
		upboundary = 0;
	    else if( sten_comp[0] == sten_comp[1] && sten_comp[0] + sten_comp[1] == 2*3)
		upboundary = 1;	    
	    else if( sten_comp[3] == 3 && sten_comp[0] != 3 && sten_comp[1] != 3 && sten_comp[4] != 3)
		upboundary = 0;
	    else if( sten_comp[1] == 3 && sten_comp[0] != 3 && sten_comp[3] != 3 && sten_comp[4] != 3)
		upboundary = 1;
	    else
		printf("ERROR in the componet of the point_FD! comp = [%d  %d  %d  %d  %d]\n", sten_comp[0], sten_comp[1], sten_comp[2], sten_comp[3], sten_comp[4]);
		
	    int pdir;
	    if(upboundary) pdir = -1; else pdir = 1;
	    comp = Rect_comp(icoords[0],wave);
	    for( i = 0; ;i++)
	    {
		icoordsp[idirs[0]] = icoords[0][idirs[0]] + i*pdir;
		if(comp != Rect_comp(icoordsp,wave) || icoordsp[idirs[0]] < imin[idirs[0]] || icoordsp[idirs[0]] >= imax[idirs[0]])
			break;
	    }
	    leng = 2 + i;

	if(upboundary)
	{
	    for(i = 0; i < 7; i++)
	    {
	        int ldex;
		if(leng < 5)
		{
			j = 0;
			int first_interior_cell = 0;
			for(j = 0; j < 5; j++)
			{
			    if(sten_comp[j] == 3)
				break;
			    first_interior_cell++;
			}			
		    	state[i] = st = sten->ststore[first_interior_cell];
	    		coords[i] = Coords(sten->pstore[first_interior_cell]);//TO BE FIXED
	    		rho[i] = Dens(st);
	    		en_den[i] = Energy(st);
	    		for (j = 0; j < dim; ++j)
	    			m[j][i] = Mom(st)[idirs[j]];
	    		for (j = 0; j < dim; ++j)
	    			vel_field[j][i] = Vel_Field(st)[idirs[j]];

	    		if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            		{
            	    		/* Include partial density into vst */
            	    		Gas_param *params = Params(st);
            	    		int    num_comps;
            	    		float  *prho;
            	    		float  **rho0 = vst->rho0;
            	    		if((num_comps = params->n_comps) != 1)
            	    		{
            	        		prho = pdens(st);
            	        		for(j = 0; j < num_comps; j++)
           	           		rho0[j][i] = prho[j];
            	    		}
            		}
			#if !defined(UNRESTRICTED_THERMODYNAMICS)
	    		vacuum_dens[i] = Vacuum_dens(st);
	    		min_pressure[i] = Min_pressure(st);
	    		min_energy[i] = Min_energy(st);
			#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */  
		}
		else
		{
	        	if(i >= 9 - (leng - 5)) ldex = icoords[0][idirs[0]] - 11 + i;
		        else ldex = icoords[0][idirs[0]] + 3 - leng;
			icoordsp[idirs[0]] = ldex;
			state[i] = st = Rect_state(icoordsp,wave);
			coords[i] = Rect_coords(icoordsp,wave);
	   		rho[i] = Dens(st);
	    		en_den[i] = Energy(st);
	    		for (j = 0; j < dim; ++j)
	    		m[j][i] = Mom(st)[idirs[j]];
	    		for (j = 0; j < dim; ++j)
	    		vel_field[j][i] = Vel_Field(st)[idirs[j]];

	    		if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            		{
                		/* Include partial density into vst */
                		Gas_param *params = Params(st);
                		int    num_comps;
                		float  *prho;
                		float  **rho0 = vst->rho0;
                		if((num_comps = params->n_comps) != 1)
                		{
               	     	    	    prho = pdens(st);
                	    	    for(j = 0; j < num_comps; j++)
                	        	rho0[j][i] = prho[j];
                		}
            		}
			#if !defined(UNRESTRICTED_THERMODYNAMICS)
	    		vacuum_dens[i] = Vacuum_dens(st);
	    		min_pressure[i] = Min_pressure(st);
	    		min_energy[i] = Min_energy(st);
			#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	     	}
	     }
	     for (i = 7; i < 7 + npts; ++i)
	     {
	    	state[i] = st = sten->ststore[i-7];
	    	coords[i] = Coords(sten->pstore[i-7]);
	    	rho[i] = Dens(st);
	    	en_den[i] = Energy(st);
	    	for (j = 0; j < dim; ++j)
	    		m[j][i] = Mom(st)[idirs[j]];
	    	for (j = 0; j < dim; ++j)
	    		vel_field[j][i] = Vel_Field(st)[idirs[j]];

	    	if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            	{
            	    /* Include partial density into vst */
            	    Gas_param *params = Params(st);
            	    int    num_comps;
            	    float  *prho;
            	    float  **rho0 = vst->rho0;
            	    if((num_comps = params->n_comps) != 1)
            	    {
            	        prho = pdens(st);
            	        for(j = 0; j < num_comps; j++)
           	           rho0[j][i] = prho[j];
            	    }
            	}
		#if !defined(UNRESTRICTED_THERMODYNAMICS)
	    	vacuum_dens[i] = Vacuum_dens(st);
	    	min_pressure[i] = Min_pressure(st);
	    	min_energy[i] = Min_energy(st);
		#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	     }
	}
	else
	{
	    for (i = 0; i < npts; ++i)
	    {
	    	state[i] = st = sten->ststore[i];
	    	coords[i] = Coords(sten->pstore[i]);
	    	rho[i] = Dens(st);
	    	en_den[i] = Energy(st);
	    	for (j = 0; j < dim; ++j)
	    		m[j][i] = Mom(st)[idirs[j]];
	    	for (j = 0; j < dim; ++j)
	    		vel_field[j][i] = Vel_Field(st)[idirs[j]];

	    	if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            	{
            	    /* Include partial density into vst */
            	    Gas_param *params = Params(st);
            	    int    num_comps;
            	    float  *prho;
            	    float  **rho0 = vst->rho0;
            	    if((num_comps = params->n_comps) != 1)
            	    {
            	        prho = pdens(st);
            	        for(j = 0; j < num_comps; j++)
            	            rho0[j][i] = prho[j];
            	    }
            	}
		#if !defined(UNRESTRICTED_THERMODYNAMICS)
	    	vacuum_dens[i] = Vacuum_dens(st);
	    	min_pressure[i] = Min_pressure(st);
	    	min_energy[i] = Min_energy(st);
		#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	    }
	    for(i = npts; i < npts + 7; i++)
	    {
	        int ldex;
	        if(i < leng) ldex = icoords[0][idirs[0]] + i - 2 ; else ldex = icoords[0][idirs[0]] + leng - 3;
	        icoordsp[idirs[0]] = ldex;
		state[i] = st = Rect_state(icoordsp,wave);
		coords[i] = Rect_coords(icoordsp,wave);
	    	rho[i] = Dens(st);
	    	en_den[i] = Energy(st);
	    	for (j = 0; j < dim; ++j)
	    	    m[j][i] = Mom(st)[idirs[j]];
	    	for (j = 0; j < dim; ++j)
	    	    vel_field[j][i] = Vel_Field(st)[idirs[j]];

	    	if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            	{
                    /* Include partial density into vst */
                    Gas_param *params = Params(st);
                    int    num_comps;
                    float  *prho;
                    float  **rho0 = vst->rho0;
                    if((num_comps = params->n_comps) != 1)
                    {
                    	prho = pdens(st);
                    	for(j = 0; j < num_comps; j++)
                        rho0[j][i] = prho[j];
                    }
            	}
		#if !defined(UNRESTRICTED_THERMODYNAMICS)
	    	vacuum_dens[i] = Vacuum_dens(st);
	    	min_pressure[i] = Min_pressure(st);
	    	min_energy[i] = Min_energy(st);
		#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */		
	    }		
	}

	Vec_Gas_field_set(vst,state) = YES;
	Vec_Gas_field_set(vst,coords) = YES;
	Vec_Gas_field_set(vst,rho) = YES;
	Vec_Gas_field_set(vst,en_den) = YES;
	Vec_Gas_field_set(vst,m) = YES;
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            if(Params(sten->ststore[0])->n_comps != 1)
                Vec_Gas_field_set(vst,rho0) = YES;
        }

#if !defined(UNRESTRICTED_THERMODYNAMICS)
	Vec_Gas_field_set(vst,vacuum_dens) = YES;
	Vec_Gas_field_set(vst,min_pressure) = YES;
	Vec_Gas_field_set(vst,min_energy) = YES;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	set_params_jumps(vst,0,npts);

	/* If using cylindrical coordinates and this is the
	 *  r-sweep include radius information for source
	 *  computation.  
	 */

	if (is_rotational_symmetry() && alpha > 0.0 && iperm[swp_num]==0)   
	{
	    float *radii = src->radii;

	    src->rmin = fabs(pos_radius(0.0,gr));
	    for (i = 0; i < npts; ++i)
	    	radii[i] = pos_radius(Coords(sten->p[0])[0]+(i-endpt)*dh,gr);
	}

	oned_interior_scheme(swp_num,iperm,sten->icoords[0],wave,newwave,
		             fr,newfr,sten,0,npts,vst,src,dt,dh,dim);

	int update_point;
	if(upboundary) update_point = 9;
	else update_point = 2;
	Dens(ans) = vst->rho[update_point];
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* Compute differences for mass fraction */
            Gas_param *params = Params(vst->state[update_point]);
            if(params->n_comps != 1)
            {
                for(j = 0; j < params->n_comps; j++)
                    pdens(ans)[j] = vst->rho0[j][update_point];
            }
        }
	Energy(ans) = vst->en_den[update_point];
	for (i = 0; i < dim; ++i)
	    Mom(ans)[idirs[i]] = vst->m[i][update_point];
	for (i = 0; i < dim; ++i)
	    Vel_Field(ans)[idirs[i]] = vst->vel_field[i][update_point];

	/* March 2012 Merged cell at the ramp.*/
	Locstate state_old = vst->state[update_point];
	
	for( j = 0; j < 4; j++)
		Boundary_flux(ans)[j] = Boundary_flux(state_old)[j];
	int k;
	Merge_cell(ans) = Merge_cell(state_old);
	Reflect_wall(ans) = Reflect_wall(state_old); Ramp_cellx(ans) = Ramp_cellx(state_old);Ramp_cellz(ans) = Ramp_cellz(state_old); //Nozzle(ans) = Nozzle(state_old);	
	flux_flag(ans)  = flux_flag(state_old);

	for(i = 0; i < dim; i++)
	{
	    for( j = 0; j < dim; j++)
		Tau(ans)[i][j] = Tau(state_old)[i][j];
	}
	for( i = 0; i < dim; i++)
	    Qh(ans)[i] = Qh(state_old)[i];
	for( i = 0; i < dim; i++)
	{
	    for( j = 0; j < 10; j++)
	    Qc(ans)[i][j] = Qc(state_old)[i][j];
	}

        if (vst->wall_unit != NULL)
            Wallunit(ans) = vst->wall_unit[update_point];
	Set_params(ans,vst->state[update_point]);
	set_type_of_state(ans,GAS_STATE);

	check_change_inflow_st(ans,Coords(sten->p[0]));

	//TMP check
	if(NO && pressure(ans) > 110.0)
	{
	    verbose_print_state("ans",ans);
	    print_Stencil(sten);
	    for (i = 0; i < npts; ++i)
	    {
		(void) printf("state[%d]\n",i);
		(void) printf("icoords %d %d\n",sten->icoords[i-endpt][0],sten->icoords[i-endpt][1]);
	        //fprint_raw_gas_data(stdout,state[i],dim);
	    }
	    for (i = 0; i < npts; ++i)
	    {
		(void) printf("state[%d]",i);
	        verbose_print_state("",state[i]);
	    }
	    clean_up(ERROR);
	}

	smooth_vel = NO;
	if(NO &&  (fabs(vel(0,ans)) > 100.0 || fabs(vel(1,ans)) > 100.0 
		|| fabs(vel(2,ans)) > 100.0))
	{
	    int		idir;
	    float	speed;

	    printf("smooth_vel is on in point_FD\n");
	    print_general_vector("crds=", Coords(sten->p[0]), 3, "\n");
	    
	    smooth_vel = YES;

	    idir =  iperm[swp_num];
	    ft_assign(ans, sten->st[0], sten->fr->sizest);
	    speed = fabs(vel(idir,ans)) + sound_speed(ans);
	    set_max_wave_speed(idir, speed, ans, Coords(sten->p[0]),sten->wave);
	}

	if(NO)
	{
	    if(sten->nc[-1] != 0 && sten->nc[1] != 0)
		ft_assign(ans, sten->st[0], sten->fr->sizest);
	}
	if(temperature(ans) > 18500.0)
	{
	    printf("#temp  %15.8e \n", temperature(ans));
	    smooth_vel = YES;
	}

	if(is_bad_state(ans,NO,"point_FD") || smooth_vel)
	{
	    //printf("#TK point_FD in ghyperbolic.c smooth_vel=%d\n",smooth_vel);
	    LF(dh,dt,ans,dir,swp_num,iperm,NULL,sten);
	}
	reset_gamma(ans);
#if defined(CHECK_FOR_BAD_STATES)
	if (debugging("bad_state") && is_bad_state(ans,YES,"point_FD"))
	{
	    screen("ERROR in point_FD(), bad state found\n");
	    verbose_print_state("ans",ans);
	    print_Stencil(sten);
	    for (i = 0; i < npts; ++i)
	    {
		(void) printf("state[%d]\n",i);
		(void) printf("icoords %d %d\n",sten->icoords[i-endpt][0],sten->icoords[i-endpt][1]);
	        fprint_raw_gas_data(stdout,state[i],dim);
	    }
	    for (i = 0; i < npts; ++i)
	    {
		(void) printf("state[%d]",i);
	        verbose_print_state("",state[i]);
	    }
	    clean_up(ERROR);
	}
#endif /* defined(CHECK_FOR_BAD_STATES) */
      printf("leave oned_Weno\n");


}		/*end point_FD*/


/*
*			g_two_side_npt_tang():
*	
*	Calls one_side_npt_tang_solver() first on the left side and then on the
*	right side of the curve curve.
*/


EXPORT void g_two_side_npt_tang_solver(
	float		ds,
	float		dt,
	Tan_stencil	*sten,
	Locstate	ansl,
	Locstate	ansr,
	Front		*fr)
{
	INTERFACE	*intfc;
	COMPONENT	lcomp, rcomp;
	HYPER_SURF	*hs = sten->newhs;
	float		vn, v[MAXD], n[MAXD];
	int		i;
	static	size_t	sizest = 0;
	static	int	dim = 0;
        // Test code
        int             local_switch = NO;
        void (*one_side_npt_tang_solver)(float,float,Tan_stencil*,Locstate,
                        struct _Front*);
        static int      on_switch = 0, off_switch = 0;

	if (hs == NULL)
	    return;
	intfc = hs->interface;
	lcomp = negative_component(hs);
	rcomp = positive_component(hs);
	if (dim == 0)
	{
	    sizest = fr->sizest;
	    dim = intfc->dim;
	}

	if ((dim==2) && (sten->hs[0] != NULL))
	{
	    Locstate	sl, sr;

	    slsr(sten->p[0],sten->hse[0],sten->hs[0],&sl,&sr);
	    copy_state(sten->leftst[0],sl);
	    copy_state(sten->rightst[0],sr);
	}

        // Test code
        if(debugging("local_LF"))
        {
            local_switch = local_LF_npt_tang_solver_switch(ds,sten,fr);
            if(local_switch == YES)
            {
                one_side_npt_tang_solver =
                     fr->_one_side_npt_tang_solver;
                fr->_one_side_npt_tang_solver = LFoblique;
            }
        }
	switch (wave_type(hs))
	{
	case SUBDOMAIN_BOUNDARY:
	case PASSIVE_BOUNDARY:
	    g_obstacle_state(ansl,sizest);
	    g_obstacle_state(ansr,sizest);
	    break;

	/* Should this be applied in general? */
	case	NEUMANN_BOUNDARY:
	    dim = fr->rect_grid->dim;
	    normal(sten->p[0],sten->hse[0],sten->hs[0],n,fr);
	    if (is_excluded_comp(lcomp,intfc) == YES)
	    {
		sten->comp = rcomp;
		sten->states = sten->rightst;
	    	one_side_npt_tang_solver(ds,dt,sten,ansr,fr);
		if (no_slip(hs))
		{
		    float alpha = 1.0 - adherence_coeff(hs);
		    alpha_state_velocity(alpha,ansr,dim);
		}
	    	zero_normal_velocity(ansr,n,dim);
	    	g_obstacle_state(ansl,sizest);
	    }
	    else
	    {
		sten->comp = lcomp;
		sten->states = sten->leftst;
	    	one_side_npt_tang_solver(ds,dt,sten,ansl,fr);
		if (no_slip(hs))
		{
		    float alpha = 1.0 - adherence_coeff(hs);
		    alpha_state_velocity(alpha,ansl,dim);
		}
	    	zero_normal_velocity(ansl,n,dim);
	    	g_obstacle_state(ansr,sizest);
	    }
	    break;

	case	DIRICHLET_BOUNDARY:
	    if (is_excluded_comp(lcomp,intfc) == YES)
	    	g_obstacle_state(ansl,sizest);
	    else
	    {
		sten->comp = lcomp;
		sten->states = sten->leftst;
	    	one_side_npt_tang_solver(ds,dt,sten,ansl,fr);
	    }
	    if (is_excluded_comp(rcomp,intfc) == YES)
	    	g_obstacle_state(ansr,sizest);
	    else
	    {
		sten->comp = rcomp;
		sten->states = sten->rightst;
	    	one_side_npt_tang_solver(ds,dt,sten,ansr,fr);
	    }
	    break;

	default:
	    sten->comp = lcomp;
	    sten->states = sten->leftst;
	    one_side_npt_tang_solver(ds,dt,sten,ansl,fr);
	    sten->comp = rcomp;
	    sten->states = sten->rightst;
	    one_side_npt_tang_solver(ds,dt,sten,ansr,fr);
	    break;
	}
        // Test code
        if(debugging("local_LF"))
        {
            if(local_switch == YES)
                fr->_one_side_npt_tang_solver =
                    one_side_npt_tang_solver;
        }
#if defined(CHECK_FOR_BAD_STATES)
	if (debugging("bad_state") && 
	    ((is_bad_state(ansl,YES,"g_two_side_npt_tang_solver")) ||
	     (is_bad_state(ansr,YES,"g_two_side_npt_tang_solver"))))
	{
	    int	    i, nrad = sten->npts/2;
	    char    s[80];
	    screen("ERROR in g_two_side_npt_tang_solver(), bad state found\n");
	    (void) printf("ansl - ");
	    fprint_raw_gas_data(stdout,ansl,current_interface()->dim);
	    (void) printf("ansr - ");
	    fprint_raw_gas_data(stdout,ansr,current_interface()->dim);

	    (void) printf("ds = %g, dt = %g\n",ds,dt);
	    print_general_vector("dir = ",sten->dir,dim,"\n");
	    print_Tan_stencil(fr,sten);

	    for (i = -nrad; i <= nrad; ++i)
	    {
	        (void) sprintf(s,"sten->leftst[%d]\n",i);
	        verbose_print_state(s,sten->leftst[i]);
	        (void) sprintf(s,"sten->rightst[%d]\n",i);
	        verbose_print_state(s,sten->rightst[i]);
	    }
	    verbose_print_state("ansl",ansl);
	    verbose_print_state("ansr",ansr);
	    (void) printf("Input hypersurface\n");
	    print_hypersurface(hs);
	    print_interface(intfc);
	    clean_up(ERROR);
	}
#endif /* defined(CHECK_FOR_BAD_STATES) */
}		/*end g_two_side_npt_tang_solver*/


/*
*			check_ans():
*
*	Used in this file and ggodunov.c only.
*	Returns YES if states are physical,  NO otherwise.
*/

EXPORT	bool check_ans(
	const char	*function,
	float		ds,
	float		dt,
	Locstate	ans,
	COMPONENT	comp,
	Stencil		*sten,
	int		increment_count)
{
	bool		bad;
	int		nrad = sten->npts/2;
	int		*icoords = sten->icoords[0];
	int		i, dim = sten->fr->interf->dim;

	bad = is_bad_state(ans,YES,function);
	if (bad == NO)
	{
	    for (i = -nrad; (bad==NO) && i <= nrad; ++i)
	    	bad = is_bad_state(sten->st[i],YES,function);
	}
	if (bad == YES)
	{
	    if (debugging("bad_gas"))
	    {
 	        screen("WARNING - check_ans() detects bad gas state in %s\n",
		       function);
	        print_general_vector("Position = ",Coords(sten->p[0]),dim,"");
		print_int_vector(", icoords = ",icoords,dim,"\n");
		(void) printf("ds = %g, dt = %g, comp = %d\n",ds,dt,comp);
		print_Stencil(sten);
		for (i = -nrad; i <= nrad; ++i)
		{
		    (void) printf("state[%d]\n",i);
		    (*sten->fr->print_state)(sten->st[i]);
		}
		(void) printf("ANSWER\n");
		(*sten->fr->print_state)(ans);
		(void) printf("\n");
	    }

	    if (too_many_bad_gas_states("check_ans()",increment_count,sten->fr))
		clean_up(ERROR);
	    return NO;
	}
	return YES;
}		/*end check_ans*/

/*
*			check_gas():
*
*	Returns YES if states are physical,  NO otherwise.
*/

EXPORT bool check_gas(
	const char      *function,
	Locstate	*sts,
	Locstate	ans,
	Tan_stencil	*sten,
	int		increment_count,
	Front		*fr)
{
	bool		bad;
	int		i, j, nrad = sten->npts/2;
	int		dim = fr->rect_grid->dim;
	char		title[20];

	bad = is_bad_state(ans,YES,function);
	if (bad == NO)
	{
	    for (i = -nrad; (bad==NO) && i <= nrad; ++i)
	    	bad = is_bad_state(sts[i],YES,function);
	}
	if (bad == YES)
	{
	    if (debugging("bad_gas"))
	    {
	        (void) printf("WARNING in check_gas(),  bad gas state\n");
	        for (i = -nrad; i <= nrad; ++i)
	        {
	            (void) printf("curve[%d] = %llu, point[%d] = %llu",i,
	    		          curve_number(Curve_of_hs(sten->hs[i])),i,
	                          point_number(sten->p[i]));
	            if (sten->p[i] != NULL)
	            {
	                for (j = 0; j < dim; ++j)
	                    (void) printf(" %g",Coords(sten->p[i])[j]);
	            }
	            (void) printf("\n");
	            (void) sprintf(title,"state[%d]",i);
	            verbose_print_state(title,sts[i]);
	        }
	        verbose_print_state("NEW",ans);
	        (void) printf("\n");
	    }

	    if (too_many_bad_gas_states("check_gas()",increment_count,fr))
	    	clean_up(ERROR);
		
	    return NO;
	}
	return YES;
}		/*end check_gas*/


LOCAL	bool too_many_bad_gas_states(
	const char	*mesg,
	int		increment_count,
	Front		*front)
{
	static const  int MAX_WARNINGS = 20;/*TOLERANCE*/
	static int total_warnings = 0;
	static int warnings_per_step = 0;
	static int timestep_of_last_warning = -1;

	if (!increment_count)
	    return NO;
	++total_warnings;
	if (front->step != timestep_of_last_warning)
	{
	    timestep_of_last_warning = front->step;
	    warnings_per_step = 0;
	}
	if (++warnings_per_step > MAX_WARNINGS)
	{
	    screen("Fatal ERROR in %s\nERROR - too many (%d) bad gas states "
	           "in a single time step\n",
		   mesg,warnings_per_step);
	    screen("ERROR - Total warnings = %d\n",total_warnings);
	    return YES;
	}
	return NO;
}		/*end too_many_bad_gas_states*/

/*
*			g_load_state_vectors():
*
*	This function loads the conservative state variables from the wave into
*	a Vec_Gas suitable for use by the vector solvers.  This function is
*	currently used only by the vector solvers, thus we make the assumption
*	that imin == 0.	 This means there are offset artificial states 
*	on either end of the Vec_Gas -- 0 to (offset-1) and (imax-offset) to
*	(imax-1).  
*
*	To simply load a Vec_Gas from arbitray imin to arbitray imax,
*	just pass offset = 0.
*/

EXPORT TIME_DYNAMICS g_load_state_vectors(
	int		swp_num,
	int		*iperm,
	Vec_Gas		*vst,
	int		imin,
	int		imax,
	Wave		*wv,
	Wave		*newwv,
	int		*icoords,
	int		pbuf)
{
	Locstate	rst;
	Locstate	*state = vst->state;
	float		**coords = vst->coords;
	float		*rho = vst->rho;
	float		*en_den = vst->en_den;
	float		*wall_unit = vst->wall_unit;
	float		*m[MAXD];
	float		*vel_field[MAXD];
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	float		*vacuum_dens = vst->vacuum_dens;
	float		*min_pressure = vst->min_pressure;
	float		*min_energy = vst->min_energy;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	int		i, j;
	int		dim = wv->rect_grid->dim;
	int		idirs[MAXD];
	static float    **Q = NULL;
	static Locstate	inst = NULL;

	clear_Vec_Gas_set_flags(vst);
	for (j = 0; j < dim; ++j)
	{
	    idirs[j] = iperm[(j+swp_num)%dim];
	    m[j] = vst->m[j];
	    vel_field[j] = vst->vel_field[j];
	}
	if (Q == NULL)
	    bi_array(&Q,3,3,FLOAT);
	
	if(inst == NULL)
	    scalar(&inst, wv->sizest);

	Q[0][0] = Q[0][1] = Q[0][2] = 0.0;
	Q[1][0] = Q[1][1] = Q[1][2] = 0.0;
	Q[2][0] = Q[2][1] = Q[2][2] = 0.0;
	for (i = 0; i < dim; ++i)
	    Q[i][idirs[i]] = 1.0;
	for (; i < 3; ++i)
	    Q[i][i] = 1.0;
	vst->Q = (const float* const*)Q;
	for (i = imin; i < imax; ++i)
	{
	    icoords[idirs[0]] = i + pbuf;
	    coords[i] = Rect_coords(icoords,wv);
	    
	    rst = Rect_state(icoords,wv);
	    
	    //hard wird inflow
	    {
	    float	jetradius=0.1, jetcenter[] = {0.,0.,4.};
	    
	    if(Rect_comp(icoords,wv) == 2 && 
	       coords[i][2] > jetcenter[2] &&
	       sqr(coords[i][0]-jetcenter[0])
	      +sqr(coords[i][1]-jetcenter[1]) < sqr(jetradius))
	    {
		ft_assign(inst,Rect_state(icoords,wv),wv->sizest);
		get_constant_state(inst,2,coords[i],wv->time);
		rst = inst;
	    }
	    }

	    wall_unit[i] = Wallunit(rst);
	    state[i] = rst;
	    rho[i] = Dens(rst);
	    en_den[i] = Energy(rst);
	    for (j = 0; j < dim; ++j)
		m[j][i] = Mom(rst)[idirs[j]];
	    for (j = 0; j < dim; ++j)
		vel_field[j][i] = Vel_Field(rst)[idirs[j]];

            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
                /* Include partial density into vst */
                Gas_param *params = Params(rst);
                int    num_comps;
                float  *prho;
                float  **rho0 = vst->rho0;
                if((params != NULL) &&
                   ((num_comps = params->n_comps) != 1))
                {
                    prho = pdens(rst);
                    for(j = 0; j < num_comps; j++)
                        rho0[j][i] = prho[j];
                }
            }
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	    if (Params(rst) != NULL)
	    {
	    	vacuum_dens[i] = Vacuum_dens(rst);
	    	min_pressure[i] = Min_pressure(rst);
	    	min_energy[i] = Min_energy(rst);
	    }
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	}
	Vec_Gas_field_set(vst,state) = YES;
	Vec_Gas_field_set(vst,coords) = YES;
	Vec_Gas_field_set(vst,rho) = YES;
	Vec_Gas_field_set(vst,en_den) = YES;
	Vec_Gas_field_set(vst,wall_unit) = YES;
	Vec_Gas_field_set(vst,m) = YES;
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            float  **rho0 = vst->rho0;
            if((Params(state[imin]) != NULL) &&
               (Params(state[imin])->n_comps != 1))
                Vec_Gas_field_set(vst,rho0) = YES;
        }
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	Vec_Gas_field_set(vst,min_pressure) = YES;
	Vec_Gas_field_set(vst,min_energy) = YES;
	Vec_Gas_field_set(vst,vacuum_dens) = YES;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	set_params_jumps(vst,imin,imax);
	if ((Params(state[0]) == NULL) && (vst->nprms <= 1))
	{
	    int	nrad = vsten_radius(wv);
	    copy_states_to_new_time_level(idirs,wv,newwv,imin+nrad,
				          imax-nrad,icoords,pbuf);
	    return CONSTANT_IN_TIME;
	}
	g_init_obstacle_states(vst,imin,dim);
	return DYNAMIC_IN_TIME;
}		/*end g_load_state_vectors*/

/*
*			g_init_obstacle_states():
*
*	This function copies a valid state into any index which is an
*	obstacle state.	 This allows the interior solver to ignore 
*	obstacle states completely.  See assign_wave_state_vectors().
*/

LOCAL void g_init_obstacle_states(
	Vec_Gas		*vst,
	int		imin,
	int		dim)
{
	int		start, end;     /*indices of obstacle state*/
	int		non_obst;     	/*index of state to be copied*/
	int		*prms_jmp = vst->prms_jmp;
	int		i, j, k;

	for (i = 0; i < vst->nprms; ++i)
	{
	    if (is_obstacle_state(vst->state[prms_jmp[i]]))
	    {
	        start = prms_jmp[i];
	        end = prms_jmp[i + 1];

	        non_obst = (start == imin) ? end : start - 1;

	        for (j = start; j < end; ++j)
	        {
	            vst->rho[j] = vst->rho[non_obst];
	            vst->en_den[j] = vst->en_den[non_obst];
	            for (k = 0; k < dim; ++k)
	            	vst->m[k][j] = vst->m[k][non_obst];
	            vst->state[j] = vst->state[non_obst];
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	            vst->vacuum_dens[j] = vst->vacuum_dens[non_obst];
	            vst->min_pressure[j] = vst->min_pressure[non_obst];
	            vst->min_energy[j] = vst->min_energy[non_obst];
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	        }
	    }
	}
}		/*end g_init_obstacle_states*/

EXPORT	void	set_rotation(
	float	    **Q,
	const float *dir,
	int	    dim)
{
	float	dmin, mag;
	int	i, imin;

	for (i = 0; i < dim; ++i)
	    Q[0][i] = dir[i];
	switch (dim)
	{
	case 1:
	    break;
	case 2:
	    Q[1][0] = -Q[0][1];
	    Q[1][1] =  Q[0][0];
	    break;
	case 3:
	    dmin = fabs(dir[dim-1]);	imin = dim-1;
	    for (i = 0; i < (dim-1); ++i)
	    {
	    	if (fabs(dir[i]) < dmin)
	    	{
	    	    dmin = fabs(dir[i]);
	    	    imin = i;
	    	}
	    }
	    Q[1][imin] = 0.0;
	    Q[1][(imin+1)%3] = -Q[0][(imin+2)%3];
	    Q[1][(imin+2)%3] =  Q[0][(imin+1)%3];
	    mag = mag_vector(Q[1],dim);
	    Q[1][0] /= mag; Q[1][1] /= mag; Q[1][2] /= mag;
	    mag = vector_product(Q[0],Q[1],Q[2],dim);
	    Q[2][0] /= mag; Q[2][1] /= mag; Q[2][2] /= mag;
	    break;
	}
}		/*end set_rotation*/

EXPORT	bool	g_detect_and_load_mix_state(
	int	idir,
	Stencil *sten,
	int	indx)
{
	COMPONENT  *comp = sten->comp;
	if(g_composition_type() != MULTI_COMP_NON_REACTIVE)
	    return NO;
	if (detect_ambient_vapor_mix(idir,sten) == YES)
	{
	    if(comp[indx] == VAPOR_COMP || comp[indx] == AMBIENT_COMP)
	    {
	        ft_assign(sten->worksp[indx],Rect_state(sten->icoords[indx],sten->wave),sten->fr->sizest);
	        sten->st[indx] = sten->worksp[indx];
	        return YES;
	    }
	}
	return NO;
}	/* end g_detect_and_load_mix_state */

LOCAL	int	detect_ambient_vapor_mix(
	int	idir,
	Stencil *sten)
{
        int     i;
        int     find_2nd = NO;
	int	endpt = stencil_radius(sten->wave);
        if(sten->newcomp == AMBIENT_COMP &&
           sten->comp[0] == VAPOR_COMP)
        {
            for(i = -endpt; i <= endpt; ++i)
            {
                if(sten->newcomp == AMBIENT_COMP && sten->comp[i] == VAPOR_COMP &&
                   Find_rect_comp(idir,sten->icoords[i],sten->newwave)
                   == AMBIENT_COMP)
                {
                    find_2nd = YES;
                    return YES;
                }
            }
        }
        else
            return NO;	
}	/* end detect_ambient_vapor_mix */

LOCAL int local_LF_npt_tang_solver_switch(
        float        ds,
        Tan_stencil  *sten,
        Front        *fr)
{
        COMPONENT  ncp, pcp;
        CURVE      *c = NULL;
        Locstate   *lsts = sten->leftst;
        Locstate   *rsts = sten->rightst;
        float      Tl, Tr;
        int        i, nrad = sten->npts/2;

        if(sten->newhs == NULL)
            return NO;
        if(wave_type(sten->newhs) <= FIRST_USER_BOUNDARY_TYPE)
            return NO;
        ncp = negative_component(sten->newhs);
        pcp = positive_component(sten->newhs);
        c = Curve_of_hs(sten->newhs);

        /* Size effact */
        if (c->num_points < 30)
        {
            on_what = 1;
            return YES;
        }

        if(wave_type(sten->newhs) == NEUMANN_BOUNDARY ||
           wave_type(sten->newhs) == DIRICHLET_BOUNDARY)
        {
            if (is_excluded_comp(ncp,fr->interf) == YES)    
            {
                /* Temperature effact */
                Tr = temperature(rsts[-nrad]);
                for (i = -nrad+1; i <= nrad; ++i)
                {
                    Tl = temperature(rsts[i]);
                    if ((Tr > Tl*1.2) ||
                        (Tl > Tr*1.2))
                    {
                        on_what = 2;
                        return YES;
                    }
                }
            }
            else
            {
                /* Temperature effact */
                Tl = temperature(lsts[-nrad]);
                for (i = -nrad+1; i <= nrad; ++i)
                {
                    Tr = temperature(lsts[i]);
                    if ((Tr > Tl*1.2) ||
                        (Tl > Tr*1.2))
                    {
                        on_what = 2;
                        return YES;
                    }
                }
            }
            return NO;
        }

        /* Size effact */
        if (c->num_points < 30)
        {
            on_what = 1;
            return YES;
        }
        /* Temperature effact */
        for (i = -nrad; i <= nrad; ++i)
        {
            Tl = temperature(lsts[i]);
            Tr = temperature(rsts[i]);
            if ((Tr > Tl*1.2) ||
                (Tl > Tr*1.2))
            {
                on_what = 2;
                return YES;
            }
        }
        /* Curvature effact */
        if(fabs(1.0/sten->curvature/ds) < 2.0)
        {
            on_what = 3;
            return YES;
        }
        return NO;
}	/* end local_LF_npt_tang_solver_switch */
