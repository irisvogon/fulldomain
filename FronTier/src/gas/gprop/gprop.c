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
*
*				gprop.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	g_point_propagate() propagates a given point on the front
*	in the normal direction.
*/

#define DEBUG_STRING "g_prop"

#include <gdecs/gdecs.h>

#define	find_propagation_direction(p,hse,hs,nor,fr)	normal(p,hse,hs,nor,fr)

	/* LOCAL Function Declarations */
LOCAL 	void 	calculate_frame_velocity(Locstate,Locstate,int,float*);
LOCAL	void	g_point_propagate(Front*,POINTER,POINT*,POINT*,
                                  HYPER_SURF_ELEMENT*,HYPER_SURF*,float,float*);
LOCAL	float	endpt_limiter(float,float,float);
LOCAL	float	van_leer_limiter(float,float,float);
LOCAL	void	SetWSStenData(WSSten*,Front*,Wave*,POINT*,
			      HYPER_SURF_ELEMENT*,HYPER_SURF*,float);
LOCAL	void	g_ws_interpolate(Locstate,float,SIDE,int,WSSten*);
LOCAL	void	set_slopes_of_states(float,Locstate*,Locstate*,WSSten*);
LOCAL	void	set_state_slopes(Locstate,Locstate*);
LOCAL	void	g_ClearWSStenData(WSSten*);
LOCAL	void	g_set_ws_slopes(WSSten*);
LOCAL	void    interface_crosses_stencil(WSSten*,int*,int*,
                                          Locstate*,Locstate*);
LOCAL	void	set_end_state_slopes(float,Locstate,Locstate,
				    Locstate,Locstate,WSSten*);
LOCAL	void	interface_crosses_stencil2d(WSSten*,float**,float**,
                                            Locstate,Locstate);
LOCAL	void	g_unsplit_point_propagate(Front*,POINTER,POINT*,POINT*, 
                                          HYPER_SURF_ELEMENT*,HYPER_SURF*,
					  float,float*);
LOCAL	void	g_set_USWSSten2d(USWSSten2d*,POINT*,HYPER_SURF_ELEMENT*,
			         HYPER_SURF*,Front*,Wave*,float);
LOCAL   bool    is_bdry_node(POINT*,CURVE*);
LOCAL	void 	compute_vcm_and_omega(Front*,POINTER,CURVE*,float);
LOCAL	void	interface_crosses_stencil3d(WSSten*,float**,float**,
                                            Locstate,Locstate);
#if defined(DEBUG_POINT_PROPAGATE)
LOCAL	void	entered_point_propagate(char*,HYPER_SURF*);
LOCAL	void	left_point_propagate(char*,POINT*,float*,float,int, WSSten*);
#else /* defined(DEBUG_POINT_PROPAGATE) */
#define	entered_point_propagate(fname,oldhs)
#define	left_point_propagate(fname,newp,V,dt,dim,wssten)
#endif /* defined(DEBUG_POINT_PROPAGATE) */


LOCAL void calculate_frame_velocity(
	Locstate	sl,
	Locstate	sr,
	int		dim,
	float		*V)
{
	float denl = Dens(sl);
	float denr = Dens(sr);
	int i;

        if(is_obstacle_state(sr))
        {
            if(!is_obstacle_state(sl))          
                for(i = 0; i < dim; ++i)
                    V[i] = -vel(i,sl);
                return;
        }
        else if(is_obstacle_state(sl))
        {
            for(i = 0; i < dim; ++i)
                V[i] = -vel(i,sr);
            return;
        }
	else
	    for(i = 0; i < dim; ++i)
		V[i] = -(denl*vel(i,sl) + denr*vel(i,sr))/
			(denl + denr);
}	/* end calculate_frame_velocity */

/*
*			g_point_propagate():
*
*	This function propagates a point on the front in the normal direction.
*
*	The front's contribution to the timestep is updated here.
*/

bool	bdry_type(SURFACE*);

bool	bdry_type(SURFACE	*s)
{
	return no_slip(Hyper_surf(s));
}

float	surf_ten(SURFACE*);

float	surf_ten(SURFACE	*s)
{
	return surface_tension(Hyper_surf(s));
}

bool	block_cross_test(
	float		*coords,
	INTERFACE	*intfc_in)
{
	INTERFACE	*intfc = prev_interface(intfc_in);
	RECT_GRID	*gr = &topological_grid(intfc);
	struct Table	*T = table_of_interface(intfc);
	int		i, j, k, m, nc, ic[3], ib[3], *gmax;
	GRID_DIRECTION	dir[3] = {EAST, NORTH, UPPER};

	gmax = gr->gmax;
	if(!rect_in_which(coords, ib, gr))
	    return NO;

	//TMP use for shear flow
	return YES;
	
	//print_general_vector("cords=", coords, 3, "\n");
	//print_int_vector("ib", ib, 3, "\n");

	for(m=0; m<3; m++)
	    for(i=0; i<2; i++)
		for(j=0; j<2; j++)
		{
		    ic[m] = ib[m];
		    ic[(m+1)%3] = ib[(m+1)%3] + i;
		    ic[(m+2)%3] = ib[(m+2)%3] + j;
		    k = seg_index3d(ic[0],ic[1],ic[2],dir[m],gmax);
		    nc = T->seg_crx_count[k];
		    if(nc > 1)
			return NO;
		}

	return YES;
}

#define	corner_index(p,i,gr)	irint(floor(((p)-(gr)->L[i])/(gr)->h[i]-0.5))
bool	find_block(
	float		*f,
	int		*icrds,
	float		*p,
	RECT_GRID	*gr)
{
	int	i;

	for(i=0; i<3; i++)
	{
	    icrds[i] = corner_index(p[i],i,gr);
	    if(icrds[i] < -gr->lbuf[i] || icrds[i] >= gr->gmax[i]+gr->ubuf[i]-1)
		return  NO;
	    f[i] = p[i] - (gr->L[i]+(0.5+icrds[i])*gr->h[i]);
	    f[i] /= gr->h[i];
	}
	return  YES;
}

bool	Peskin_sample_speed(
	float	*v,
	float	*f,
	int	*icrds,
	Wave	*wave)
{
	RECT_GRID	*gr = wave->rect_grid;
	int		i,ix,iy,iz,ic[3], imin[3], imax[3];
	float		coef[3][4], w;
	Locstate	st;

	for(i=0; i<3; i++)
	{
	    imin[i] = icrds[i] - 1;
	    if(imin[i] < -gr->lbuf[i])
	        return NO;
	    imax[i] = icrds[i] + 2;
	    if(imax[i] > gr->gmax[i] + gr->ubuf[i] - 1)
		return NO;
	}
	for(i=0; i<3; i++)
	{
	    coef[i][0] = 1.0 + cos(PI*(1.0 + f[i])/2.0);
	    coef[i][1] = 1.0 + cos(PI*f[i]/2.0);
	    coef[i][2] = 1.0 + cos(PI*(1.0 - f[i])/2.0);
	    coef[i][3] = 1.0 + cos(PI*(2.0 - f[i])/2.0);
	    v[i] = 0.0;
	}

	for(iz=0; iz<4; iz++)
	    for(iy=0; iy<4; iy++)
		for(ix=0; ix<4; ix++)
		{
		    ic[0] = ix + imin[0];
		    ic[1] = iy + imin[1];
		    ic[2] = iz + imin[2];

		    st = Rect_state(ic,wave);
		    w = coef[0][ix]*coef[1][iy]*coef[2][iz];

		    for(i=0; i<3; i++)
			v[i] += vel(i,st)*w;
		}
	for(i=0; i<3; i++)
	    v[i] /= 64.0;

	//print_general_vector("pkv", v, 3, "\n");

	return YES;
}

bool	sample_speed(
	float	*v,
	float	*f,
	int	*icrds,
	Wave	*wave)
{
	int		i, j, k, ic[3];
	Locstate	s1, s2;
	float		velz[2][2][3], dens[2][2], den;
	float		ncor, v1[3], v2[3];

	if(debugging("Pessam"))
	    if(Peskin_sample_speed(v, f, icrds, wave))
		return YES;
	
	for(i=0; i<2; i++)
	    for(j=0; j<2; j++)
	    {
		ic[0] = icrds[0] + i;
		ic[1] = icrds[1] + j;
		ic[2] = icrds[2];
		s1 = Rect_state(ic,wave);	//bottom state
		//print_general_vector("p ", Rect_coords(ic,wave), 3, "\n");
if(NO && debugging("flgasexp"))
{
		if(Rect_comp(ic,wave) == 2)
		    ft_assign(v1, Gvel(s1)[0], 3*FLOAT);
		else
		{
		    ncor = 0.0;
		    for(k=0; k<3; k++)
			ncor += (Gvel(s1)[1][k] - Gvel(s1)[0][k])*Gn(s1)[k];
		    for(k=0; k<3; k++)
			v1[k] = Gvel(s1)[0][k] + ncor*Gn(s1)[k];
		}
}	
		ic[2]++;
		s2 = Rect_state(ic,wave);	//top state
		//print_general_vector("p ", Rect_coords(ic,wave), 3, "\n");
if(NO && debugging("flgasexp"))
{
		if(Rect_comp(ic,wave) == 2)
		    ft_assign(v2, Gvel(s2)[0], 3*FLOAT);
		else
		{
		    ncor = 0.0;
		    for(k=0; k<3; k++)
			ncor += (Gvel(s2)[1][k] - Gvel(s2)[0][k])*Gn(s2)[k];
		    for(k=0; k<3; k++)
			v2[k] = Gvel(s2)[0][k] + ncor*Gn(s2)[k];
		}
}
	
if(debugging("flgasexp"))
{
if(debugging("momvel"))
{
		for(k=0; k<3; k++)
		    velz[i][j][k] = vel(k,s1)*Dens(s1)*(1.0-f[2]) + vel(k,s2)*Dens(s2)*f[2];
		dens[i][j] = Dens(s1)*(1.0-f[2]) + Dens(s2)*f[2];
}
else
{
		//use the fluid velocity as the moving velocity.
		for(k=0; k<3; k++)
		{
		    //velz[i][j][k] = vel(k,s1)*(1.0-f[2]) + vel(k,s2)*f[2];
		    //previous fluid
		    velz[i][j][k] = Gvel(s1)[0][k]*(1.0-f[2]) + Gvel(s2)[0][k]*f[2];
		    //velz[i][j][k] = v1[k]*(1.0-f[2]) + v2[k]*f[2];
		}
}
}
else
{
		for(k=0; k<3; k++)
		    velz[i][j][k] = vel(k,s1)*(1.0-f[2]) + vel(k,s2)*f[2];
}
		if(debugging("dbspd"))
		    print_general_vector("velz", velz[i][j], 3, "\n");

	    }  //for j

if(debugging("momvel"))
{
	dens[0][0] = dens[0][0]*(1.0-f[1]) + dens[0][1]*f[1];
	dens[1][0] = dens[1][0]*(1.0-f[1]) + dens[1][1]*f[1];
	den = dens[0][0]*(1.0-f[0]) + dens[1][0]*f[0];

	for(k=0; k<3; k++)
	{
	    velz[0][0][k] = velz[0][0][k]*(1.0-f[1]) + velz[0][1][k]*f[1];
	    velz[1][0][k] = velz[1][0][k]*(1.0-f[1]) + velz[1][1][k]*f[1];
	    v[k] = (velz[0][0][k]*(1.0-f[0]) + velz[1][0][k]*f[0])/den;
	}
}
else
{
	for(k=0; k<3; k++)
	{
	    velz[0][0][k] = velz[0][0][k]*(1.0-f[1]) + velz[0][1][k]*f[1];
	    velz[1][0][k] = velz[1][0][k]*(1.0-f[1]) + velz[1][1][k]*f[1];
	    v[k] = velz[0][0][k]*(1.0-f[0]) + velz[1][0][k]*f[0];
	}
}

	if(debugging("dbspd"))
	    print_general_vector("sam v", v, 3, "\n");
	
	return YES;
}

bool	sample_state(
	Locstate	ans,
	float		*v,
	int		ind,
	float		*f,
	int		*icrds,
	Wave		*wave)
{
	int		i, j, k, ic[3],iii;
	Locstate	s1, s2;
	float		presz[2][2], densz[2][2], pdensz[2][2][10], sum;
	int		st_type;

	st_type = state_type(ans);
	set_state(ans, TGAS_STATE, ans);

	for(i=0; i<2; i++)
	    for(j=0; j<2; j++)
	    {
		ic[0] = icrds[0] + i;
		ic[1] = icrds[1] + j;
		ic[2] = icrds[2];
		s1 = Rect_state(ic,wave);	//bottom state
		//print_general_vector("p ", Rect_coords(ic,wave), 3, "\n");
		
		ic[2]++;
		s2 = Rect_state(ic,wave);	//top state
		//print_general_vector("p ", Rect_coords(ic,wave), 3, "\n");
	
		presz[i][j] = Gp(s1)[ind]*(1.0-f[2]) + Gp(s2)[ind]*f[2];
		densz[i][j] = Gd(s1)[ind]*(1.0-f[2]) + Gd(s2)[ind]*f[2];
	
		if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                {
                    for(iii = 0; iii < Params(ans)->n_comps; iii++)
                        pdensz[i][j][iii] = Gpd(s1)[ind][iii]*(1.0-f[2]) + Gpd(s2)[ind][iii]*f[2];
                }
	
		//if(debugging("dbspd"))
		//    print_general_vector("velz", velz[i][j], 3, "\n");
	    }  //for j

	presz[0][0] = presz[0][0]*(1.0-f[1]) + presz[0][1]*f[1];
	presz[1][0] = presz[1][0]*(1.0-f[1]) + presz[1][1]*f[1];
	Press(ans) = presz[0][0]*(1.0-f[0]) + presz[1][0]*f[0];
	
	densz[0][0] = densz[0][0]*(1.0-f[1]) + densz[0][1]*f[1];
	densz[1][0] = densz[1][0]*(1.0-f[1]) + densz[1][1]*f[1];
	Dens(ans) = densz[0][0]*(1.0-f[0]) + densz[1][0]*f[0];
        
	if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            for(iii = 0; iii < Params(ans)->n_comps; iii++)
            {
                pdensz[0][0][iii] = pdensz[0][0][iii]*(1.0-f[1]) + pdensz[0][1][iii]*f[1];
                pdensz[1][0][iii] = pdensz[1][0][iii]*(1.0-f[1]) + pdensz[1][1][iii]*f[1];
            }
        }
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            for(iii = 0; iii < Params(ans)->n_comps; iii++)
                pdens(ans)[iii] = pdensz[0][0][iii]*(1.0-f[0]) + pdensz[1][0][iii]*f[0];
		
	    //normalize
	    sum = 0.0;
	    for(i = 0; i < Params(ans)->n_comps; i++)
		sum += pdens(ans)[i];
	    for(i = 0; i < Params(ans)->n_comps; i++)
		pdens(ans)[i] = Dens(ans)/sum*pdens(ans)[i];
        }

	Vel(ans)[0] = v[0];
	Vel(ans)[1] = v[1];
	Vel(ans)[2] = v[2];

	set_state(ans, st_type, ans);

	//if(debugging("dbspd"))
	//    print_general_vector("sam v", v, 3, "\n");
}

LOCAL	bool sample_point_state(
	Locstate	   sl,
	Locstate	   sr,
	Wave		   *wave,
	POINT		   *oldp,
	HYPER_SURF	   *oldhs,
	float		   dt,
	float		   *V,
	INTERFACE	   *intfc)
{
	RECT_GRID	*gr = wave->rect_grid;
	float		f[3], v[3], ptmp[3];
	int		k, icrds[3], w_type;
	Locstate	st;

	w_type = wave_type(oldhs);
	if(w_type != CONTACT)
	    return NO;

	remove_from_debug("dbspd");
	if(the_point(oldp))
	{
	    add_to_debug("dbspd");
	    print_general_vector("dbspd", Coords(oldp), 3, "\n");
	}
	//printf("#contact enter.\n");
	//print_general_vector("oldp", Coords(oldp), 3, "\n");
	
	if(!find_block(f,icrds,Coords(oldp),gr) || 
	   !block_cross_test(Coords(oldp), intfc))
	    return NO;
	
	//sample_speed(v, f, icrds, wave);
	//first order
	
	if(negative_component(oldhs) == 2)
	    sample_state(sl,V,0,f,icrds,wave);
	else
	    sample_state(sl,V,1,f,icrds,wave);
	
	if(positive_component(oldhs) == 3)
	    sample_state(sr,V,1,f,icrds,wave);
	else
	    sample_state(sr,V,0,f,icrds,wave);

	return YES;
}




LOCAL	bool g_point_speed(
	Wave		   *wave,
	POINT		   *oldp,
	HYPER_SURF	   *oldhs,
	float		   dt,
	float		   *V,
	INTERFACE	   *intfc)
{
	RECT_GRID	*gr = wave->rect_grid;
	float		f[3], v[3], ptmp[3];
	int		k, icrds[3], w_type;
	Locstate	st;

	w_type = wave_type(oldhs);
	if(w_type != CONTACT)
	    return NO;

	remove_from_debug("dbspd");
	if(the_point(oldp))
	{
	    add_to_debug("dbspd");
	    print_general_vector("dbspd", Coords(oldp), 3, "\n");
	}
	//printf("#contact enter.\n");
	//print_general_vector("oldp", Coords(oldp), 3, "\n");
	
	if(!find_block(f,icrds,Coords(oldp),gr) || 
	   !block_cross_test(Coords(oldp), intfc))
	    return NO;
	
	//sample_speed(v, f, icrds, wave);
	//first order
	sample_speed(V, f, icrds, wave);
	
	return YES;
	
	if(negative_component(oldhs) == 2)
	    sample_state(left_state(oldp),V,0,f,icrds,wave);
	else
	    sample_state(left_state(oldp),V,1,f,icrds,wave);
	
	if(positive_component(oldhs) == 3)
	    sample_state(right_state(oldp),V,1,f,icrds,wave);
	else
	    sample_state(right_state(oldp),V,0,f,icrds,wave);

	//print_general_vector("V=", V, 3, "\n");

	return YES;

	//print_general_vector("f= ", f, 3, "\n");
	//print_general_vector("v= ", v, 3, "\n\n");
	    
	for(k=0; k<3; k++)
	    ptmp[k] = Coords(oldp)[k] + 0.5*dt*v[k];
	    
	if(!find_block(f,icrds,ptmp,gr) ||
	   !block_cross_test(ptmp, intfc))
	    return NO;
	sample_speed(V, f, icrds, wave);
	  
	//print_general_vector("V= ", V, 3, "\n");
	return YES;
}


LOCAL	bool g_point_speed_prev(
	Wave		   *wave,
	POINT		   *oldp,
	HYPER_SURF	   *oldhs,
	float		   dt,
	float		   *V,
	INTERFACE	   *intfc)
{
	RECT_GRID	*gr = wave->rect_grid;
	float		f[3], v[3], ptmp[3];
	int		k, icrds[3], w_type;
	Locstate	st;

	w_type = wave_type(oldhs);
	if(w_type != CONTACT)
	    return NO;

	//printf("#contact enter.\n");
	//print_general_vector("oldp", Coords(oldp), 3, "\n");
	
	if(!find_block(f,icrds,Coords(oldp),gr) || 
	   !block_cross_test(Coords(oldp), intfc))
	    return NO;
	sample_speed(v, f, icrds, wave);
	
	//print_general_vector("f= ", f, 3, "\n");
	//print_general_vector("v= ", v, 3, "\n\n");
	    
	for(k=0; k<3; k++)
	    ptmp[k] = Coords(oldp)[k] + 0.5*dt*v[k];
	    
	if(!find_block(f,icrds,ptmp,gr) ||
	   !block_cross_test(ptmp, intfc))
	    return NO;
	sample_speed(V, f, icrds, wave);
	  
	//print_general_vector("V= ", V, 3, "\n");
	return YES;
}


#if defined(DEBUG_POINT_PROPAGATE)
LOCAL	bool debug_point_propagate = NO;	/* Debugging switch */
#endif /* defined(DEBUG_POINT_PROPAGATE) */

bool	check_change_inflow_st(Locstate,float*);

LOCAL	void g_point_propagate(
	Front		   *fr,
	POINTER		   p2wave,
	POINT		   *oldp,
	POINT		   *newp,
	HYPER_SURF_ELEMENT *oldhse,
	HYPER_SURF	   *oldhs,
	float		   dt,
	float		   *V)
{
	printf("enter g_points_propagate\n");
        INTERFACE	*intfc = fr->interf;
	Wave		*wave = (Wave*)p2wave;
	int		i, dim = intfc->dim;
	int		w_type;
	float		frame_velocity[MAXD], Vtmp[3];
#if defined(DEBUG_POINT_PROPAGATE)
	static	char	fname[] = "g_point_propagate";
#endif /* defined(DEBUG_POINT_PROPAGATE) */
	static WSSten	*sten = NULL;
	bool		ptspd, movfr;

/*        printf("xiaoxue_pp\n");
	entered_point_propagate(fname,oldhs);

        //RK_PRINT("******* 1\n");

        printf("xiaoxue_pp\n");
	if (sten == NULL) 
	    sten = AllocDefaultWSSten(fr);
	else
	    ClearWSStenData(sten);

        //RK_PRINT("******* 2\n");
	printf("xiaoxue_pp\n");
        w_type = (is_subdomain_boundary(oldhs)) ?
			SUBDOMAIN_BOUNDARY : wave_type(oldhs);

	if(debugging("sample_speed"))
	    fr->movingframe = YES;

	if (w_type >= FIRST_PHYSICS_WAVE_TYPE && fr->movingframe == YES)
	{
	    sten->V = frame_velocity;
	    if (dim == 2 && is_bdry_node(oldp,Curve_of_hs(oldhs)))
	    	sten->V = NULL;
	}
	else
	    sten->V = NULL;
	
        printf("xiaoxue_pp\n");
        //RK_PRINT("******* 3\n");
	if (passive_point_propagate(w_type,oldp,newp,V,intfc)==YES)
	{
	    left_point_propagate(fname,newp,V,dt,dim,NULL);
	    return;
	}

        printf("xiaoxue_pp\n");
        //RK_PRINT("******* 4\n");
	movfr = NO;
        printf("enter SetWSStenData\n");
	SetWSStenData(sten,fr,wave,oldp,oldhse,oldhs,dt);
        printf("xiaoxue:leave SetWSStenData\n");	
	if(debugging("fcrxsten") && 
	   (sten->indxl < sten->nsts || sten->indxr < sten->nsts))
	{
	    movfr = YES;
	    remove_from_debug("sample_speed");
	}

	if(debugging("sample_speed"))
	{
	    ptspd = YES;
	    if(!g_point_speed(wave, oldp, oldhs, dt, Vtmp, intfc))
	        ptspd = NO;
	}


        printf("xiaoxue_pp\n");
        //RK_PRINT("******* 5\n");
	if(!debugging("sample_speed"))
	{
            //RK_PRINT("!!!!!!! 1\n");
	    npt_w_speed(sten,left_state(newp),right_state(newp),V);
            //RK_PRINT("!!!!!!! 1.3\n");
	    check_change_inflow_st(left_state(newp), Coords(newp));
            //RK_PRINT("!!!!!!! 1.6\n");
	    check_change_inflow_st(right_state(newp), Coords(newp));
            //RK_PRINT("!!!!!!! 2\n");
	}
	else if(w_type == CONTACT)
	{
            //RK_PRINT("!!!!!!! 3\n");
	    printf("xiaoxue_pp!\n");
            static Locstate	stl=NULL, str;
	    float		Vsurf[3], tmp[3];
	    int			st_type;
	    
            
            printf("xiaoxue_pp!\n");
            //RK_PRINT("!!!!!!! 4\n");
	    if(stl == NULL)
	    {
                printf("xiaoxue_pp!\n");
                //RK_PRINT("!!!!!!! 5\n");
		alloc_state(intfc,&stl,fr->sizest);
                //RK_PRINT("!!!!!!! 6\n");
		printf("xiaoxue_pp!\n");
                alloc_state(intfc,&str,fr->sizest);
                //RK_PRINT("!!!!!!! 7\n");
                printf("xiaoxue_pp!\n");
	    }

            //RK_PRINT("!!!!!!! 8\n");
                printf("xiaoxue_pp!\n");
	    ft_assign(left_state(newp), left_state(oldp), fr->sizest);
            //RK_PRINT("!!!!!!! 9\n");
printf("xiaoxue_pp!\n");
	    ft_assign(right_state(newp), right_state(oldp), fr->sizest);
            //RK_PRINT("!!!!!!!10\n");
printf("xiaoxue_pp!\n");	

	    ft_assign(stl, left_state(oldp), fr->sizest);
printf("xiaoxue_pp!\n"); 
           //RK_PRINT("!!!!!!!11\n");
	    ft_assign(str, right_state(oldp), fr->sizest);
  printf("xiaoxue_pp!\n"); 
         //RK_PRINT("!!!!!!!12\n");
    
	    if(!is_obstacle_state(stl) && !is_obstacle_state(str) && ptspd)
	    {
                //RK_PRINT("!!!!!!!13\n");
                printf("xiaoxue_pp!\n");
		sample_point_state(stl,str,wave,oldp,oldhs,dt,Vtmp,intfc);
                //RK_PRINT("!!!!!!!14\n");
	printf("xiaoxue_pp!\n");
		w_speed(Coords(oldp),stl,str, left_state(newp),right_state(newp),
		   Vsurf, sten->pjump, sten->nor, w_type, fr);
                //RK_PRINT("!!!!!!!15\n");
	        printf("xiaoxue_pp!\n");
		
		//left and right states have the same speed.
                //RK_PRINT("!!!!!!!16\n");
        printf("xiaoxue_pp!\n");
		for(i=0; i<3; i++)
		{
		    //Vtmp[i] += Vsurf[i];
		    //Vtmp[i] = Vsurf[i];
                printf("xiaoxue_pp!\n");
                    //RK_PRINT("!!!!!!!17 : i=%d\n",i);
		    Vtmp[i] = vel(i, left_state(newp));
		}
                printf("xiaoxue_pp!\n");
                //RK_PRINT("!!!!!!!18\n");
	    }
	}
	else
	{
	    //keep the original states for boundary surfaces 
	    //do not move the boundary surfaces
            //RK_PRINT("!!!!!!!19\n");
printf("xiaoxue_pp!!\n");
	    ft_assign(left_state(newp), left_state(oldp), fr->sizest);
            //RK_PRINT("!!!!!!!20\n");
printf("xiaoxue_pp!!\n");
	    ft_assign(right_state(newp), right_state(oldp), fr->sizest);
            //RK_PRINT("!!!!!!!21\n");
printf("xiaoxue_pp!!\n");
	    for(i=0; i<3; i++)
	    {
                //RK_PRINT("!!!!!!!22 i=%d\n",i);
printf("xiaoxue_pp!!\n");
		V[i] = 0.0;
		Vtmp[i] = 0.0;
	    }
            //RK_PRINT("!!!!!!!23\n");
printf("xiaoxue_pp!!\n");
	}
        //RK_PRINT("******* 6\n");
printf("xiaoxue_pp!!\n");

	if(fr->movingframe == YES && sten->V && !debugging("sample_speed")) 
	{
printf("xiaoxue_pp!!\n");
	    for (i = 0; i < dim; i++)
	    {
		sten->V[i] *= -1.0;
	    	V[i] += sten->V[i];
	    }
	    add_velocity_to_state(left_state(newp),sten->V);
	    add_velocity_to_state(right_state(newp),sten->V);
	}
        //RK_PRINT("******* 7\n");
printf("xiaoxue_pp!!\n");

	if (is_scalar_wave(w_type) && debugging("experiment"))
	{
	    Locstate sl =  left_state(newp);
	    Locstate sr = right_state(newp);
	    float rhol = Dens(sl), rhor = Dens(sr);
	    float A = fabs(rhol - rhor)/(rhol + rhor);

	    for (i = 0; i < dim; ++i)
	    	V[i] = (1.0-A)*V[i]+A*(mom(i,sl)+mom(i,sr))/(rhol+rhor);
	}
	
	if(debugging("sample_speed"))
	{
	    if(ptspd)
		ft_assign(V, Vtmp, 3*FLOAT);
	    else
	    {
		V[0] = 0.0;
		V[1] = 0.0;
		V[2] = 0.0;
	    }
	}
        //RK_PRINT("******* 8\n");
printf("xiaoxue_pp!!\n");
        //RK_PRINT("******* 9\n");
	for (i = 0; i < dim; ++i)
	    Coords(newp)[i] = Coords(oldp)[i] + V[i]*dt;
printf("xiaoxue_pp!!\n");
        {
	    Locstate sl, sr;
	    sl = left_state(newp);
	    sr = right_state(newp);
printf("xiaoxue_pp!!\n");
	    if(!is_obstacle_state(sl) && pressure(sl) > 306170620.0)
	    {
	        printf("sl is hight\n");
	        print_WSStenData(sten);
	    }
printf("xiaoxue_pp!!\n");  
	    if(!is_obstacle_state(sr) && pressure(sr) > 306170620.0)
	    {
	        printf("sr is hight\n");
	        print_WSStenData(sten);
	    }
	}*/
        for(int i = 0; i < dim ; i++)
            V[i] = 0.0;
	/* Record speed for time step calculation */
printf("xiaoxue_pp!!\n");
        //RK_PRINT("*******11\n");
	for (i = 0; i < dim; ++i)
	    set_max_front_speed(i,fabs(V[i]),return_obst_state(),
				Coords(newp),fr);
printf("xiaoxue_pp!!!\n");
//	set_max_front_speed(dim,fabs(scalar_product(V,sten->nor,dim))/sten->dn,
//			    return_obst_state(),Coords(newp),fr);

	set_max_front_speed(dim, 0.0,
			    return_obst_state(),Coords(newp),fr);

        if(movfr)
	    add_to_debug("sample_speed");
        
        left_point_propagate(fname,newp,V,dt,dim,sten);
}		/*end g_point_propagate*/

LOCAL	void	g_ClearWSStenData(
	WSSten	*sten)
{
	int	i, j, nsts = sten->nsts;
	size_t	sizest = g_sizest();

	sten->coords = sten->coords_store;
	sten->nor = sten->nor_store;
	for (i = 0; i < 3; ++i)
	{
	    sten->coords[i] = 0.0;
	    sten->nor[i] = 0.0;
	}

	for (i = 0; i < nsts; ++i)
	{
	    sten->lcrds[i] = sten->lcrds_store[i];
	    sten->rcrds[i] = sten->rcrds_store[i];
	    for (j = 0; j < 3; ++j)
	    {
		sten->lcrds[i][j] = 0.0;
		sten->rcrds[i][j] = 0.0;
	    }
	}
	for (i = 0; i < nsts; ++i)
	{
	    sten->sl[i] = (Locstate) (sten->sl_store + i*sizest);
	    zero_scalar(sten->sl[i],sizest);
	    set_type_of_state(sten->sl[i],UNKNOWN_STATE);

	    sten->dsl[i] = (Locstate) (sten->dsl_store + i*sizest);
	    zero_scalar(sten->dsl[i],sizest);
	    set_type_of_state(sten->dsl[i],UNKNOWN_STATE);

	    sten->tsl[i] = (Locstate) (sten->tsl_store + i*sizest);
	    zero_scalar(sten->tsl[i],sizest);
	    set_type_of_state(sten->tsl[i],UNKNOWN_STATE);

	    sten->sr[i] = (Locstate) (sten->sr_store + i*sizest);;
	    zero_scalar(sten->sr[i],sizest);
	    set_type_of_state(sten->sr[i],UNKNOWN_STATE);

	    sten->dsr[i] = (Locstate) (sten->dsr_store + i*sizest);
	    zero_scalar(sten->dsr[i],sizest);
	    set_type_of_state(sten->dsr[i],UNKNOWN_STATE);

	    sten->tsr[i] = (Locstate) (sten->tsr_store + i*sizest);
	    zero_scalar(sten->tsr[i],sizest);
	    set_type_of_state(sten->tsr[i],UNKNOWN_STATE);
	}
	sten->ncomp = NO_COMP;
	sten->pcomp = NO_COMP;
	sten->w_type = UNKNOWN_WAVE_TYPE;
	sten->hs = NULL;
	sten->dn = 0.0;
	sten->dt = 0.0;
	sten->pjump = 0.0;
	sten->front = NULL;
	sten->wave = NULL;
}		/*end g_ClearWSStenData*/

LOCAL	void	SetWSStenData(
	WSSten			*sten,
	Front			*fr,
	Wave			*wave,
	POINT			*oldp,
	HYPER_SURF_ELEMENT 	*oldhse,
	HYPER_SURF		*oldhs,
	float			dt)
{
	Locstate sl, sr;	       /* states at left/right sides of front */
	int i;

        printf("xiaoxue_setwss\n");
	ClearWSStenData(sten);/*TODO: This can be a debug line once we are
			       *sure sten is set properly
			       */
#if defined(DEBUG_POINT_PROPAGATE)
	debug_print("WSSten","Entered SetWSStenData()\n");
#endif /* defined(DEBUG_POINT_PROPAGATE) */
printf("xiaoxue_setwss\n");
	if (fr->interf->modified)
	{
	    if (!make_interface_topology_lists(fr->interf))
	    {
		screen("ERROR in SetWSStenData(), "
		       "make_interface_topology_lists(), failed\n");
		clean_up(ERROR);
	    }
	}
printf("xiaoxue_setwss\n");
	sten->w_type = wave_type(oldhs);
	sten->p = oldp;
	sten->hs = oldhs;
	sten->hse = oldhse;
	sten->front = fr;
	sten->wave = wave;
	sten->dt = dt;
printf("xiaoxue_setwss\n");
	/* Calculate normal direction */
	/* and positions of normally displaced states */
	find_propagation_direction(oldp,oldhse,oldhs,sten->nor,fr);
printf("xiaoxue_setwss\n");
	/* Calculate pressure jump due to surface tension */
	 
	sten->pjump = set_pjump_at_wave(oldp,oldhse,oldhs,fr,sten->nor);
printf("xiaoxue_setwss\n");
	/* Find states near front for normal sweep */
	slsr(oldp,oldhse,oldhs,&sl,&sr);
printf("xiaoxue_setwss\n");
	states_near_location(sten,Coords(oldp),sten->nor,
				  positive_component(oldhs),
				  negative_component(oldhs),sl,sr);
printf("xiaoxue_setwss\n");
	if (sten->V)	/* subtract moving frame velocity */
	{
	    calculate_frame_velocity(sl,sr,fr->interf->dim,sten->V);
	    for (i = 0; i < sten->nsts; ++i)
	    {
	    	add_velocity_to_state(sten->sl[i],sten->V);
		add_velocity_to_state(sten->sr[i],sten->V);
	    }
	}
printf("xiaoxue_setwss\n");
	/*TMP*/
	if (debugging("tmp_init_curve"))
	    print_WSStenData(sten);
#if defined(DEBUG_POINT_PROPAGATE)
	if (debug_point_propagate == YES) 
	    print_WSStenData(sten);
#endif /* defined(DEBUG_POINT_PROPAGATE) */

#if defined(DEBUG_POINT_PROPAGATE)
	debug_print("WSSten","Left SetWSStenData()\n");
#endif /* defined(DEBUG_POINT_PROPAGATE) */

}	/*end SetWSStenData*/

//#bjet2  
EXPORT	void states_near_location_from_WallTan(
	WSSten		*sten,
	float		*coords,
	float		*nor,
	COMPONENT	pcomp,
	COMPONENT	ncomp,
	Locstate	sl,
	Locstate	sr,
	Tan_stencil     *tan_sten,
	float		dn)
{
	Front		*fr = sten->front;
	int		dim = fr->interf->dim, sizest = fr->sizest;
	int		i, j, nsts = sten->nsts;
	

	if(nsts != tan_sten->npts/2 + 1)
	{
	    printf("ERROR, states_near_location_from_WallTan, stencils are inconsistent. %d %d\n", 
	        nsts, tan_sten->npts/2 + 1);
	    clean_up(ERROR);
	}
	
	sten->pcomp = pcomp;
	sten->ncomp = ncomp;
	if (nor != sten->nor)
	{
	    for (j = 0; j < dim; ++j)
		sten->nor[j] = nor[j];
	}
	for (j = 0; j < dim; ++j)
	    sten->coords[j] = sten->lcrds[0][j] = sten->rcrds[0][j] = coords[j];

	sten->dn = dn;
	for (i = 1; i < nsts; ++i)
	{
	    for (j = 0; j < dim; ++j)
	    {
	        sten->lcrds[i][j] = coords[j] - i*nor[j]*dn;
	        sten->rcrds[i][j] = coords[j] + i*nor[j]*dn;
	    }
	}
	set_state(sten->sl[0],state_type(sl),sl);
	set_state(sten->sr[0],state_type(sr),sr);

	//here ASSUME right states are in fluid, therefore, only rightst is used
	//also rightst[0] is not used since it will use sl, sr from fluid surface
	//states
	// rightst[-2]  rightst[-1]  sl sr rightst[1]  rightst[2]
	if (is_obstacle_state(sl))
	{
	    for (i = 0; i < nsts; ++i)
	        obstacle_state(fr->interf,sten->sl[i],sizest);
	}
	else
	{
	    for (i = 1; i < nsts; ++i)
	    {
	        ft_assign(sten->sl[i], tan_sten->rightst[-i], sizest);
		set_state(sten->sl[i],state_type(sl),sten->sl[i]);
	    }
	}

	if (is_obstacle_state(sr))
	{
	    for (i = 0; i < nsts; ++i)
	        obstacle_state(fr->interf,sten->sr[i],sizest);
	}
	else
	{
	    for (i = 1; i < nsts; ++i)
	    {
	        ft_assign(sten->sr[i], tan_sten->rightst[i], sizest);
		set_state(sten->sr[i],state_type(sr),sten->sr[i]);
	    }
	}

	if (debugging("bad_state"))
	{
	    bool bad = NO;
	    for (i = 0; i < nsts; ++i)
	    {
	        if (is_bad_state(sten->sr[i],YES,"states_near_location"))
		{
		    bad = YES;
		    (void) printf("sten->sr[%d] is bad\n",i);
	            fprint_raw_gas_data(stdout,sten->sr[i],dim);
		}
	        if (is_bad_state(sten->sl[i],YES,"states_near_location"))
		{
		    bad = YES;
		    (void) printf("sten->sl[%d] is bad\n",i);
	            fprint_raw_gas_data(stdout,sten->sl[i],dim);
		}
	    }
	    if (bad)
	    {
		print_WSStenData(sten);
	        screen("ERROR in states_near_location(), bad state detected\n");
		clean_up(ERROR);
	    }
	}

}		/*end states_near_location*/

LOCAL	void	move_point_along_wall(
	Tparams			*tp,
	TRI	   		**tris,
	POINT			*oldp, 
	POINT			*newp,
	float			*V,
	Front			*fr, 
	float			dt)
{
RECT_GRID  *gr = &topological_grid(fr->interf);
float	   *h = gr->h;
int	   dim = fr->interf->dim;
int	   i, ind;
float	   magv, posn[3], minh, tol = 0.01;
const float    *nor;
TN	   tn;
bool	   found;

	magv = dt*Dot3d(V, tp[0].tan);
	minh = min3(h[0], h[1], h[2]);
	
	//if moving so small, should return. Otherwise the following arg. is unstable.
	if(fabs(magv) < tol*minh)
	{
	    ft_assign(Coords(newp), Coords(oldp), 3*FLOAT);
	    return;
	}

	found = YES;
	if(magv > 0.0)
	{
	    if(tp[1].tnr.tri != NULL)
	        ind = 1;
	    else  if(tp[0].tnr.tri != NULL)
	        ind = 0;
	    else
	    {
	        printf("WARNING  move_point_along_wall, tnr is NULL\n");
	        found = NO;
	    }
	}
	else
	{
	    if(tp[0].tnl.tri != NULL)
	        ind = 0;
	    else  if(tp[1].tnl.tri != NULL)
	        ind = 1;
	    else
	    {
	        printf("WARNING  move_point_along_wall, tnl is NULL\n");
	        found = NO;
	    }
	}

	//print_general_vector("#V=", V, 3, "\n");
	//printf("\n#pos along wall  %24.16e  %24.16e\n", magv, minh);

	if(found)
	{
	    //move along the tan direction for length magv.
	    find_position_along_wall(posn, &tn, magv, &tp[ind], fr);
	    nor = Tri_normal(tn.tri);
	}
	else
	{
	    //try to use nearest_interface_point to get the moved point
	    //WARNING it is dangerous, because the posn may jump across a gap when
	    //using nearest_interface_point.

	    HYPER_SURF		*hs, *hs_on;
	    HYPER_SURF_ELEMENT  *hse;
	    float		t[3], newpos[3];

	    printf("#move point with nearest_pt %d %d %d %d magv=%24.16e\n", 
	        tp[0].tnl.tri, tp[1].tnl.tri,
		tp[0].tnr.tri, tp[1].tnr.tri, magv);
	
	    for (i = 0; i < dim; ++i)
	        newpos[i] = Coords(oldp)[i] + magv*tp[0].tan[i];
	   
	    hs = magv > 0 ? tp[1].hs : tp[0].hs;
	    if(nearest_interface_point(newpos,positive_component(hs),
	    		      fr->interf, NO_SUBDOMAIN, hs,
			      posn, t, &hse, &hs_on) != YES)
	    {
	        printf("ERROR move_point_along_wall, nearest_interface_point fails %24.16e\n", magv);
	        clean_up(ERROR);
	    }
	    
	    print_general_vector("newpos= ", newpos, 3, "\n");
	    print_general_vector("posn= ", posn, 3, "\n");
	   
	    nor = Tri_normal(Tri_of_hse(hse));
	}

	magv = Mag3d(nor);
	//move to the back of the wall, so insert_grid_based_curve_crossings3d
	//can find the intersection points
	for (i = 0; i < dim; ++i)
	    Coords(newp)[i] = posn[i] - tol*h[i]*nor[i]/magv;
}

LOCAL	void	set_up_tangent_params_on_wall(
	Tparams			*tp,
	TRI	   		**tris,
	float			*nor,
	Front			*fr,
	POINT			*p,
	BOND			*b,
	CURVE			*c,
	HYPER_SURF		*hs)
{
SURFACE    *s, *s1, *surfs[2];
BOND_TRI   **btris;

	surfs[0] = NULL;
	surfs[1] = NULL;

	s = Surface_of_hs(hs);
	for(btris = Btris(b); btris && *btris; btris++)
	{
	    s1 = (*btris)->surface;
	    if(s1 == s)
	        continue;
	    //Here ASSUME the positive_comp of the wall is in fluid.
	    if(positive_component(s) == positive_component(s1))
	    {
	        tris[1] = (*btris)->tri;
		surfs[1] = s1;
	    }
	    if(negative_component(s) == positive_component(s1))
	    {
	        tris[0] = (*btris)->tri;
		surfs[0] = s1;
	    }
	}
	if(tris[0] == NULL || tris[1] == NULL)
	{
	    printf("ERROR set_up_tangent_params_on_wall, one side of tris is NULL\n ");
	    printf("comp: %d %d | tris: %d %d\n", 
		    negative_component(s), positive_component(s), 
		    surfs[0], surfs[1]);
	    clean_up(0);
	}

	//onwall tp[0], surfs[0], tris[0]  neg side of the fluid surface
	//onwall tp[1], surfs[1], tris[1]  pos side of the fluid surface
	set_up_wall_tangent_params(tp,p,b,c,nor,surfs,tris, fr->interf,fr);

	if(debugging("set_up_tangent_params_on_wall"))
	{
	    printf("#fluid s comp (%d  %d)\n", negative_component(s), 
	    			positive_component(s));
	    print_Tparams("tp0", &tp[0]);
	    printf("\n");
	    print_Tparams("tp1", &tp[1]);
	    printf("\n");
	    //clean_up(0);
	}
}	/*end SetWSStenData*/

//#bjet2
LOCAL	void	SetWSStenData_from_WallTan(
	WSSten			*sten,
	Tparams			*tp,
	TRI	   		**tris,
	Front			*fr,
	Wave			*wave,
	POINT			*oldp,
	BOND			*oldb,
	CURVE			*oldc,
	HYPER_SURF_ELEMENT 	*oldhse,
	HYPER_SURF		*oldhs,
	float			dt)
{
static  Tan_stencil *tan_sten = NULL;
Locstate   sl, sr;
int 	   i, nsts = sten->nsts - 1;  
//note: the nsts is the # of point FAR FROM the surface in one side.

	if (tan_sten == NULL)
	    tan_sten = alloc_tan_stencil(fr, nsts);

	ClearWSStenData(sten);
	
	if (fr->interf->modified)
	{
	    if (!make_interface_topology_lists(fr->interf))
	    {
		screen("ERROR in SetWSStenData(), "
		       "make_interface_topology_lists(), failed\n");
		clean_up(ERROR);
	    }
	}

	sten->w_type = wave_type(oldhs);
	sten->p = oldp;
	sten->hs = oldhs;
	sten->hse = oldhse;
	sten->front = fr;
	sten->wave = wave;
	sten->dt = dt;

	//tris[0]: wall tris from neg side of the fluid surface
	//tp[0]; tangent setup from tris[0]
	//tris[1]: wall tris from neg side of the fluid surface
	//tp[1]; tangent setup from tris[1]
	fill_tan_stencil_along_wall(tan_sten, fr, &tp[0], tris, -1);
	fill_tan_stencil_along_wall(tan_sten, fr, &tp[1], tris, 1);

	if(debugging("SetWSStenData_from_WallTan"))
	{
	    printf("#Tan_stencil\n");
	    test_print_Tan_stencil(tan_sten);
	    printf("#len %15.8e  %15.8e  %15.8e\n", 
	        Mag3d(sten->nor), Mag3d(tp->tan), Dot3d(tp->tan, sten->nor));
	}

	// Calculate pressure jump due to surface tension 
	sten->pjump = set_pjump_at_wave(oldp,oldhse,oldhs,fr,sten->nor);

	//prop along the tangent of the wall
	ft_assign(sten->nor, tp->tan, 3*sizeof(float));

	slsr(oldp,oldhse,oldhs,&sl,&sr);
	//fill sten from tan_sten
	states_near_location_from_WallTan(sten,Coords(oldp),sten->nor,
				  positive_component(oldhs),
				  negative_component(oldhs),sl,sr,
				  tan_sten, tp[0].ds);

	if(debugging("SetWSStenData_from_WallTan") && the_point(oldp) && NO)
	{
	    printf("#pt in tri\n");
	    print_tri(Tri_of_hse(oldhse), oldhs->interface);
	    printf("#sten from tan_sten\n");
	    print_WSStenData(sten);
	}

}	/*end SetWSStenData*/

EXPORT	void g_point_propagate_along_wall(
	Front		   *fr,
	POINTER		   p2wave,
	POINT		   *oldp,
	BOND		   *oldb,
	CURVE		   *oldc,
	HYPER_SURF_ELEMENT *oldhse,
	HYPER_SURF	   *oldhs,
	POINT		   *newp,
	float		   dt,
	float		   *V)
{
INTERFACE	*intfc = fr->interf;
Wave		*wave = (Wave*)p2wave;
int		i, dim = intfc->dim;
TRI	   	*tris[2];
static Tparams  tp[2];
static WSSten	*sten = NULL;

	//printf("#g_point_propagate_along_wall is called.\n");

	if (sten == NULL) 
	    sten = AllocDefaultWSSten(fr);
	else
	    ClearWSStenData(sten);

        sten->V = NULL;
	normal(oldp,oldhse,oldhs,sten->nor,fr);
	set_up_tangent_params_on_wall(tp,tris,sten->nor,fr,
	    oldp,oldb,oldc,oldhs);

	SetWSStenData_from_WallTan(sten,tp,tris,fr,wave,oldp,oldb,oldc,
	    oldhse,oldhs,dt);

	// Advance front point and states in wall tangent sweep 

	npt_w_speed(sten,left_state(newp),right_state(newp),V);
	
	//print_general_vector("#V=", V, 3, "\n");

	//WARNING need to be updated, so that the point is on the wall.
	//for (i = 0; i < dim; ++i)
	//    Coords(newp)[i] = Coords(oldp)[i] + V[i]*dt;

	move_point_along_wall(tp, tris, oldp, newp, V, fr, dt);
	
	/* Record speed for time step calculation */

	for (i = 0; i < dim; ++i)
	    set_max_front_speed(i,fabs(V[i]),return_obst_state(),
				Coords(newp),fr);

	set_max_front_speed(dim,fabs(scalar_product(V,sten->nor,dim))/sten->dn,
			    return_obst_state(),Coords(newp),fr);

}		/*end g_point_propagate*/

EXPORT	void pseudo_g_point_propagate_along_wall(
	Front		   *fr,
	POINTER		   p2wave,
	POINT		   *oldp,
	BOND		   *oldb,
	CURVE		   *oldc,
	HYPER_SURF_ELEMENT *oldhse,
	HYPER_SURF	   *oldhs,
	POINT		   *newp,
	float		   dt,
	float		   *V)
{
INTERFACE	*intfc = fr->interf;
Wave		*wave = (Wave*)p2wave;
int		i, dim = intfc->dim;
TRI	   	*tris[2];
static Tparams  tp[2];
static WSSten	*sten = NULL;

	if (sten == NULL) 
	    sten = AllocDefaultWSSten(fr);
	else
	    ClearWSStenData(sten);

        sten->V = NULL;
	normal(oldp,oldhse,oldhs,sten->nor,fr);
	set_up_tangent_params_on_wall(tp,tris,sten->nor,fr,
	    oldp,oldb,oldc,oldhs);

	SetWSStenData_from_WallTan(sten,tp,tris,fr,wave,oldp,oldb,oldc,
	    oldhse,oldhs,dt);
	
	// Advance front point and states in wall tangent sweep 
	//npt_w_speed(sten,left_state(newp),right_state(newp),V);
	
	fourth_order_point_propagate(fr,wave,oldp,newp,oldhse,oldhs,dt,V);
	ft_assign(left_state(newp),left_state(oldp),fr->sizest);
	ft_assign(right_state(newp),right_state(oldp),fr->sizest);
	
	for (i = 0; i < dim; ++i)
	    V[i] = (Coords(newp)[i] - Coords(oldp)[i])/dt;
	
	//WARNING need to be updated, so that the point is on the wall.
	move_point_along_wall(tp, tris, oldp, newp, V, fr, dt);
}		/*end g_point_propagate*/


EXPORT	void	set_point_propagate(
	Front		*fr,
	bool	unsplit_pt_prop)
{
	OperatorSplitPointPropagate(fr) = g_point_propagate;
	WSStenAllocator(fr) = g_AllocWSSten;
	nsts_WSSten(fr) = 3;
	stype_WSSten(fr) = TGAS_STATE;
	if (debugging("oldws"))
	{
	    /**FOR OLD CODE FUNCTIONALITY TURN ON THIS DEBUG STRING**/
	    nsts_WSSten(fr) = 2;
	    stype_WSSten(fr) = GAS_STATE;
	}
	if (fr->rect_grid->dim == 2)
	{
	    UnsplitPointPropagate(fr) = g_unsplit_point_propagate;
	    SetUSWSStencil(fr) = g_set_USWSSten2d;
	}
	if ((unsplit_pt_prop == YES) && (UnsplitPointPropagate(fr) != NULL))
	{
	    fr->surface_propagate = NULL;
	    fr->_point_propagate = UnsplitPointPropagate(fr);
	    fr->snd_node_propagate = NULL;
	    fr->tan_curve_propagate = NULL;
	}
	else
	{
	    fr->surface_propagate = NULL;
	    fr->_point_propagate = OperatorSplitPointPropagate(fr);
	    //#bjet2
	    fr->_point_propagate_along_wall = g_point_propagate_along_wall;
	}
}		/*end set_point_propagate*/

LOCAL	void g_unsplit_point_propagate(
	Front		*fr,
	POINTER		p2wave,
	POINT		*oldp,
	POINT		*newp,
	HYPER_SURF_ELEMENT *oldhse,
	HYPER_SURF	*oldhs,
	float		dt,
	float		*V)
{
	INTERFACE          *intfc = fr->interf;
	Wave		   *wave = (Wave*)p2wave;
	float	           **normals;   /* unit normals at front points */
	int	           i, j, dim = intfc->dim;
	int	           w_type;
	float	           pjump;
#if defined(DEBUG_POINT_PROPAGATE)
	static char        fname[] = "g_unsplit_point_propagate";
#endif /* defined(DEBUG_POINT_PROPAGATE) */

	static Locstate       sl = NULL, sr = NULL;
	static int	      nor_rad, tan_rad;
	static	USWSSten2d    *uswssten = NULL;
	static float	      alpha, vzbar[MAXD];

	entered_point_propagate(fname,oldhs);

	w_type = (is_subdomain_boundary(oldhs)) ? SUBDOMAIN_BOUNDARY :
			wave_type(oldhs);
	if (passive_point_propagate(w_type,oldp,newp,V,intfc) == YES)
	{
	    left_point_propagate(fname,newp,V,dt,dim,NULL);
	    return;
	}

		/* Allocate storage */

	if (uswssten == NULL) 
	{
	    nor_rad = nsts_WSSten(fr);
	    tan_rad = fr->npts_tan_sten/2;
	    uswssten = AllocUSWSSten2d(fr,nor_rad,tan_rad);
	    alloc_state(fr->interf,&sl,fr->sizest);
	    alloc_state(fr->interf,&sr,fr->sizest);
	    if (is_rotational_symmetry())
		alpha = rotational_symmetry();
	}
	normals = nor_vec(uswssten);

#if defined(DEBUG_POINT_PROPAGATE)
	if (debug_point_propagate == YES) 
	    print_point_propagate_data(oldp,oldhse,oldhs,dim);
#endif /* defined(DEBUG_POINT_PROPAGATE) */

	set_USWSSten2d(uswssten,oldp,oldhse,oldhs,fr,wave,dt);

	/* Calculate pressure jump due to surface tension */
	for (i = -tan_rad; i <= tan_rad; ++i)
	{
	    pjump = set_pjump_at_wave(oldp,oldhse,oldhs,fr,normals[i]);
	    nor_fr(uswssten)[i]->pjump = pjump;
	}


#if defined(DEBUG_POINT_PROPAGATE)
	if (debug_point_propagate == YES) 
	{
	    static const char *fmt[3] = {"%-22s%s","%-22s%s","%-33s%s"};
	    static const char *vname[3] = {"X", "Y", "Z"};

	    g_PrintUSWSSten2d(uswssten);
	    for (j = -tan_rad; j <= tan_rad ; ++j)
	    {
		static char label[21];

	        (void) sprintf(label,"\nnor[%d] = ",j);
	        print_general_vector(label,normals[j],dim,"\n");
	    }
	    print_general_vector(", tngt = ",uswssten->tngt,dim,"\n");
	    (void) printf("%-9s","INDEX");
	    (void) printf(fmt[dim-1],"ON FRONT POINT","    ");
	    (void) printf(fmt[dim-1],"LEFT OFF FRONT POINT","    ");
	    (void) printf(fmt[dim-1],"RIGHT OFF FRONT POINT","\n");
	    (void) printf("%-9s","I");
	    for (j = 0; j < 3; ++j)
	    {
	    	for (i = 0; i < dim; ++i)
	    	    (void) printf("%-11s",vname[i]);
	    	(void) printf("%s",(j < 2) ? "    " : "\n");
	    }
	    for (j = -tan_rad; j <= tan_rad ; ++j)
	    {
	        (void) printf("%-9d",j);
	        for (i = 0; i < dim; ++i)
	    	    (void) printf("%-11g",Coords(tan_fr(uswssten)[0]->p[j])[i]);
	        (void) printf("    ");
	        for (i = 0; i < dim; ++i)
	    	    (void) printf("%-11g",Coords(tan_fr(uswssten)[-1]->p[j])[i]);
	        (void) printf("    ");
	        for (i = 0; i < dim; ++i)
	    	    (void) printf("%-11g",Coords(tan_fr(uswssten)[1]->p[j])[i]);
	        (void) printf("\n");
	    }
	    (void) printf("\n");
	    for (j = -tan_rad; j <= tan_rad ; ++j)
	    {
	         (void) printf("right_state(tan_fr(uswssten)[-1]->p[%d])\n",j);
	         (*fr->print_state)(right_state(tan_fr(uswssten)[-1]->p[j]));
	         (void) printf("left_state(tan_fr(uswssten)[0]->p[%d])\n",j);
	         (*fr->print_state)(left_state(tan_fr(uswssten)[0]->p[j]));
	         (void) printf("right_state(tan_fr(uswssten)[0]->p[%d])\n",j);
	         (*fr->print_state)(right_state(tan_fr(uswssten)[0]->p[j]));
	         (void) printf("left_state(tan_fr(uswssten)[1]->p[%d])\n",j);
	         (*fr->print_state)(left_state(tan_fr(uswssten)[1]->p[j]));
	    }
	}
#endif /* defined(DEBUG_POINT_PROPAGATE) */


	/* Translate to a frame where wave has approximately zero vertical
	*  velocity.
	*  NOTE:  Should this (or something similar) be applied in the
	*  	  non-symmetric runs also?
	*/

	if (is_rotational_symmetry() && alpha > 0.0 && dim == 2)
	{
	    /* Get a first guess at wave speed */
	    if (is_obstacle_state(tan_fr(uswssten)[0]->rightst[0]))
	    	vzbar[1] = -vel(1,tan_fr(uswssten)[0]->leftst[0]);
	    else if (is_obstacle_state(tan_fr(uswssten)[0]->leftst[0]))
	    	vzbar[1] = -vel(1,tan_fr(uswssten)[0]->rightst[0]);
	    else
	    	vzbar[1] = -0.5*(vel(1,tan_fr(uswssten)[0]->rightst[0]) +
	    			 vel(1,tan_fr(uswssten)[0]->leftst[0]));

	    add_velocity_to_state(tan_fr(uswssten)[0]->leftst[0],vzbar);
	    add_velocity_to_state(tan_fr(uswssten)[0]->rightst[0],vzbar);

	    V[1] = 0.0;

	    /* Improve this guess by a call to w_speed */
	    if (w_type >= FIRST_PHYSICS_WAVE_TYPE)
	    	w_speed(Coords(tan_fr(uswssten)[0]->p[0]),
	    		tan_fr(uswssten)[0]->leftst[0],
	    		tan_fr(uswssten)[0]->rightst[0],
	    		sl,sr,V,0.0,normals[0],w_type,fr);

	    V[0] = 0.0;
	    V[1] = -V[1];

	    /* Remove remaining vertical wave speed */
	    add_velocity_to_state(tan_fr(uswssten)[0]->leftst[0],V);
	    add_velocity_to_state(tan_fr(uswssten)[0]->rightst[0],V);
			
	    vzbar[1] += V[1];

	    for (j = -tan_rad; j<=tan_rad; ++j)
	    {
	        add_velocity_to_state(tan_fr(uswssten)[-1]->rightst[j],vzbar);
	        add_velocity_to_state(tan_fr(uswssten)[1]->leftst[j],vzbar);

	        if (j == 0) continue;

	        add_velocity_to_state(tan_fr(uswssten)[0]->leftst[j],vzbar);
	        add_velocity_to_state(tan_fr(uswssten)[0]->rightst[j],vzbar);
	    }
	}


	unsplit_w_speed2d(uswssten,left_state(newp),right_state(newp),V);

	/*
	 * !!!TODO!!!TODO!!!TODO!!!TODO!!!TODO!!!
	 *
	 * Make this a point propagation sub-option when unsplit operator
	 * has been selected.
	 */

	if (is_scalar_wave(w_type) && debugging("experiment"))
	{
	    Locstate sl =  left_state(newp);
	    Locstate sr = right_state(newp);
	    float rhol = Dens(sl), rhor = Dens(sr);
	    float A = fabs(rhol - rhor)/(rhol + rhor);

	    for (i = 0; i < dim; ++i)
	        V[i] = (1.0-A)*V[i]+A*(mom(i,sl)+mom(i,sr))/(rhol+rhor);
	}
	for (i = 0; i < dim; ++i)
	    Coords(newp)[i] = Coords(oldp)[i] + V[i]*dt;
	for (i = 0; i < dim; ++i)
	{
	    set_max_front_speed(i,fabs(V[i]),return_obst_state(),
				Coords(newp),fr);
	}
	set_max_front_speed(dim,
			    fabs(scalar_product(V,normals[0],dim))/
				uswssten->dn,
			    return_obst_state(),Coords(newp),fr);

	if (is_rotational_symmetry() && alpha > 0.0 && dim == 2)
	{
	    /* Add back vertical velocity removed earlier */

	    vzbar[1] = -vzbar[1];

	    add_velocity_to_state(left_state(newp),vzbar);
	    add_velocity_to_state(right_state(newp),vzbar);
	    if (w_type >= FIRST_PHYSICS_WAVE_TYPE)
	    	Coords(newp)[1] += vzbar[1]*dt;

	    V[1] += vzbar[1];
	    set_max_front_speed(1,fabs(V[1]),return_obst_state(),
	    			Coords(newp),fr);
	    set_max_front_speed(dim,
		                fabs(scalar_product(V,normals[0],dim))/
				    uswssten->dn,
		                return_obst_state(),Coords(newp),fr);
	}

	left_point_propagate(fname,newp,V,dt,dim,NULL);
}		/*end g_unsplit_point_propagate*/


/*
*			g_set_USWSSten2d():
*
*	Sets up the unsplit wave speed stencil to be used by unsplit_w_speed2d()
*	for the computation of the unsplit point propagate.
*
*	The states and coords in a USWSSten2D are represented twice:
*	1) as a set of WSSten(s) normal to the curve
*	2) as a set of Tan_stencil(s) parallel to the curve.
*	Current version is only valid for TWO space dimensions.
*       
*	Orientation with respect to curve is as indicated.
*	uswssten->tngt = unit tangent in the POSITIVE direction for the curve.
*	uswssten->nor  = unit normal 90 deg. CLOCKWISE from tangent.
*	Nor always points from LEFT to RIGHT or from NEG comp to POS comp
*	and (nor,tngt) are RIGHT-HANDED local coordinates.
*	The stencil returned is indicated in the diagram(s) below for
*	nor_rad == tan_rad == 2.
*
*       Indexing of the WSSten(s): stencil index then point index.
*
*	PL[-2][2]    PL[-1][2]     PL[0][2]       PL[1][2]   PL[2][2]  
*	  |           |             |              |          |       
*	  |           |             |              |          |       
*	  |           |             |              |          |       
*	PL[-2][1]    PL[-1][1]     PL[0][1]       PL[1][1]   PL[2][1]  
*	  |           |             |              |          |       
*	  |           |             |              |          |       
*	  |           |             |              |          |       
*	 P[-2]--->--->P[-1]---->----P[0] --->----->P[1]--->---P[2]---->
*	  |           |             |      ds      |          |       
*	  |           |             |              |          |       
*	  |           |             |              |          |       
*	PR[-2][1]    PR[-1][1]     PR[0][1]       PR[1][1]   PR[2][1]     
*	  |           |             |              |          |       
*	 \|/         \|/ dn        \|/            \|/        \|/      
*	  |           |             |              |          |       
*	PR[-2][2]    PR[-1][2]     PR[0][2]       PR[1][2]   PR[2][2]     
*

*	Indexing of the Tan_stencil(s): point index then stencil index.
*
*	P[-2][-2]   P[-1][-2]     P[0][-2]       P[1][-2]   P[2][-2] 
*	  |           |             |              |          |       
*	  |           |             |              |          |       
*	  |           |             |              |          |       
*	P[-2][-1]   P[-1][-1]     P[0][-1]       P[1][-1]   P[2][-1] 
*	  |           |             |              |          |       
*	  |           |             |              |          |       
*	  |           |             |              |          |       
*	 P[-2]--->--->P[-1]---->----P[0]--->----->P[1]--->---P[2]---->
*	  |           |             |     ds       |          |       
*	  |           |             |              |          |       
*	  |           |             |              |          |       
*	P[-2][1]    P[-1][1]      P[0][1]        P[1][1]    P[2][1]     
*	  |           |             |              |          |       
*	 \|/         \|/ dn        \|/            \|/        \|/      
*	  |           |             |              |          |       
*	P[-2][2]    P[-1][2]      P[0][2]        P[1][2]    P[2][2]     
*
*/

LOCAL	void	g_set_USWSSten2d(
	USWSSten2d		*uswssten,
	POINT			*p,
	HYPER_SURF_ELEMENT	*hse,
	HYPER_SURF		*hs,
	Front			*fr,
	Wave			*wave,
	float			dt)
{
	CURVE		*c = Curve_of_hs(hs);/*ASSUMES TWOD*/
	BOND		*b = Bond_of_hse(hse);/*ASSUMES TWOD*/
	Tan_stencil	*on_fr = tan_fr(uswssten)[0];
	WSSten		*nsten;
	WSSten		*nsten0 = nor_fr(uswssten)[0];
	WSSten		*tan_ans = uswssten->tan_ans;
	float           **normals = nor_vec(uswssten);

	int		i, j, dim, tan_rad = uswssten->tan_rad;
 
	debug_print("USWS","Entered g_set_USWSSten2d()\n");

	uswssten->fr = fr;
	uswssten->wave = wave;
	uswssten->hs = hs;
	set_uswssten_geometry_and_center(uswssten,p,hse,hs,fr,dt);

	/*
	 * Load data at tan_rad points along curve from stencil center
	 * in the indicated direction at separation uswssten->ds.
	 */
	states_at_distance_along_curve(p,b,c,NEGATIVE_ORIENTATION,uswssten->ds,
				       tan_rad,on_fr->leftst-1,on_fr->rightst-1,
				       on_fr->hs-1,on_fr->hse-1,on_fr->t-1,
				       on_fr->p-1,fr);
	states_at_distance_along_curve(p,b,c,POSITIVE_ORIENTATION,uswssten->ds,
				       tan_rad,on_fr->leftst+1,on_fr->rightst+1,
				       on_fr->hs+1,on_fr->hse+1,on_fr->t+1,
				       on_fr->p+1,fr);

	dim = fr->rect_grid->dim;
	for (i = -tan_rad; i <= tan_rad; ++i)
	{
	    Locstate	sl = on_fr->leftst[i];
	    Locstate	sr = on_fr->rightst[i];

	    nsten = nor_fr(uswssten)[i];
	    nsten->hs = hs;
	    if (i != 0)
	    {
		find_propagation_direction(on_fr->p[i],on_fr->hse[i],
					   on_fr->hs[i],normals[i],fr);
	    }
	    for (j = 0; j < dim; ++j)
	    {
		nsten->coords[j] = Coords(on_fr->p[i])[j];
		nsten->lcrds[0][j] = nsten->rcrds[0][j] = nsten->coords[j];
		nsten->nor[j] = normals[i][j];
	    }
	    nsten->dn = uswssten->dn;
	    nsten->pcomp = positive_component(hs);
	    nsten->ncomp = negative_component(hs);
	    nsten->w_type = wave_type(hs);
	    nsten->dt = dt;
	    nsten->pjump = 0.0;
	    nsten->front = fr;
	    nsten->wave = wave;
	    set_state(nsten->sl[0],state_type(sl),sl);
	    set_state(nsten->sr[0],state_type(sr),sr);
	}
	tan_ans->dn = nsten0->dn;
	tan_ans->dt = nsten0->dt;
	tan_ans->pcomp = nsten0->pcomp;
	tan_ans->ncomp = nsten0->ncomp;
	tan_ans->w_type = nsten0->w_type;
	tan_ans->pjump = nsten0->pjump;
	tan_ans->front = nsten0->front;
	tan_ans->wave = nsten0->wave;
	tan_ans->hs = nsten0->hs;

	set_USWSSten2d_off_front_states(uswssten,hs,fr,wave);
	for (j = 0; j < dim; ++j)
	{
	    tan_ans->coords[j] = nsten0->coords[j];
	    tan_ans->nor[j] = nsten0->nor[j];
	    for (i = 0; i < nsten0->nsts; ++i)
	    {
		tan_ans->lcrds[i][j] = nsten0->lcrds[i][j];
		tan_ans->rcrds[i][j] = nsten0->rcrds[i][j];
	    }
	}
	debug_print("USWS","Left g_set_USWSSten2d()\n");

}	/*end g_set_USWSSten2d*/

EXPORT	void	set_uswssten_geometry_and_center(
	USWSSten2d		*uswssten,
	POINT			*p,
	HYPER_SURF_ELEMENT	*hse,
	HYPER_SURF		*hs,
	Front			*fr,
	float			dt)
{
	BOND		*b = Bond_of_hse(hse);/*ASSUMES TWOD*/
	CURVE		*c = Curve_of_hs(hs);/*ASSUMES TWOD*/
	Tan_stencil	*on_fr = tan_fr(uswssten)[0];
	float           **normals = nor_vec(uswssten);
	float		*tngt = uswssten->tngt;
	float		*h = fr->rect_grid->h;
	int             tan_rad = uswssten->tan_rad;
	int		i, j, dim = fr->rect_grid->dim;

	/*
	 * Calculate normal, tangent directions and
	 * set state at stencil center.
	 */

	uswssten->dt = dt;
	find_propagation_direction(p,hse,hs,normals[0],fr);

	/* By default set single normal, uswssten is FLAT. */
	for (i = -tan_rad; i <= tan_rad; ++i)
	    for (j = 0; j < dim; ++j)
		normals[i][j] = normals[0][j];

	tangent(p,b,c,tngt,fr);
	uswssten->dn = grid_size_in_direction(normals[0],h,dim);
	uswssten->ds = grid_size_in_direction(tngt,h,dim);

		/* find the stencil states */

	for (i = 0; i < dim; ++i)
	    Coords(on_fr->p[0])[i] = Coords(p)[i];
	ft_assign(left_state(on_fr->p[0]),
		left_state_at_point_on_curve(p,b,c),fr->sizest);
	ft_assign(right_state(on_fr->p[0]),
		right_state_at_point_on_curve(p,b,c),fr->sizest);
	if (p == b->start)
	    on_fr->t[0] = 0.0;
	else if (p == b->end)
	    on_fr->t[0] = 1.0;
	else
	    on_fr->t[0] = separation(p,b->start,dim)/bond_length(b);
	on_fr->hse[0] = Hyper_surf_element(b);
	on_fr->hs[0] = Hyper_surf(c);
}	/*end set_stencil_geometry_and_center*/

EXPORT	void	set_USWSSten2d_off_front_states(
	USWSSten2d		*uswssten,
	HYPER_SURF		*hs,
	Front			*fr,
	Wave			*wave)
{
	COMPONENT      pos_comp = positive_component(hs);
	COMPONENT      neg_comp = negative_component(hs);
	INTERFACE      *intfc = fr->interf;
	Tan_stencil    *on_fr = tan_fr(uswssten)[0];
	Tan_stencil    *off_frl, *off_frr;
	size_t	       sizest = fr->sizest;
	float          **normals = nor_vec(uswssten);
	int            i, j, l, dim = fr->rect_grid->dim;
	int            tan_rad = uswssten->tan_rad;
	int	       nor_rad = uswssten->nor_rad;

	debug_print("USWS","Entered set_USWSSten2d_off_front_states()\n");

	/*
	 * Load coords and states into Tan_stencil(s), then copy to WSSten(s).
	 */
	for (l = 1; l < nor_rad; ++l)
	{
	    off_frl = tan_fr(uswssten)[-l];
	    off_frr = tan_fr(uswssten)[l];
	    for (j = -tan_rad; j <= tan_rad; ++j)
	    {
	        off_frl->hse[j] = off_frr->hse[j] = NULL;
	        off_frl->hs[j] = off_frr->hs[j] = NULL;
	        off_frl->t[j] = off_frr->t[j] = ERROR_FLOAT;
	        for (i = 0; i < dim; ++i)
	        {
	            Coords(off_frl->p[j])[i] =
	                Coords(on_fr->p[j])[i] - l*normals[j][i]*uswssten->dn;
	            Coords(off_frr->p[j])[i] =
	                Coords(on_fr->p[j])[i] + l*normals[j][i]*uswssten->dn;
	        }
	        if (is_obstacle_state(left_state(on_fr->p[j])))
		{
	            obstacle_state(intfc,right_state(off_frl->p[j]),sizest);
		}
	        else
	        {
	            hyp_solution(Coords(off_frl->p[j]),neg_comp,hs,
				 NEGATIVE_SIDE,fr,wave,
				 right_state(off_frl->p[j]),
				 left_state(on_fr->p[j]));
	        }
	        if (is_obstacle_state(right_state(on_fr->p[j])))
		{
	            obstacle_state(intfc,left_state(off_frr->p[j]),sizest);
		}
	        else
	        {
	            hyp_solution(Coords(off_frr->p[j]),pos_comp,hs,
				 POSITIVE_SIDE,fr,wave,
				 left_state(off_frr->p[j]),
				 right_state(on_fr->p[j]));
	        }
	    }
	}
	for (j = -tan_rad; j <= tan_rad; ++j)
	{
	    WSSten *nsten = nor_fr(uswssten)[j];
	    for (l = 1; l < nor_rad; ++l)
	    {
		POINT	    *pl = tan_fr(uswssten)[-l]->p[j];
		Locstate    sl = right_state(pl);
		POINT	    *pr = tan_fr(uswssten)[l]->p[j];
		Locstate    sr = left_state(pr);
	        for (i = 0; i < dim; ++i)
		{
		    nsten->lcrds[l][i] = Coords(pl)[i];
		    nsten->rcrds[l][i] = Coords(pr)[i];
		    set_state(nsten->sl[l],state_type(sl),sl);
		    set_state(nsten->sr[l],state_type(sr),sr);
		}
	    }
	}
	debug_print("USWS","Left set_USWSSten2d_off_front_states()\n");
}	/*end set_USWSSten2d_off_front_states*/


EXPORT	void	g_unsplit_w_speed2d(
	USWSSten2d   	*uswssten,
	Locstate	ansl,
	Locstate	ansr,
	float           *V)
{
	HYPER_SURF      *oldhs = uswssten->hs;
	Front           *fr = uswssten->fr;
	COMPONENT       neg_comp = negative_component(oldhs);
	COMPONENT       pos_comp = positive_component(oldhs);
	CURVE		*oldc = Curve_of_hs(tan_fr(uswssten)[0]->hs[0]);
	Tan_stencil	*on_fr, *off_frl, *off_frr;
	float		ds = uswssten->ds;
	float		dt = uswssten->dt;
	float           *tngt = uswssten->tngt;
	int		j, tan_rad = uswssten->tan_rad;
	int		nor_rad = uswssten->nor_rad;
	int		i, dim = fr->rect_grid->dim;
	static float	**W = NULL;	/* Wave speeds */

	debug_print("USWS","Entered g_unsplit_w_speed2d()\n");

	if (W == NULL)
	{
	    bi_array(&W,tan_fr(uswssten)[0]->npts,MAXD,FLOAT);
	    W += tan_rad;
	}

	on_fr = tan_fr(uswssten)[0];

	if ((fr->step % 2) == 0)
	{
	    Tan_stencil	*nor_ans = uswssten->nor_ans;
	    for (j = -tan_rad; j <= tan_rad ; ++j)
	    {
		WSSten	*sten = nor_fr(uswssten)[j];
	    	nor_ans->hse[j] = on_fr->hse[j];
	    	nor_ans->hs[j] = on_fr->hs[j];
	    	nor_ans->t[j] = on_fr->t[j];
	    	npt_w_speed(sten,left_state(nor_ans->p[j]),
			    right_state(nor_ans->p[j]),W[j]);
	    	for (i = 0; i < dim; ++i)
	    	    Coords(nor_ans->p[j])[i] =
	    	    	Coords(on_fr->p[j])[i] + W[j][i]*dt;
	    }
	    if (debugging("USWS"))
	    {
		(void) printf("uswssten->nor_ans after normal sweeps\n");
		print_Tan_stencil(uswssten->fr,uswssten->nor_ans);
	    }
	    for (i = 0; i < dim; ++i)
	    	V[i] = W[0][i];
	    nor_ans->newhs = Hyper_surf(oldc);
	    nor_ans->dir = tngt;
	    npt_tang_solver(ds,dt,nor_ans,ansl,ansr,fr);
	}
	else
	{
	    WSSten    *tan_ans = uswssten->tan_ans;

	    on_fr->newhs = Hyper_surf(oldc);
	    on_fr->dir = tngt;
	    npt_tang_solver(ds,dt,on_fr,tan_ans->sl[0],tan_ans->sr[0],fr);
	    for (i = 1; i < nor_rad; ++i)
	    {
	        off_frl = tan_fr(uswssten)[-i];
	        off_frr = tan_fr(uswssten)[i];
		off_frl->newhs = NULL;
		off_frl->comp = neg_comp;
		off_frl->states = off_frl->rightst;
		off_frl->dir = tngt;
	        one_side_npt_tang_solver(ds,dt,off_frl,tan_ans->sl[i],fr);
		off_frr->newhs = NULL;
		off_frr->comp = pos_comp;
		off_frr->states = off_frr->leftst;
		off_frr->dir = tngt;
	        one_side_npt_tang_solver(ds,dt,off_frr,tan_ans->sr[i],fr);
	    }
	    if (debugging("USWS"))
	    {
		(void) printf("uswssten->tan_ans after tangential sweeps\n");
		PrintWSSten(uswssten->tan_ans);
	    }
	    npt_w_speed(tan_ans,ansl,ansr,V);
	}
	debug_print("USWS","Left g_unsplit_w_speed2d()\n");
}		/*end g_unsplit_w_speed2d*/


/*
*		g_pt_prop_by_w_speed():
*
*	Basically the same as g_point_propagate(),  but uses
*	w_speed() instead of npt_w_speed().  This function
*	is primarily intended for use in estimating node velocities
*	for nodes produced by bifurcations where the hyp solution
*	has not been generated yet.
*/



EXPORT void g_pt_prop_by_w_speed(
	Front		*fr,
	POINT		*oldp,
	POINT		*newp,
	HYPER_SURF_ELEMENT *oldhse,
	HYPER_SURF	*oldhs,
	float		dt,
	float		*V)
{
	int		i, dim = fr->interf->dim;
	int		w_type;
	float		nor[MAXD];  /* components of unit normal to front */
	float		dn;
	float		pjump;
	Locstate	sl,sr;	    /* states at left, right sides of front */
#if defined(DEBUG_POINT_PROPAGATE)
	static char	fname[] = "g_pt_prop_by_w_speed";
#endif /* defined(DEBUG_POINT_PROPAGATE) */

	entered_point_propagate(fname,oldhs);

	w_type = (is_subdomain_boundary(oldhs)) ? SUBDOMAIN_BOUNDARY :
			wave_type(oldhs);
	if (passive_point_propagate(w_type,oldp,newp,V,fr->interf) == YES)
	{
	    left_point_propagate(fname,newp,V,dt,dim,NULL);
	    return;
	}

#if defined(DEBUG_POINT_PROPAGATE)
	if (debug_point_propagate == YES) 
		print_point_propagate_data(oldp,oldhse,oldhs,dim);
#endif /* defined(DEBUG_POINT_PROPAGATE) */

	/* Calculate normal direction */
	/* and positions of normally displaced states */

	find_propagation_direction(oldp,oldhse,oldhs,nor,fr);

	/* Find states near front for normal sweep */

	slsr(oldp,oldhse,oldhs,&sl,&sr);

	/* Calculate pressure jump due to surface tension */
	 
	pjump = set_pjump_at_wave(oldp,oldhse,oldhs,fr,nor);


#if defined(DEBUG_POINT_PROPAGATE)
	if (debug_point_propagate == YES) 
	{
	    for (i = 0; i < dim; ++i)
	    {
	    	(void) printf("nor[%d] = %g",i,nor[i]);
	    	(void) printf("%s",(i<dim-1) ? ", " : "\n");
	    }
	    verbose_print_state("sl",sl);
	    verbose_print_state("sr",sr);
	}
#endif /* defined(DEBUG_POINT_PROPAGATE) */

	/* Advance front point and states in normal sweep */

	w_speed(Coords(oldp),sl,sr,left_state(newp),right_state(newp),
		V,pjump,nor,w_type,fr);

	for (i = 0; i < dim; ++i)
	    Coords(newp)[i] = Coords(oldp)[i] + V[i]*dt;

	/* Record speed for time step calculation */

	dn = grid_size_in_direction(nor,fr->rect_grid->h,dim);
	for (i = 0; i < dim; ++i)
	{
	    set_max_front_speed(i,fabs(V[i]),return_obst_state(),
	    			Coords(newp),fr);
	}
	set_max_front_speed(dim,fabs(scalar_product(V,nor,dim))/dn,
			    return_obst_state(),Coords(newp),fr);

	left_point_propagate(fname,newp,V,dt,dim,NULL);
}		/*end g_pt_prop_by_w_speed*/

/*
* Note on dynamic surface tension:
*  Dynamic surface tension is meant to stabilize contacts against shear
*  instabilities.  According to the incompressible linear theory,
*  if the wave number k of an instability is greater than a certain
*  critical value (determined in part by the surface tension) then this
*  instability will be suppressed.  Given that we want to suppress
*  instabilities on the order of M grid zones (say 2 or 3 grid zones)
*  the surface tension T at a point must satisfy
*
*	T >= M * rho1 * rho2 * U^2 * dx/(2Pi * (rho1 + rho2))
*
*  where rho1 and rho2 are the densities on each side of the contact,
*  U is the shear velocity (u1 - u2) and dh is the grid spacing.
*  See Chandrasekhar, Hydrodynamic and Hydromagnetic Stability, 1961
*
*/

EXPORT	float set_pjump_at_wave(
	POINT		   *p,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF	   *hs,
	Front		   *fr,
	float		   *n)
{
	Locstate	sl, sr;
	float		pjump = 0.0;
	int		dim = fr->rect_grid->dim;
	
	slsr(p,hse,hs,&sl,&sr);

	if (    (is_scalar_wave(wave_type(hs))) &&
		( (surface_tension(hs) != 0.0) ||
			(Params(sl) != NULL &&
			Params(sl)->avisc.dynamic_st != 0.0 ))     )
	{
	    float M, tension;

	    tension = surface_tension(hs);

	    M = Params(sl)->avisc.dynamic_st;
	    if (n && M > 0.0)
	    {
	    	float	rho1, rho2, U[MAXD];
	    	float	vl[MAXD], vr[MAXD];
	    	float	vln, vrn;
	    	float dh, *h = fr->rect_grid->h;
	    	int	i;

	    	for (i = 0; i < dim; ++i)
	    	{
		    vl[i] = vel(i,sl);
		    vr[i] = vel(i,sr);
		}
		vln = scalar_product(vl,n,dim);
		vrn = scalar_product(vl,n,dim);
		for (i = 0; i < dim; ++i)
		{
		    vl[i] -= vln*n[i];
		    vr[i] -= vrn*n[i];
		    U[i] = vl[i] - vr[i];
		}

		rho1 = Dens(sl);
		rho2 = Dens(sr);

		dh = grid_size_in_direction(n,h,dim);

		tension += M*rho1*rho2*scalar_product(U,U,dim)*dh / 
						(2*PI*(rho1+rho2));
	    }
	    switch (dim)
	    {
	    case 2: 
	    	pjump = tension*mean_curvature_at_point(p,hse,hs,fr);
		break;
	    case 3:
	    	pjump = 2.0*tension*(p->curvature);
		break;
	    }
	    
	    //printf("#pjump %15.8e %15.8e %15.8e\n", tension, p->curvature, pjump);
	    
	    if(debugging("surface_tension"))
	    {
	    	int i, dim = fr->rect_grid->dim;

	    	(void) printf("position =");
	    	for (i = 0; i < dim; ++i)
	    	    (void) printf(" %g",Coords(p)[i]);
		(void) printf("\n");
		(void) printf("curvature = %g\n",
			      mean_curvature_at_point(p,hse,hs,fr));
		(void) printf("surface tension = %g\n",tension);
		(void) printf("pressure jump = %g\n",pjump);
	    }
	}
#if defined(COMBUSTION_CODE)
	else if (is_shock_wave(wave_type(hs)) &&
		 Composition_type(sl) != PURE_NON_REACTIVE)
	{
	    //pjump = mean_curvature_at_point(p,hse,hs,fr);
	    pjump = p->curvature;
	}
#endif /* defined(COMBUSTION_CODE) */
	
	return pjump;
}		/*end set_pjump_at_wave*/


EXPORT	bool passive_point_propagate(
	int		w_type,
	POINT		*oldp,
	POINT		*newp,
	float		*V,
	INTERFACE	*intfc)
{
	int		i;
	int		dim = intfc->dim;

	switch (w_type)
	{
	case PASSIVE_BOUNDARY:
	case SUBDOMAIN_BOUNDARY:
	    obstacle_state(intfc,left_state(newp),size_of_state(intfc));
	    obstacle_state(intfc,right_state(newp),size_of_state(intfc));
	    for (i = 0; i < dim; ++i)
	    {
	    	V[i] = 0.0;
	    	Coords(newp)[i] = Coords(oldp)[i];
	    }
	    return YES;
	default:
	    break;
	}
	return NO;
}		/*end passive_point_propagate*/

EXPORT	USWSSten2d	*AllocUSWSSten2d(
	Front	*fr,
	int	nor_rad,
	int	tan_rad)
{
	INTERFACE	*intfc;
	byte		*store;
	USWSSten2d	*uswssten;
	int		i, j;
	int		nts, nns;
	size_t		size, sizest;
	size_t		tan_sten_offset, nor_sten_offset;

	nts = 2*nor_rad - 1;	/* Number of Tan_stencil(s) in 2D stencil */
	nns = 2*tan_rad + 1;	/* Number of WSSten(s) in 2D stencil */
	size = sizeof(USWSSten2d);
	tan_sten_offset = (size%sizeof(Tan_stencil*) != 0) ?
	    (size/sizeof(Tan_stencil*) + 1)*sizeof(Tan_stencil*) : size;
	size = tan_sten_offset + nts*sizeof(Tan_stencil*);
	nor_sten_offset = (size%sizeof(WSSten*) != 0) ?
	    (size/sizeof(WSSten*) + 1)*sizeof(WSSten*) : size;
	size = nor_sten_offset + nns*sizeof(WSSten*);

	scalar(&store,size);
	uswssten = (USWSSten2d*)store;
	uswssten->tan_rad = tan_rad;
	uswssten->nor_rad = nor_rad;
	tan_fr(uswssten) = (Tan_stencil**)(store + tan_sten_offset +
					 (nts/2)*sizeof(Tan_stencil*));
	nor_fr(uswssten) = (WSSten**)(store + nor_sten_offset +
					 (nns/2)*sizeof(WSSten*));
	for (i = -nts/2; i <= nts/2; ++i)
	    tan_fr(uswssten)[i] = alloc_tan_stencil(fr,tan_rad);
	uswssten->nor_ans = alloc_tan_stencil(fr,tan_rad);
	sizest = fr->sizest;
	intfc = fr->interf;
	for (i = 1; i <= nts/2; ++i)
	{
	    for (j = -tan_rad; j <= tan_rad ; ++j)
	    {
	        obstacle_state(intfc,left_state(tan_fr(uswssten)[-i]->p[j]),
			       sizest);
	        obstacle_state(intfc,right_state(tan_fr(uswssten)[i]->p[j]),
			       sizest);
	    }
	}
	for (i = -nns/2; i <= nns/2; ++i)
	    nor_fr(uswssten)[i] = AllocWSSten(nor_rad,stype_WSSten(fr),fr);
	uswssten->tan_ans = AllocWSSten(nor_rad,stype_WSSten(fr),fr);
	bi_array(&nor_vec(uswssten),nns,MAXD,FLOAT);
	nor_vec(uswssten) += nns/2;
	return uswssten;
}	/*end AllocUSWSSten2d*/


EXPORT	WSSten	*g_AllocWSSten(
	int	nsts,
	int	stype,
	Front	*front)
{
	int	i;
	size_t	sizest = front->sizest;
	WSSten	*sten;

	scalar(&sten,sizeof(WSSten));
	sten->nsts = nsts;
	sten->coords = sten->coords_store;
	sten->nor = sten->nor_store;
	sten->_ws_interpolate = g_ws_interpolate;
	sten->_set_ws_slopes = g_set_ws_slopes;
	sten->_ClearWSStenData = g_ClearWSStenData;
	sten->_FPrintWSSten = g_FPrintWSSten;
	sten->ncomp = NO_COMP;
	sten->pcomp = NO_COMP;
	sten->w_type = UNKNOWN_WAVE_TYPE;
	sten->stype = stype;

	uni_array(&sten->lcrds,nsts,sizeof(float*));
	bi_array(&sten->lcrds_store,nsts,3,sizeof(float));

	uni_array(&sten->rcrds,nsts,sizeof(float*));
	bi_array(&sten->rcrds_store,nsts,3,sizeof(float));

	uni_array(&sten->sl,nsts,sizeof(Locstate));
	uni_array(&sten->sl_store,nsts,sizest);

	uni_array(&sten->dsl,nsts,sizeof(Locstate));
	uni_array(&sten->dsl_store,nsts,sizest);

	uni_array(&sten->tsl,nsts,sizeof(Locstate));
	uni_array(&sten->tsl_store,nsts,sizest);

	uni_array(&sten->sr,nsts,sizeof(Locstate));
	uni_array(&sten->sr_store,nsts,sizest);

	uni_array(&sten->dsr,nsts,sizeof(Locstate));
	uni_array(&sten->dsr_store,nsts,sizest);

	uni_array(&sten->tsr,nsts,sizeof(Locstate));
	uni_array(&sten->tsr_store,nsts,sizest);

	for (i = 0; i < nsts; ++i)
	{
	    sten->lcrds[i] = sten->lcrds_store[i];
	    sten->rcrds[i] = sten->rcrds_store[i];

	    sten->sl[i] = (Locstate) (sten->sl_store + i*sizest);
	    zero_scalar(sten->sl[i],sizest);
	    set_type_of_state(sten->sl[i],UNKNOWN_STATE);

	    sten->dsl[i] = (Locstate) (sten->dsl_store + i*sizest);
	    zero_scalar(sten->dsl[i],sizest);
	    set_type_of_state(sten->dsl[i],UNKNOWN_STATE);

	    sten->tsl[i] = (Locstate) (sten->dsl_store + i*sizest);
	    zero_scalar(sten->tsl[i],sizest);
	    set_type_of_state(sten->tsl[i],UNKNOWN_STATE);

	    sten->sr[i] = (Locstate) (sten->sr_store + i*sizest);
	    zero_scalar(sten->sr[i],sizest);
	    set_type_of_state(sten->sr[i],UNKNOWN_STATE);

	    sten->dsr[i] = (Locstate) (sten->dsr_store + i*sizest);
	    zero_scalar(sten->dsr[i],sizest);
	    set_type_of_state(sten->dsr[i],UNKNOWN_STATE);

	    sten->tsr[i] = (Locstate) (sten->tsr_store + i*sizest);
	    zero_scalar(sten->tsr[i],sizest);
	    set_type_of_state(sten->tsr[i],UNKNOWN_STATE);
	}

	return sten;
}	/*end g_AllocWSSten*/

/*
*			g_set_ws_slopes():
*
*	Main driver for the setting of interpolation slopes in the npt wave
*	speed stencil.
*/

LOCAL	void	g_set_ws_slopes(
	WSSten	*wssten)
{
	int		i, nsts = wssten->nsts;
	int		st_type = wssten->stype;
	size_t		sizest = wssten->front->sizest;
	Locstate	*sl = wssten->sl;
	Locstate	*dsl = wssten->dsl;
	Locstate	*tsl = wssten->tsl;
	Locstate	*sr = wssten->sr;
	Locstate	*dsr = wssten->dsr;
	Locstate	*tsr = wssten->tsr;


	for (i = 0; i < nsts; ++i)
	{
	    tsl[i] = (Locstate) (wssten->tsl_store + i*sizest);
	    set_state(tsl[i],st_type,sl[i]);
	    tsr[i] = (Locstate) (wssten->tsr_store + i*sizest);
	    set_state(tsr[i],st_type,sr[i]);
	}
	if (st_type == OBSTACLE_STATE)
	{
	    for (i = 0; i < nsts; ++i)
	    {
	        g_obstacle_state(dsl[i],sizest);
	        g_obstacle_state(dsr[i],sizest);
	    }
	    return;
	}
	if (is_obstacle_state(tsl[0]))
	{
	    for (i = 0; i < nsts; ++i)
	        g_obstacle_state(dsl[i],sizest);
	}
	else
	    set_slopes_of_states(-1.0,tsl,dsl,wssten);

	if (is_obstacle_state(sr[0]))
	{
	    for (i = 0; i < nsts; ++i)
	        g_obstacle_state(dsr[i],sizest);
	}
	else
	    set_slopes_of_states(1.0,tsr,dsr,wssten);

}	/*end g_set_ws_slopes*/

/*
*			set_slopes_of_states():
*
*	Compute the interpolation slopes for a npt wave speed stencil.
*	Slopes are computed using the Van Leer limiter.
*/

LOCAL	void	set_slopes_of_states(
	float		dir,
	Locstate	*s,
	Locstate	*ds,
	WSSten		*wssten)
{
	int	dim, n;
	int	i, nsts = wssten->nsts;
	int	st_type = wssten->stype;

	dim = Params(s[0])->dim;
	switch (nsts)
	{
	case 1:
	    g_obstacle_state(ds[0],Params(s[0])->sizest);/*Zero all fields*/
	    set_type_of_state(ds[0],st_type);
	    break;
	case 2:
	    set_type_of_state(ds[0],st_type);
	    set_type_of_state(ds[1],st_type);
	    Dens(ds[0]) = Dens(ds[1]) = Dens(s[1]) - Dens(s[0]);
            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
                int   k;
                if(Params(s[1])->n_comps != 1 ||
                   Params(s[0])->n_comps != 1)
                {
                    if(Params(s[0])->n_comps != Params(s[1])->n_comps)
                    {
                        printf("ERROR: set_slopes_of_states\n");
                        printf("MULTI_COMP_NON_REACTIVE, s[0] and s[1] have different # of comps.\n");
                        clean_up(ERROR);
                    }
                    for(k = 0; k < Params(s[0])->n_comps; k++)
                        pdens(ds[0])[k] = pdens(ds[1])[k] = pdens(s[1])[k] - pdens(s[0])[k];
                }
            }
	    Press(ds[0]) = Press(ds[1]) =  Press(s[1]) - Press(s[0]);
	    for (i = 0; i < dim; ++i)
	        Vel(ds[0])[i] = Vel(ds[1])[i] = Vel(s[1])[i] - Vel(s[0])[i];
	    reset_gamma(ds[0]);
#if defined(COMBUSTION_CODE)
	    switch(Composition_type(s[0]))
	    {
	    case ZND:
	        React(ds[0]) = React(ds[1])= React(s[1]) - React(s[0]);
	        break;
	    case TWO_CONSTITUENT_REACTIVE:
	        React(ds[0]) = React(ds[1])= React(s[1]) - React(s[0]);
	        Dens1(ds[0]) = Dens1(ds[1])= Dens1(s[1]) - Dens1(s[0]);
	        break;
	    case PURE_NON_REACTIVE:
	    case PTFLAME:
	    case THINFLAME:
	    default:
	        break;
	    }
#endif /* defined(COMBUSTION_CODE) */
	    break;
	default:
	    n = nsts-1;
	    for (i = 1; i < n; ++i)
	        set_state_slopes(ds[i],s+i);
	    set_end_state_slopes(1.0,s[0],ds[0],s[1],ds[1],wssten);
	    set_end_state_slopes(-1.0,s[n],ds[n],s[n-1],ds[n-1],wssten);
	    break;
	}
	if ((dim > 1) && is_scalar_wave(wssten->w_type) && wssten->p)
	{
	    /*
	    float curvature =
	      mean_curvature_at_point(wssten->p,wssten->hse,wssten->hs,
	                              wssten->front);
	    */
	    float curvature = (wssten->p)->curvature;
	    float dh = mag_vector(wssten->front->rect_grid->h,dim);
	    dir *= 1.0/(1.0 + curvature_factor(wssten->front)*curvature*dh);
	}
	for (n = 0; n < nsts; ++n)
	{
	    if(debugging("tst_nan"))
	    {
	        printf("#n bf %d\n", n);
		fprint_raw_gas_data(stdout, ds[n], dim);
		//verbose_print_state("ds", ds[n]);
	    }

	    Dens(ds[n]) *= dir;
            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
                int   k;
                if(Params(s[0])->n_comps != 1)
                {
                    for(k = 0; k < Params(s[0])->n_comps; k++)
                        pdens(ds[n])[k] *= dir;
                }
            }
	    Press(ds[n]) *= dir;
	    for (i = 0; i < dim; ++i)
	        Vel(ds[n])[i] *= dir;
	    reset_gamma(ds[n]);
#if defined(COMBUSTION_CODE)
	    switch(Composition_type(s[0]))
	    {
	    case ZND:
	        React(ds[n]) *= dir;
	        break;
	    case TWO_CONSTITUENT_REACTIVE:
	        React(ds[n]) *= dir;
	        Dens1(ds[n]) *= dir;
	        break;
	    case PURE_NON_REACTIVE:
	    case PTFLAME:
	    case THINFLAME:
	    default:
	        break;
	    }
#endif /* defined(COMBUSTION_CODE) */

	    if(debugging("tst_nan"))
	    {
	        printf("#n af %d %24.16e\n", n, dir);
		fprint_raw_gas_data(stdout, ds[n], dim);
		//verbose_print_state("ds", ds[n]);
	    }
	}
}	/*end set_slopes_of_states*/

/*
*			set_end_state_slopes():
*
*	Set the state slopes at the endpoint of a wave speed stencil.
*	The value of sgn is 1.0 at the front end,  and -1.0 at the far
*	field end.
*/

/*ARGSUSED*/
LOCAL	void	set_end_state_slopes(
	float		sgn,
	Locstate	se,	/*End point state*/
	Locstate	dse,	/*End point slope*/
	Locstate	si,	/*Interior state*/
	Locstate	dsi,	/*Interior slope*/
	WSSten		*wssten)
{
	int	i, dim;

	dim = Params(si)->dim;

	/*
	 *   This version uses increased upwinding to set the slope at the
	 *   end point of the npt wave speed stencil.  The state at the end
	 *   is replaced by the extrapolated value from the interior and the
	 *   slope at the interior.  The slope is copied from interior to
	 *   the endpoint.
	 */

	Dens(dsi) = endpt_limiter(Dens(se),Dens(si),Dens(dsi));
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            int   k;
            /* The density slopes of components are not consistent */
            if(Params(se)->n_comps != 1 ||
               Params(si)->n_comps != 1)
            {
                if(Params(se)->n_comps != Params(si)->n_comps)
                {
                    printf("ERROR: set_end_state_slopes\n");
                    printf("se and si have different # of comps.\n");
                    clean_up(ERROR);
                }
                for(k = 0; k < Params(se)->n_comps; k++)
                    pdens(dsi)[k] = endpt_limiter(pdens(se)[k],pdens(si)[k],pdens(dsi)[k]);
            }
        }
	Press(dsi) = endpt_limiter(Press(se),Press(si),Press(dsi));
	for (i = 0; i < dim; ++i)
	    Vel(dsi)[i] = endpt_limiter(Vel(se)[i],Vel(si)[i],Vel(dsi)[i]);
	reset_gamma(dsi);
#if defined(COMBUSTION_CODE)
	switch(Composition_type(si))
	{
	case ZND:
	    React(dsi) = endpt_limiter(React(se),React(si),React(dsi));
	    break;
	case TWO_CONSTITUENT_REACTIVE:
	    React(dsi) = endpt_limiter(React(se),React(si),React(dsi));
	    Dens1(dsi) = endpt_limiter(Dens1(se),Dens1(si),Dens1(dsi));
	    break;
	case PURE_NON_REACTIVE:
	case PTFLAME:
	case THINFLAME:
	default:
	    break;
	}
#endif /* defined(COMBUSTION_CODE) */

	Dens(se) = Dens(si) - sgn*Dens(dsi);
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            int    k;
            float  sum;
            /* scaling pdens(se) to be consistent */
            if(Params(se)->n_comps != 1)
            {
                sum = 0.0;
                for(k = 0; k < Params(se)->n_comps; k++)
                    pdens(se)[k] = pdens(si)[k] - sgn*pdens(dsi)[k];
                for(k = 0; k < Params(se)->n_comps; k++)
                    sum += pdens(se)[k];
                /**
                printf("IN set_end_state_slopes, before scaling\n");
                printf("se pdens[0] %g, pdens[0] %g, pdens[0] %g, dens %g, diff %g\n",
                     pdens(se)[0], pdens(se)[1], pdens(se)[2], Dens(se),
                     pdens(se)[0]+ pdens(se)[1]+ pdens(se)[2]- Dens(se));
                printf("si pdens[0] %g, pdens[0] %g, pdens[0] %g, dens %g, diff %g\n",
                     pdens(si)[0], pdens(si)[1], pdens(si)[2], Dens(si),
                     pdens(si)[0]+ pdens(si)[1]+ pdens(si)[2]- Dens(si));
                **/

                for(k = 0; k < Params(se)->n_comps; k++)
                    pdens(se)[k] = Dens(se)/sum*pdens(se)[k];

                /**
                printf("IN set_end_state_slopes\n");
                printf("se pdens[0] %g, pdens[0] %g, pdens[0] %g, dens %g, diff %g\n\n",
                     pdens(se)[0], pdens(se)[1], pdens(se)[2], Dens(se),
                     pdens(se)[0]+ pdens(se)[1]+ pdens(se)[2]- Dens(se));
                **/
            }
        }
	Press(se) = Press(si) - sgn*Press(dsi);
	for (i = 0; i < dim; ++i)
	    Vel(se)[i] = Vel(si)[i] - sgn*Vel(dsi)[i];
	set_type_of_state(se,state_type(si));
	reset_gamma(si);
#if defined(COMBUSTION_CODE)
	switch(Composition_type(si))
	{
	case ZND:
	    React(se) = React(si) - sgn*React(dsi);
	    break;
	case TWO_CONSTITUENT_REACTIVE:
	    React(se) = React(si) - sgn*React(dsi);
	    Dens1(se) = Dens1(si) - sgn*Dens1(dsi);
	    break;
	case PURE_NON_REACTIVE:
	case PTFLAME:
	case THINFLAME:
	default:
	    break;
	}
#endif /* defined(COMBUSTION_CODE) */
	ft_assign(dse,dsi,Params(si)->sizest);
}		/*end set_end_state_slopes*/



LOCAL	float endpt_limiter(
	float	xe,
	float	xi,
	float	dx)
{
	float	fdie = fabs(xi - xe);
	float	fdx = fabs(dx);
	float	smin = min(fdx,fdie);

	return (dx < 0.0) ? -smin : smin;
}		/*end endpt_limiter*/


/*
*			set_state_slopes():
*
*	Uses the van Leer limiter to compute the slope of a state at an
*	interior point of a wave speed stencil.
*/

LOCAL	void	set_state_slopes(
	Locstate	slope,
	Locstate	*s)
{
	Locstate	sl, s0, sr;
	int i, dim;

	sl = s[-1];
	s0 = s[0];
	sr = s[1];

	set_type_of_state(slope,state_type(s0));
	Dens(slope) = van_leer_limiter(Dens(sl),Dens(s0),Dens(sr));
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            int   k;
            if(Params(sl)->n_comps != 1 ||
               Params(sr)->n_comps != 1 ||
               Params(s0)->n_comps != 1)
            {
                if(Params(s0)->n_comps != Params(sl)->n_comps ||
                   Params(s0)->n_comps != Params(sr)->n_comps ||
                   Params(sl)->n_comps != Params(sr)->n_comps)
                {
                    printf("ERROR: set_state_slopes\n");
                    printf(" states have different # of comps.\n");
                    clean_up(ERROR);
                }
                for(k = 0; k < Params(s0)->n_comps; k++)
                    pdens(slope)[k] = van_leer_limiter(pdens(sl)[k],pdens(s0)[k],pdens(sr)[k]);
            }
        }
	Press(slope) = van_leer_limiter(Press(sl),Press(s0),Press(sr));
	dim = Params(s[0])->dim;
	for (i = 0; i < dim; ++i)
	    Vel(slope)[i] = van_leer_limiter(Vel(sl)[i],Vel(s0)[i],Vel(sr)[i]);
	reset_gamma(slope);
#if defined(COMBUSTION_CODE)
	switch(Composition_type(s0))
	{
	case ZND:
	    React(slope) = van_leer_limiter(React(sl),React(s0),React(sr));
	    break;
	case TWO_CONSTITUENT_REACTIVE:
	    React(slope) = van_leer_limiter(React(sl),React(s0),React(sr));
	    Dens1(slope) = van_leer_limiter(Dens1(sl),Dens1(s0),Dens1(sr));
	    break;
	case PURE_NON_REACTIVE:
	case PTFLAME:
	case THINFLAME:
	default:
	    break;
	}
#endif /* defined(COMBUSTION_CODE) */
}	/*end set_state_slopes*/


LOCAL	float	van_leer_limiter(
	float	sl,
	float	s0,
	float	sr)
{
	float	sgn;
	float	ds, dsl, dsr, ds0;

	dsl = 2.0*(s0 - sl);
	ds0 = 0.5*(sr - sl);
	dsr = 2.0*(sr - s0);
	if (dsl*dsr <= 0.0)
	    return 0.0;
	sgn = (ds0 >= 0.0) ? 1.0 : -1.0;
	dsl = fabs(dsl);
	dsr = fabs(dsr);
	ds0 = fabs(ds0);
	ds = min(dsl,dsr);
	ds = min(ds,ds0);
	return sgn*ds;
}	/*end van_leer_limiter*/
	

/*
*			g_ws_interpolate():
*
*	Interpolates values near a front using the state data provided by the
*	WSSten structure.
*
*	Input:
*		x - fractional distance between the front and the off front
*		    location.  The interpolation is done at the position
*		    (1 - x)*crds0 + x*crds1.  The coordinate system runs
*		    in the direction of the front normal.  Thus for negative
*		    sides x = 0 corresponds to the off front location while
*                   x = 1 corresponds to the on front location.  For positive
*		    side x = 0 corresponds to the on front location while
*	            x = 1 is the off front location.  The interpolation is
*		    based on the formula
*
*				f(x) = f + x*df
*
*		    Note that when computing the interpolation using the state
*		    slopes,  this requires that x be replaced by 1 - x for the
*		    negative side.
*
*		side - POSITIVE_SIDE or NEGATIVE_SIDE depending on the side of
*                      of the interface upon which the interpolation is done.
*
*		stype - the state type representation of the returned answer
*                       state.
*
*		wssten - the normal sweep stencil containing the state data
*	                 and reconstructed slopes
*
*	Output:
*		ans - the interpolated state at position (1-x)*crds0 + x*crds1.
*/

/*ARGSUSED*/
LOCAL	void	g_ws_interpolate(
	Locstate	ans,
	float		x,
	SIDE            side,
	int             stype,
	WSSten		*wssten)
{
	float	         *crds0;
	Locstate         st0;
	Locstate         slope0;
	float	         *crds1;
	Locstate         st1;
	Locstate         slope1;
	Locstate         st, slope;
	int	         i, dim;
	int	         st_type;
	static	Locstate wrkst = NULL;

	if (side == POSITIVE_SIDE)
	{
	    crds0  = wssten->coords;
	    st0    = wssten->tsr[0];
	    slope0 =  wssten->dsr[0];
	    crds1  = wssten->rcrds[1];
	    st1    = wssten->tsr[1];
	    slope1 = wssten->dsr[1];

	}
	else if (side == NEGATIVE_SIDE)
	{
	    crds0  =  wssten->lcrds[1];
	    st0    = wssten->tsl[1];
	    slope0 = wssten->dsl[1];
	    crds1  = wssten->coords;
	    st1    = wssten->tsl[0];
	    slope1 = wssten->dsl[0];

	}
	else
	{
	    screen("ERROR in g_ws_interpolate(), invalid side %d\n",side);
	    clean_up(ERROR);
	}

	if(debugging("tst_nan"))
	//if (debugging("interpolate"))
	{
	    dim = wssten->front->rect_grid->dim;
	    (void) printf("Data into g_ws_interpolate\n");
	    (void) printf("x = %g, side = %s, stype = %s\n",
			  x,side_name(side),state_type_name(stype));
	    print_general_vector("crds0 = ",crds0,dim,"\n");
	    verbose_print_state("st0",st0);
	    if (slope0 != NULL)
	    {
	        (void) printf("slope0\n");
	        fprint_raw_gas_data(stdout,slope0,dim);
	    }
	    else
	        (void) printf("slope0 = NULL\n");
	    print_general_vector("crds1 = ",crds1,dim,"\n");
	    verbose_print_state("st1",st1);
	    if (slope1 != NULL)
	    {
	        (void) printf("slope1\n");
	        fprint_raw_gas_data(stdout,slope1,dim);
	    }
	    else
	        (void) printf("slope1 = NULL\n");
	}
	if ((slope0 == NULL) || (slope1 == NULL) || (wssten->nsts <= 2))
	{
	    if (debugging("interpolate"))
		(void) printf("g_ws_interpolate using linear interpolation\n");
	    interpolate_states(wssten->front,1.0-x,x,crds0,st0,crds1,st1,ans);
	    if (debugging("interpolate"))
	        verbose_print_state("ans",ans);
	    return;
	}

 	if (x < 0.5)
    	{
    	    st = st0;
    	    slope = slope0;
    	}
    	else
    	{
    	    x = x - 1.0;
    	    st = st1;
    	    slope = slope1;
    	}

	if (is_obstacle_state(st))
	{
	    g_obstacle_state(ans,g_sizest());
	    return;
	}
	(*Params(st)->_clear_state)(ans,Params(st)->sizest);
	Set_params(ans,st);
	st_type = state_type(st);
	if (st_type != state_type(slope))
	{
	    if (wrkst == NULL)
		scalar(&wrkst,Params(st)->sizest);
	    set_state(wrkst,state_type(slope),st);
	    st = wrkst;
	}
	set_type_of_state(ans,state_type(st));

	/*
	*  NOTE: The following code is independent of the state type
	*  since the lines simply represent linear interpolation of
	*  floating point fields in the state.
	*/

	Dens(ans) = Dens(st) + x*Dens(slope);
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            int    k;
            float  sum = 0.0;
            // printf("In g_ws_interpolate\n");
            // verbose_print_state("st", st);
            // verbose_print_state("slope", slope);
            /* scaling partial density */
            if(Params(ans)->n_comps != 1)
            {
                for(k = 0; k < Params(ans)->n_comps; k++)
                {
                    pdens(ans)[k] = pdens(st)[k] + x*pdens(slope)[k];
                    /* New 051005 */
                    if(fabs(pdens(ans)[k]) < 10.0*MACH_EPS && pdens(ans)[k] < 0.0)
                        pdens(ans)[k] = 0.0;
  //                  else if(fabs(pdens(ans)[k]) > 10.0*MACH_EPS && pdens(ans)[k] < 0.0)
  //                  {
    //                    printf("ERROR in g_ws_interpolate() x = %g\n", x);
      //                  printf("partial density < 0.0\n");
        //                clean_up(ERROR);
          //          }
                    /* End of New 051005 */
                }
                for(k = 0; k < Params(ans)->n_comps; k++)
                    sum += pdens(ans)[k];
                for(k = 0; k < Params(ans)->n_comps; k++)
                    pdens(ans)[k] = Dens(ans)/sum*pdens(ans)[k];

                /***
                printf("IN  g_ws_interpolate() x = %g\n", x);
                printf("ans: pdens[0] %g, pdens[1] %g, pdens[2] %g, dens %g, diff %g\n",
                     pdens(ans)[0], pdens(ans)[1], pdens(ans)[2], Dens(ans),
                     pdens(ans)[0]+ pdens(ans)[1]+ pdens(ans)[2]- Dens(ans));
                printf("st: pdens[0] %g, pdens[1] %g, pdens[2] %g, dens %g, diff %g\n",
                     pdens(st)[0], pdens(st)[1], pdens(st)[2], Dens(st),
                     pdens(st)[0]+ pdens(st)[1]+ pdens(st)[2]- Dens(st));
                printf("slope: pdens[0] %g, pdens[1] %g, pdens[2] %g, dens %g, diff %g\n",
                     pdens(slope)[0], pdens(slope)[1], pdens(slope)[2], Dens(slope),
                     pdens(slope)[0]+ pdens(slope)[1]+ pdens(slope)[2]- Dens(slope));
                ***/
            }
        }
	Press(ans) = Press(st) + x*Press(slope);
	dim = Params(st)->dim;
	for (i = 0; i < dim; ++i)
	    Vel(ans)[i] = Vel(st)[i] + x*Vel(slope)[i];
	reset_gamma(ans);
#if defined(COMBUSTION_CODE)
	switch(Composition_type(ans))
	{
	case ZND:
	    React(ans) = React(st) + x*React(slope);
	    break;
	case TWO_CONSTITUENT_REACTIVE:
	    React(ans) = React(st) + x*React(slope);
	    Dens1(ans) = Dens1(st) + x*Dens1(st);
	    break;
	case PURE_NON_REACTIVE:
	case PTFLAME:
	case THINFLAME:
	default:
	    break;
	}
#endif /* defined(COMBUSTION_CODE) */

	set_state(ans,stype,ans);

	if (debugging("interpolate"))
	    verbose_print_state("ans",ans);
}	/*end g_ws_interpolate*/




/*
*			states_near_location():
*
*	Find states near front for normal sweep.
*/

EXPORT	void states_near_location(
	WSSten		*sten,
	float		*coords,
	float		*nor,
	COMPONENT	pcomp,
	COMPONENT	ncomp,
	Locstate	sl,
	Locstate	sr)
{
        printf("enter states_near_location\n");
       	Front		*fr = sten->front;
	HYPER_SURF	*hs = sten->hs;
	Wave		*wave = sten->wave;
	Locstate	stmpl, stmpr;
	float		dn = sten->dn;
	float		*h = fr->rect_grid->h;
	int		dim = fr->interf->dim;
	int		indxl, indxr, i, j, nsts = sten->nsts;

#if defined(DEBUG_POINT_PROPAGATE)
	debug_print("point_propagate","Entered states_near_location()\n");
#endif /* defined(DEBUG_POINT_PROPAGATE) */

printf("xiaoxue_statesnear\n");
	sten->pcomp = pcomp;
	sten->ncomp = ncomp;
	if (nor != sten->nor)
	{
	    for (j = 0; j < dim; ++j)
		sten->nor[j] = nor[j];
	}
printf("xiaoxue_statesnear\n");
	for (j = 0; j < dim; ++j)
	    sten->coords[j] = sten->lcrds[0][j] = sten->rcrds[0][j] = coords[j];
printf("xiaoxue_statesnear\n");
	sten->dn = dn = grid_size_in_direction(sten->nor,h,dim);
	for (i = 1; i < nsts; ++i)
	{
	    for (j = 0; j < dim; ++j)
	    {
	        sten->lcrds[i][j] = coords[j] - i*nor[j]*dn;
	        sten->rcrds[i][j] = coords[j] + i*nor[j]*dn;
	    }
	}
printf("xiaoxue_statesnear\n");	
	set_state(sten->sl[0],state_type(sl),sl);
	set_state(sten->sr[0],state_type(sr),sr);
printf("xiaoxue_statesnear\n");
	interface_crosses_stencil(sten,&indxl,&indxr,&stmpl,&stmpr);
printf("xiaoxue_statesnear\n");
	sten->indxl = indxl;
	sten->indxr = indxr;

	if (is_obstacle_state(sl))
	{
        printf("xiaoxue_statesnear_isobs\n");
#if defined(DEBUG_POINT_PROPAGATE)
	    if (debug_point_propagate)
	        (void) printf("sl is an obstacle state\n");
#endif /* defined(DEBUG_POINT_PROPAGATE) */
	    for (i = 0; i < nsts; ++i)
	        obstacle_state(fr->interf,sten->sl[i],fr->sizest);
	}
	else
	{
                printf("xiaoxue_statesnear not obstacle\n");
	    for (i = 1; i < indxl; ++i)
	    {
#if defined(DEBUG_POINT_PROPAGATE)
	        if (debug_point_propagate)
		{
	            (void) printf("Using hyp_solution to "
		                  "evaluate sten->sl[%d]\n",i);
		    (void) printf("sten->lcrds[%d] = ",i);
		    print_general_vector("",sten->lcrds[i],dim,", ");
		    (void) printf("sten->ncomp = %d\n",sten->ncomp);
	        }
#endif /* defined(DEBUG_POINT_PROPAGATE) */
	        hyp_solution(sten->lcrds[i],sten->ncomp,hs,NEGATIVE_SIDE,
		             fr,wave,sten->sl[i],sl);
	    }
	    for (i = indxl; i < nsts; ++i)
	        set_state(sten->sl[i],state_type(sl),stmpl);
	}
	if (is_obstacle_state(sr))
	{
                printf("xiaoxue_statesnear is_obstacle_state(sr)\n");
#if defined(DEBUG_POINT_PROPAGATE)
	    if (debug_point_propagate)
	        (void) printf("sr is an obstacle state\n");
#endif /* defined(DEBUG_POINT_PROPAGATE) */
	    for (i = 0; i < nsts; ++i)
	        obstacle_state(fr->interf,sten->sr[i],fr->sizest);
	}
	else
	{
                printf("xiaoxue_statesnear sr not obstale\n");
	    for (i = 1; i < indxr; ++i)
	    {
#if defined(DEBUG_POINT_PROPAGATE)
	        if (debug_point_propagate)
		{
	            (void) printf("Using hyp_solution to "
		                  "evaluate sten->sr[%d]\n",i);
		    (void) printf("sten->rcrds[%d] = ",i);
		    print_general_vector("",sten->rcrds[i],dim,", ");
		    (void) printf("sten->pcomp = %d\n",sten->pcomp);
	        }
#endif /* defined(DEBUG_POINT_PROPAGATE) */
	        hyp_solution(sten->rcrds[i],sten->pcomp,hs,POSITIVE_SIDE,
			     fr,wave,sten->sr[i],sr);
	    }
	    for (i = indxr; i < nsts; ++i)
	        set_state(sten->sr[i],state_type(sr),stmpr);
	}
	if (debugging("bad_state"))
	{
	    bool bad = NO;
	    for (i = 0; i < nsts; ++i)
	    {
	        if (is_bad_state(sten->sr[i],YES,"states_near_location"))
		{
		    bad = YES;
		    (void) printf("sten->sr[%d] is bad\n",i);
	            fprint_raw_gas_data(stdout,sten->sr[i],dim);
		}
	        if (is_bad_state(sten->sl[i],YES,"states_near_location"))
		{
		    bad = YES;
		    (void) printf("sten->sl[%d] is bad\n",i);
	            fprint_raw_gas_data(stdout,sten->sl[i],dim);
		}
	    }
	    if (bad)
	    {
		print_WSStenData(sten);
	        screen("ERROR in states_near_location(), bad state detected\n");
		clean_up(ERROR);
	    }
	}
printf("leave xiaoxue_statesnear\n");
#if defined(DEBUG_POINT_PROPAGATE)
	debug_print("point_propagate","Left states_near_location()\n");
#endif /* defined(DEBUG_POINT_PROPAGATE) */
}		/*end states_near_location*/


#if defined(DEBUG_POINT_PROPAGATE)
LOCAL	void entered_point_propagate(
	char		*fname,
	HYPER_SURF	*oldhs)
{
	int		w_type;

	debug_print("point_propagate","\nEntered %s()\n",fname);

	w_type = wave_type(oldhs);
	switch(w_type) 
	{
	case PASSIVE_BOUNDARY:
	    break;
	case NEUMANN_BOUNDARY:
	    if (debugging("nbdry")) debug_point_propagate = YES;
	    break;
	case MOVABLE_NEUMANN_BOUNDARY:
	    if (debugging("mnbdry")) debug_point_propagate = YES;
	    break;
	case DIRICHLET_BOUNDARY:
	    if (debugging("dbdry")) debug_point_propagate = YES;
	    break;
	case TIME_DIRICHLET_BOUNDARY:
	    if (debugging("tdbdry")) debug_point_propagate = YES;
	    break;
	case VELOCITY_SPECIFIED:
	    if (debugging("vsbdry")) debug_point_propagate = YES;
	    break;
	case SUBDOMAIN_BOUNDARY:
	    if (pp_numnodes() > 1 && debugging("sdbdry"))
	    	debug_point_propagate = YES;
	    else if (debugging("pbdry")) /*TODO REMOVE*/
	    	debug_point_propagate = YES;
	    break;
	case BACKWARD_SHOCK_WAVE:
	    if (debugging("shock") || debugging("bshock"))
	    	debug_point_propagate = YES;
	    break;
	case BACKWARD_SOUND_WAVE_LE:
	    if (debugging("raref") || debugging("le_raref") ||
						debugging("ble_raref"))
		debug_point_propagate = YES;
	    break;
	case BACKWARD_SOUND_WAVE_TE:
	    if (debugging("raref") || debugging("te_raref") ||
						debugging("bte_raref"))
	    	debug_point_propagate = YES;
	    break;
	case CONTACT:
	    if (debugging("contact"))
	    	debug_point_propagate = YES;
	    break;
	case FORWARD_SHOCK_WAVE:
	    if (debugging("shock") || debugging("fshock"))
	    	debug_point_propagate = YES;
	    break;
	case FORWARD_SOUND_WAVE_LE:
	    if (debugging("raref") || debugging("le_raref") ||
						debugging("fle_raref"))
	    	debug_point_propagate = YES;
	    break;
	case FORWARD_SOUND_WAVE_TE:
	    if (debugging("raref") || debugging("te_raref") ||
						debugging("fte_raref"))
	    	debug_point_propagate = YES;
	    break;
	default:
	    break;
	}
	if (debugging("point_propagate"))
	    debug_point_propagate = YES;
	else if (debug_point_propagate == YES)
	    (void) printf("\nEntered %s()\n",fname);
}		/*end entered_point_propagate*/

LOCAL	void left_point_propagate(
	char   *fname,
	POINT  *newp,
	float  *V,
	float  dt,
	int    dim,
	WSSten *sten)
{
#if defined(CHECK_FOR_BAD_STATES)
	if (debugging("bad_state") &&
	    (is_bad_state(left_state(newp),YES,"left_point_propagate") ||
	     is_bad_state(right_state(newp),YES,"left_point_propagate")))
	{
	    screen("ERROR in left_point_propagate(), bad state detected\n");
	    (void) printf("left_state(newp) - ");
	    fprint_raw_gas_data(stdout,left_state(newp),dim);
	    (void) printf("right_state(newp) - ");
	    fprint_raw_gas_data(stdout,right_state(newp),dim);
	    if (sten)
	    {
		print_point_propagate_data(sten->p,sten->hse,sten->hs,dim);
	        print_WSStenData(sten);
	    }
	    clean_up(ERROR);
	}
#endif /* defined(CHECK_FOR_BAD_STATES) */

	if (debug_point_propagate == NO)
	    return;

	print_general_vector("velocity = ",V,dim,"\n");
	(void) printf("dt = %g\n",dt);
	print_general_vector("new position = ",Coords(newp),dim,"\n");
	verbose_print_state("answer state left",left_state(newp));
	verbose_print_state("answer state right",right_state(newp));

	debug_print("point_propagate","\nLeft %s()\n",fname);
	if (debugging("point_propagate"))
	    return;
	(void) printf("\nLeft %s()\n",fname);

	debug_point_propagate = NO;
}		/*end left_point_propagate*/
#endif /* defined(DEBUG_POINT_PROPAGATE) */


EXPORT	void print_point_propagate_data(
	POINT		   *p,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF	   *hs,
	int		   dim)
{
	(void) printf(" hs = %llu, hse = %llu, Boundary_hs(hs) = %d\n",
		      hypersurface_number(hs),
		      hypersurface_element_number(hse,hs->interface),
		      Boundary_hs(hs));
	print_wave_type(" wave type(hs) = ",wave_type(hs),"\n",hs->interface);
	(void) printf(" positive component = %d, negative_component = %d\n",
	              positive_component(hs),negative_component(hs));
	switch (dim)
	{
	case 1:
	    (void) printf("propagate point 0x%p: %g\n",(POINTER)p,Coords(p)[0]);
	    break;
	case 2:
	{
	    CURVE	*c = Curve_of_hs(hs);
	    BOND	*b = Bond_of_hse(hse);

	    if (p == c->end->posn) 
	    {
	    	(void) printf("Propagating node of curve %llu\n",
			      curve_number(c));
	        (void) printf("at end of bond %llu:",
			      bond_number(c->last,c->interface));
	    }
	    else if (p == c->start->posn) 
	    {
	    	(void) printf("Propagating node of curve %llu\n",
			      curve_number(c));
	    	(void) printf("at start of bond %llu:",
			      bond_number(c->first,c->interface));
	    }
	    else if (p == b->start) 
	    {
	    	(void) printf("Propagating point on curve %llu\n",
			      curve_number(c));
	       	(void) printf("at start of bond %llu:",
			      bond_number(b,c->interface));
	    }
	    else 
	    {
	    	(void) printf("Propagating point on curve %llu\n",
			      curve_number(c));
		(void) printf("at end of bond %llu:",
			      bond_number(b,c->interface));
	    }
	    (void) printf(" %g %g -> %g %g\n",
			  Coords(b->start)[0],Coords(b->start)[1],
			  Coords(b->end)[0],Coords(b->end)[1]);
	    break;
	}
	case 3:
	{
	    (void) printf("propagate point %llu: %g %g %g\n", point_number(p),
			  Coords(p)[0],Coords(p)[1],Coords(p)[2]);
	    print_tri(Tri_of_hse(hse),hs->interface);
	    break;
	}
	default:
	    screen("ERROR in print_point_propagate_data(), "
		   "invalid dimension %d\n",dim);
	    clean_up(ERROR);
	}
}		/*end print_point_propagate_data*/


EXPORT	void	reflect_wssten(
	WSSten *wssten,
	SIDE   side,
	Front  *front)
{
	Locstate *sr = wssten->sr, *sl = wssten->sl;
	float    *pt = wssten->coords;
	float    *nor = wssten->nor;
	float    **rcrds = wssten->rcrds, **lcrds = wssten->lcrds;
	int      nsts = wssten->nsts;
	int      dim = front->rect_grid->dim;
	float    vwall = 0.0;
	float    ds, dp[3];
	int      stypel, styper;
	int      i, j;

	sr = wssten->sr;
	sl = wssten->sl;
	switch (side)
	{
	case POSITIVE_SIDE:
	    if (include_wall_normal_velocity(front))
	    {
	        vwall = scalar_product(VelocityVector(sr[0],NULL),nor,dim);
		for (i = 0; i < nsts; ++i)
		{
		    styper = state_type(sr[i]);
		    set_state(sr[i],EGAS_STATE,sr[i]);
		    for (j = 0; j < dim; ++j)
		        Vel(sr[i])[j] -= vwall*nor[j];
		    set_state(sr[i],styper,sr[i]);
		}
	    }
	    for (i = 0; i < nsts; ++i)
	    {
		for (j = 0; j < dim; ++j)
		    dp[j] = rcrds[i][j] - pt[j];
		ds = 2.0*scalar_product(dp,nor,dim);
		for (j = 0; j < dim; ++j)
		    lcrds[i][j] = rcrds[i][j] - ds*nor[j];
		set_state(sl[i],state_type(sr[i]),sr[i]);
		reflect_state(sl[i],front->interf,rcrds[i],pt,nor);
	    }
	    break;
	case NEGATIVE_SIDE:
	    if (include_wall_normal_velocity(front))
	    {
	        vwall = scalar_product(VelocityVector(sl[0],NULL),nor,dim);
		for (i = 0; i < nsts; ++i)
		{
		    stypel = state_type(sl[i]);
		    set_state(sl[i],EGAS_STATE,sl[i]);
		    for (j = 0; j < dim; ++j)
		        Vel(sl[i])[j] -= vwall*nor[j];
		    set_state(sl[i],stypel,sl[i]);
		}
	    }
	    for (i = 0; i < nsts; ++i)
	    {
		for (j = 0; j < dim; ++j)
		    dp[j] = lcrds[i][j] - pt[j];
		ds = 2.0*scalar_product(dp,nor,dim);
		for (j = 0; j < dim; ++j)
		    rcrds[i][j] = lcrds[i][j] - ds*nor[j];
		set_state(sr[i],state_type(sl[i]),sl[i]);
		reflect_state(sr[i],front->interf,lcrds[i],pt,nor);
	    }
	    break;
	default:
	    screen("ERROR in reflect_wssten(), "
	           "invalid side %s\n",side_name(side));
	    clean_up(ERROR);
	}
	if (vwall != 0.0)
	{
	    for (i = 0; i < nsts; ++i)
	    {
		stypel = state_type(sl[i]);
		styper = state_type(sr[i]);
		set_state(sl[i],EGAS_STATE,sl[i]);
		set_state(sr[i],EGAS_STATE,sr[i]);
		for (j = 0; j < dim; ++j)
		{
		    Vel(sl[i])[j] += vwall*nor[j];
		    Vel(sr[i])[j] += vwall*nor[j];
	        }
		set_state(sr[i],styper,sr[i]);
		set_state(sl[i],stypel,sl[i]);
	    }
	}
	set_ws_slopes(wssten);
}		/*end reflect_wssten*/

LOCAL	void interface_crosses_stencil(
	WSSten   *sten,
	int      *indxl,
	int      *indxr,
	Locstate *stmpl,
	Locstate *stmpr)
{
	Front      *front = sten->front;
	INTERFACE  *intfc = front->interf;
	HYPER_SURF *hs = sten->hs;
	int        nsts = sten->nsts;
	int        indx;
	int        i;
	int        dim = front->rect_grid->dim;
	POINT      *p, **pts;
	float      *crxl, *crxr;
	static Locstate sstorel, sstorer;

	*stmpl = *stmpr = NULL;
	*indxl = *indxr = nsts;

	if (debugging("noicrsten"))
	    return;

	if (sstorel == NULL)
	{
	    alloc_state(intfc,&sstorel,front->sizest);
	    alloc_state(intfc,&sstorer,front->sizest);
	}

	switch (dim)
	{
	case 1:
	    p = Point_of_hs(hs);
	    pts = intfc->points;
	    i = index_of_pointer((POINTER*)pts,(POINTER)p);
	    if ((!is_obstacle_state(sten->sl[0]))  &&
	        (i > 0)                            &&
		(!is_subdomain_boundary(Hyper_surf(pts[i-1]))) &&
		(!is_passive_boundary(pts[i-1])))

	    {
	    	if (sten->lcrds[nsts-1][0] <= Coords(pts[i-1])[0])
		{
		    for (indx = 1; indx < nsts; ++indx)
	    	        if (sten->lcrds[indx][0] <= Coords(pts[i-1])[0])
			    break;
		    *indxl = indx;
		    *stmpl = sstorel;
		    set_state(sstorel,state_type(sten->sl[0]),
		              right_state(pts[i-1]));
		}
	    }
	    if ((!is_obstacle_state(sten->sr[0]))  &&
	        (i<(intfc->num_points-1))          &&
		(!is_subdomain_boundary(Hyper_surf(pts[i+1]))) &&
		(!is_passive_boundary(pts[i+1])))
	    {
	    	if (Coords(pts[i+1])[0] <= sten->rcrds[nsts-1][0])
		{
		    for (indx = 1; indx < nsts; ++indx)
	    	        if (Coords(pts[i+1])[0] <= sten->rcrds[indx][0])
			    break;
		    *indxr = indx;
		    *stmpr = sstorer;
		    set_state(sstorer,state_type(sten->sr[0]),
		              left_state(pts[i+1]));
		}
	    }
	    break;
	case 2:
	    interface_crosses_stencil2d(sten,&crxl,&crxr,sstorel,sstorer);
	    break;
	case 3:
	    interface_crosses_stencil3d(sten,&crxl,&crxr,sstorel,sstorer);
	    break;
	default:
	    screen("ERROR in interface_crosses_stencil(), invalid dim %d\n",
	          dim);
	    clean_up(ERROR);
	}
	if ((dim == 2) || (dim == 3))
	{
	    if (crxl != NULL)
	    {
	        *stmpl = sstorel;
		*indxl = (int) ceil(
		    distance_between_positions(crxl,sten->coords,dim)/sten->dn);
	    }
	    if (crxr != NULL)
	    {
	        *stmpr = sstorer;
		*indxr = (int) ceil(
		    distance_between_positions(crxr,sten->coords,dim)/sten->dn);
	    }
	}
}		/*end interface_crosses_stencil*/

LOCAL	void	interface_crosses_stencil2d(
	WSSten  *sten,
	float    **crxl,
	float    **crxr,
	Locstate stl,
	Locstate str)
{
	COMPONENT    pcomp = sten->pcomp;
	COMPONENT    ncomp = sten->ncomp;
	Front        *front = sten->front;
	INTERFACE    *intfc = front->interf;
	POINT        *p = sten->p;
	CURVE        *c = Curve_of_hs(sten->hs);
	CURVE        *cl, *cr;
	BOND         *b = Bond_of_hse(sten->hse);
	BOND         *bl, *br;
	RECT_GRID    *tgr = &topological_grid(intfc);
	float        *nor = sten->nor;
	float        dn = sten->dn;
	float        *crds = sten->coords;
	float        *crds0, *crds1;
	float        dx[3];
	float        dstart[3], dend[3];
	float        sp, sp_start, sp_end;
	float        tn[3];
	float        alpha, beta;
	float        alpha_l, alpha_r;
	float        beta_l, beta_r;
	float        corner[3];
	int          i0[3], i1[3], di[3];
	int          i, j, k, l, dim = tgr->dim;
	int          Nchecked;
	int          nsts = sten->nsts;
	struct Table *T = intfc->table;
	static BOND  **checked;
	static int   size_checked;
	static float crxl_store[3], crxr_store[3];
	static SIDE  **corner_side;
	static int   corner_side_size[2];

	*crxl = *crxr = NULL;

	if ((p == NULL) || (c == NULL) || (b == NULL))
	    return;

	crds0 = sten->lcrds[sten->nsts-1];
	crds1 = sten->rcrds[sten->nsts-1];

	if (checked == NULL)
	{
	    size_checked = 20;
	    uni_array(&checked,size_checked,sizeof(BOND*));
	}
	Nchecked = 0;
	if (p == b->start)
	{
	    if (b->prev)
	        checked[Nchecked++] = b->prev;
	    checked[Nchecked++] = b;
	}
	else if (p == b->end)
	{
	    checked[Nchecked++] = b;
	    if (b->next)
	        checked[Nchecked++] = b->next;
	}
	for (i = 0; i < dim; ++i)
	{
	    int imax = tgr->gmax[i]-1;
	    int imin = 0;

	    j = cell_index(crds0[i],i,tgr);
	    k = cell_index(crds1[i],i,tgr);
	    i0[i] = min(j,k);
	    i0[i] = max(i0[i],imin);
	    i0[i] = min(i0[i],imax);

	    i1[i] = max(j,k);
	    i1[i] = max(i1[i],imin);
	    i1[i] = min(i1[i],imax);

	    di[i] = i1[i] - i0[i] + 1;
	}

	if ((corner_side_size[0] <= di[0]) || (corner_side_size[1] <= di[1]))
	{
	    if (corner_side != NULL)
	        free(corner_side);

	    corner_side_size[0] = 2*di[0] + 1;
	    corner_side_size[1] = 2*di[1] + 1;
	    bi_array(&corner_side,corner_side_size[0],corner_side_size[1],
	           sizeof(SIDE));
	}
	tn[0] =  nor[1];
	tn[1] = -nor[0];
	for (i = 0; i <= di[0]; ++i)
	{
	    corner[0] = cell_edge(i0[0]+i,0,tgr); 
	    dx[0] = corner[0] - crds[0];
	    for (j = 0; j <= di[1]; ++j)
	    {
	        corner[1] = cell_edge(i0[1]+j,1,tgr); 
	        dx[1] = corner[1] - crds[1];
		sp = dx[0]*tn[0] + dx[1]*tn[1];
		if (sp > 0)
		    corner_side[i][j] = POSITIVE_SIDE;
		else if (sp < 0)
		    corner_side[i][j] = NEGATIVE_SIDE;
		else
		    corner_side[i][j] = UNKNOWN_SIDE;
	    }
	}
	alpha_l = -HUGE_VAL;
	alpha_r =  HUGE_VAL;
	beta_l = -HUGE_VAL;
	beta_r = -HUGE_VAL;
	bl = br = NULL;
	cl = cr = NULL;
	for (i = 0; i < di[0]; ++i)
	{
	  for (j = 0; j < di[1]; ++j)
	  {
	    if (
	        (corner_side[i][j] != corner_side[i+1][j]) ||
	        (corner_side[i][j] != corner_side[i][j+1]) ||
	        (corner_side[i][j] != corner_side[i+1][j+1])
	    )
	    { /* Block crosses normal line,  check for crosses */
	      BOND  **b;
	      CURVE **c;
	      int   N;
	      N = T->num_of_bonds[i0[1]+j][i0[0]+i];
	      b = T->bonds[i0[1]+j][i0[0]+i];
	      c = T->curves[i0[1]+j][i0[0]+i];
	      for (k = 0; k < N; ++k, ++b, ++c)
	      {
		if ((!is_subdomain_boundary(Hyper_surf(*c))) &&
		    (!is_passive_boundary(*c)))
		{
		  for (l = 0; l < Nchecked; ++l)
		  {
		    if (*b == checked[l])
		      break;
		  }
		  if (l == Nchecked) /* New bond */
		  {
		    dstart[0] = Coords((*b)->start)[0] - crds[0];
		    dstart[1] = Coords((*b)->start)[1] - crds[1];
		    dend[0] = Coords((*b)->end)[0] - crds[0];
		    dend[1] = Coords((*b)->end)[1] - crds[1];
		    sp_start = dstart[0]*tn[0] + dstart[1]*tn[1];
		    sp_end = dend[0]*tn[0] + dend[1]*tn[1];
		    if (sp_start*sp_end < 0.0) /* Bond Crosses normal line */
		    {
		      beta = sp_start/(sp_start - sp_end);
		      sp_start = dstart[0]*nor[0] + dstart[1]*nor[1];
		      sp_end = dend[0]*nor[0] + dend[1]*nor[1];
		      alpha = sp_start + beta*(sp_end - sp_start);
		      if (fabs(alpha) <= dn*(nsts-1))
		      {
			if ((alpha <= 0.0) && (alpha_l < alpha))
			{
			  alpha_l = alpha;
			  beta_l = beta;
	                  *crxl = crxl_store;
			  crxl_store[0] = crds[0] + alpha*nor[0];
			  crxl_store[1] = crds[1] + alpha*nor[1];
			  bl = *b;
			  cl = *c;
			}
			if ((alpha >= 0.0) && (alpha < alpha_r))
			{
			  alpha_r = alpha;
			  beta_r = beta;
			  *crxr = crxr_store;
			  crxr_store[0] = crds[0] + alpha*nor[0];
			  crxr_store[1] = crds[1] + alpha*nor[1];
			  br = *b;
			  cr = *c;
			}
		      }
		    }
		    if (Nchecked >= size_checked)
	            {
		      BOND **btmp;
	              uni_array(&btmp,2*size_checked,sizeof(BOND*));
		      for (l = 0; l < Nchecked; ++l)
		        btmp[l] = checked[l];
		      free(checked);
		      size_checked *= 2;
		      checked = btmp;
		    }
		    checked[Nchecked++] = *b;
		  }
	        }
	      }
	    }
	  }
	}
	if (is_obstacle_state(sten->sl[0]))
	    *crxl = NULL;
	else if (bl != NULL)
	{
	    if (negative_component(cl) == ncomp)
	    {
	        left_state_along_bond(beta_l,bl,cl,stl);
	    }
	    else if (positive_component(cl) == ncomp)
	    {
	        right_state_along_bond(beta_l,bl,cl,stl);
	    }
	    else
	    {
		float coords[MAXD];
		int i;
		for (i = 0; i < dim; ++i)
		    coords[i] = Coords(bl->start)[0] + beta_l*
			    (Coords(bl->end)[0] - Coords(bl->start)[0]);
		nearest_intfc_state(coords,ncomp,intfc,stl,NULL,NULL);
	        screen("WARNING in interface_crosses_stencil2d(), "
		       "inconsistent negative component\n"
		       "Calling nearest_intfc_state() to resolve\n");
		if (debugging("intfc_crx_stencil"))
		{
		    print_general_vector("Input p = ",Coords(p),dim,"\n");
		    print_general_vector("Input nor = ",nor,dim,"\n");
		    (void) printf("Input curve %llu\n",curve_number(c));
		    (void) printf("Input Bond b\n");
		    print_bond(b);
		    (void) printf("alpha_l = %g\n",alpha_l);
		    (void) printf("beta_l = %g\n",beta_l);
		    (void) printf("Curve cl %llu\n",curve_number(cl));
		    (void) printf("Bond bl\n");
		    print_bond(bl);
		    printf("Propagating curve:\n");
		    print_curve(c);
		    printf("Crossing curve:\n");
		    print_curve(cl);
		}
	    }
	}
	if (is_obstacle_state(sten->sr[0]))
	    *crxr = NULL;
	else if (br != NULL)
	{
	    if (negative_component(cr) == pcomp)
	    {
	        left_state_along_bond(beta_r,br,cr,str);
	    }
	    else if (positive_component(cr) == pcomp)
	    {
	        right_state_along_bond(beta_r,br,cr,str);
	    }
	    else
	    {
		float coords[MAXD];
		int i;
		for (i = 0; i < dim; ++i)
		    coords[i] = Coords(br->start)[0] + beta_r*
			    (Coords(br->end)[0] - Coords(br->start)[0]);
		nearest_intfc_state(coords,pcomp,intfc,str,NULL,NULL);
	        screen("WARNING in interface_crosses_stencil2d(), "
		       "inconsistent positive component\n"
		       "Calling nearest_intfc_state() to resolve\n");
		if (debugging("intfc_crx_stencil"))
		{
		    print_general_vector("Input p = ",Coords(p),dim,"\n");
		    print_general_vector("Input nor = ",nor,dim,"\n");
		    (void) printf("Input curve %llu\n",curve_number(c));
		    (void) printf("Input Bond b\n");
		    print_bond(b);
		    (void) printf("alpha_r = %g\n",alpha_r);
		    (void) printf("beta_r = %g\n",beta_r);
		    (void) printf("Curve cr %llu\n",curve_number(cr));
		    (void) printf("Bond br\n");
		    print_bond(br);
		    printf("Propagating curve:\n");
		    print_curve(c);
		    printf("Crossing curve:\n");
		    print_curve(cr);
		}
	    }
	}
}		/*end interface_crosses_stencil2d*/

EXPORT void g_curve_propagate2d(
	Front           *fr,
        POINTER         wave,
        CURVE           *oldc,
        CURVE           *newc,
        float           dt)
{
	if (wave_type(oldc) == MOVABLE_NEUMANN_BOUNDARY)
	{
	    /* Compute velocity of center of mass and angular velocity */
	    compute_vcm_and_omega(fr,wave,oldc,dt);
	}
	f_curve_propagate2d(fr,wave,oldc,newc,dt);
	if (wave_type(oldc) == MOVABLE_NEUMANN_BOUNDARY)
	{
	    /* Propagate center of mass */
	    int i,dim;
	    dim = fr->rect_grid->dim;
	    for (i = 0; i < dim; ++i)
	    	center_of_mass(newc)[i] = center_of_mass(oldc)[i] +
	    		dt*center_of_mass_velo(oldc)[i];
	}
	return;
}	/* end g_curve_propagate2d */

LOCAL   bool is_bdry_node(
        POINT *p,
        CURVE *c)
{
        if (is_closed_curve(c)) return NO;
        if ((p == c->start->posn && Boundary_hsb(Hyper_surf_bdry(c->start))) ||
            (p == c->end->posn && Boundary_hsb(Hyper_surf_bdry(c->end))))
            return YES;
        return NO;
}       /* end is_bdry_node */

/* 	Function computes velocity of center of mass and
 *  	angular velocity of a regid body, must be a closed curve. 
*/
LOCAL	void compute_vcm_and_omega(
	Front *fr,
	POINTER wave,
	CURVE *curve,
	float dt)
{
	float f[MAXD],total_force[MAXD],rr[MAXD];
	float torque,t,pres;
	BOND *b;
	bool is_right_rigid;
	float nor_in[MAXD];
	int i,dim = fr->rect_grid->dim;
	Locstate ss,se;

	if (is_obstacle_state(right_start_state(curve)))
	    is_right_rigid = YES;
	else 
	    is_right_rigid = NO;
	if (!is_closed_curve(curve))
	{
	    screen("ERROR: in compute_vcm_and_omega(), insufficient ");
	    screen("information to calculate force and torque\n");
	    clean_up(ERROR);
	}
	for (i = 0; i < dim; ++i)
	{
	    total_force[i] = 0.0;
	}
	torque = 0.0;
	for (b = curve->first; b != NULL; b = b->next)
	{
	    if (is_right_rigid)
	    {
		nor_in[0] = Coords(b->end)[1] - Coords(b->start)[1];
		nor_in[1] = Coords(b->start)[0] - Coords(b->end)[0];
	    	if (b == curve->first) 
	    	    ss = left_start_state(curve);
	    	else
	    	    ss = left_state(b->start);
	    	if (b == curve->last) 
	    	    se = left_end_state(curve);
	    	else
	    	    se = left_state(b->end);
	    }
	    else
	    {
		nor_in[1] = Coords(b->end)[0] - Coords(b->start)[0];
		nor_in[0] = Coords(b->start)[1] - Coords(b->end)[1];
	    	if (b == curve->first) 
	    	    ss = right_start_state(curve);
	    	else
	    	    ss = right_state(b->start);
	    	if (b == curve->last) 
	    	    se = right_end_state(curve);
	    	else
	    	    se = right_state(b->end);
	    }
	    pres = 0.5*(pressure(ss) + pressure(se));
	    for (i = 0; i < dim; ++i)
	    {
		f[i] = pres*nor_in[i];
	    	rr[i] = 0.5*(Coords(b->start)[i] + Coords(b->end)[i]);
				- center_of_mass(curve)[i];
	    	total_force[i] += f[i];
	    }
	    Cross2d(f,rr,t);
	    torque += t;
	}
	angular_velo(curve) += dt*torque/mom_inertial(curve);
	for (i = 0; i < dim; ++i)
	{
	    center_of_mass_velo(curve)[i] += 
	    		dt*total_force[i]/total_mass(curve);
	}
	if (debugging("rigid_body"))
	{
	    printf("In compute_vcm_and_omega()\n");
	    printf("total_force = %f %f\n",total_force[0],total_force[1]);
	    printf("torque = %f\n",torque);
	    printf("mom_inertial = %f\n",mom_inertial(curve));
	    printf("total_mass = %f\n",total_mass(curve));
	    printf("center_of_mass = %f %f\n",center_of_mass(curve)[0],
			center_of_mass(curve)[1]);
	    printf("center_of_mass_velo(curve) = %f %f\n",
		center_of_mass_velo(curve)[0],center_of_mass_velo(curve)[1]);
	    printf("angular_velo(curve) = %f\n",angular_velo(curve));
	}
}	/* end compute_vcm_and_omega */

//	#bjet2 interface cross stencil of 3d
typedef  struct {
	POINT	*p;
	float	val;
}	PSET;

typedef  struct {
	TRI	*t;
}	TSET;

#define	MAX_TSET	1000
#define	MAX_PSET	1000

LOCAL	int interpolate_tri_state(
	TRI		*t,
	HYPER_SURF	*hs,
	float		*p,
	Locstate	leftst,
	Locstate	rightst,
	Front		*fr)
{
	Locstate	sl1, sl2, sl3;
	Locstate	sr1, sr2, sr3;
	POINT		*p1,*p2,*p3;
	float		f[3];
	INTERFACE	*intfc = fr->interf;

	set_weight_for_tri_interpolation(p,t,f,NULL,intfc);

	p1 = Point_of_tri(t)[0];
	p2 = Point_of_tri(t)[1];
	p3 = Point_of_tri(t)[2];
	slsr(p1,Hyper_surf_element(t),hs,&sl1,&sr1);
	slsr(p2,Hyper_surf_element(t),hs,&sl2,&sr2);
	slsr(p3,Hyper_surf_element(t),hs,&sl3,&sr3);

	if ((tri_interpolate_states(fr,f[0],f[1],f[2],Coords(p1),sl1,
		                    Coords(p2),sl2,Coords(p3),sl3,
				    leftst)
					!= FUNCTION_SUCCEEDED)
	    ||
	    (tri_interpolate_states(fr,f[0],f[1],f[2],Coords(p1),sr1,
		                    Coords(p2),sr2,Coords(p3),sr3,
				    rightst)
					!= FUNCTION_SUCCEEDED))
	{
	    printf("ERROR in interpolate_tri_state(), "
		   "tri_interpolate_states() failed\n");
	}
}

//calc the range of blocks of a line p1, p2
bool	set_blocks_of_line(
	int		*imin,
	int		*imax,
	float		*p1,
	float		*p2,
	INTERFACE	*intfc)
{
	int		i, ip1[4], ip2[4];
	RECT_GRID	*gr = &topological_grid(intfc);

	for(i=0; i<3; i++)
	{
	    imin[i] = gr->gmax[i] + gr->ubuf[i] + 10;
	    imax[i] = -gr->lbuf[i] - 10;
	}

	ip1[3] = rect_in_which(p1, ip1, gr);
	if(ip1[3])
	{
	    ft_assign(imin, ip1, 3*INT);
	    ft_assign(imax, ip1, 3*INT);
	}
	
	ip2[3] = rect_in_which(p2, ip2, gr);
	if(ip2[3])
	{
	    for(i=0; i<3; i++)
	    {
		imin[i] = min(imin[i], ip2[i]);
		imax[i] = max(imax[i], ip2[i]);
	    }
	}
	if(ip1[3] || ip2[3])
	    return YES;
	return NO;
}

void	reset_point_tri_flag(
	TSET		*tset,
	int		*tnt,
	PSET		*pset,
	int		*pnp,
	int		*imin,
	int		*imax,
	SURFACE		*surf,
	INTERFACE	*intfc)
{
	SURFACE		**s;
	TRI		**t;
	POINT		**p;
	int		i,k,ix,iy,iz,nt;
	struct Table	*T = intfc->table;
	
	//component3d
	for (iz = imin[2]; iz <= imax[2]; ++iz)
	{
	    for (iy = imin[1]; iy <= imax[1]; ++iy)
	    {
	        for (ix = imin[0]; ix <= imax[0]; ++ix)
		{
	            if (T->compon3d[iz][iy][ix] != ONFRONT)
			continue;
	            t = T->tris[iz][iy][ix];
	            s = T->surfaces[iz][iy][ix];
	            nt = T->num_of_tris[iz][iy][ix];

	            for (k = 0; k < nt; k++)
	            {
			if(s[k] != surf)
			    continue;
			if(pset == NULL)
			{
			    //init all the flags.
			    Tri_order(t[k]) = -1;
			    p = Point_of_tri(t[k]);
			    Point_order(p[0]) = -1;
			    Point_order(p[1]) = -1;
			    Point_order(p[2]) = -1;
			}
			else if(Tri_order(t[k]) == -1)
			{
			    //insert tri and pt to sets
			    tset[*tnt].t = t[k];
			    Tri_order(t[k]) = *tnt;
			    (*tnt)++;
			    if(*tnt >= MAX_TSET)
			    {
				printf("ERROR set_point_tri_flag "
				       "tnt is too large.\n");
				clean_up(ERROR);
			    }

			    p = Point_of_tri(t[k]);
			    for(i=0; i<3; i++)
			    {
				if(Point_order(p[i]) != -1)
				    continue;
				pset[*pnp].p = p[i];
				Point_order(p[i]) = *pnp;
				(*pnp)++;
			    }
			    if(*pnp >= MAX_PSET-3)
			    {
				printf("ERROR set_point_tri_flag "
				       "pnp is too large.\n");
				clean_up(ERROR);
			    }
			}
	            }   //for k : tris on a block.
		}
	    }
	}
}

//void    tecplot_sten_points(const char*, WSSten*);

LOCAL	void	interface_crosses_stencil3d(
	WSSten  *sten,
	float    **crxl,
	float    **crxr,
	Locstate stl,
	Locstate str)
{
	Front		*fr = sten->front;
	INTERFACE	*intfc = fr->interf;
	HYPER_SURF	*hs = sten->hs;
	SURFACE		*surf = Surface_of_hs(hs);
	POINT		**p; 
	TRI		*tri, **ptris, *tri1, *trin, *trip;
	TSET		tset[MAX_TSET];
	PSET		pset[MAX_PSET];
	int		nt, i, j, k, pnp, tnt;
	int		nsts = sten->nsts;
	float		dn, t[2][3], v[3],crx[3],*nor,*coords;
	float		dist,distl,distr,tol,n0,n1,n2;
	static float 	crxl_store[3], crxr_store[3];
	static Locstate	tmpst = NULL;
	static float	**pt = NULL;
	static int	**imin, **imax;

	DEBUG_ENTER(interface_crosses_stencil3d)
	
	*crxl = NULL;
	*crxr = NULL;

	if(pt == NULL)
	{
	    bi_array(&pt, 2*nsts, 3, FLOAT);
	    bi_array(&imin, 2*nsts, 3, INT);
	    bi_array(&imax, 2*nsts, 3, INT);
	}
	if(tmpst == NULL)
	    scalar(&tmpst, fr->sizest);

	if(sten->w_type != CONTACT)
	{
	    DEBUG_LEAVE(interface_crosses_stencil3d)
	    return;
	}

	dn = sten->dn;
	nor = sten->nor;
	//set positions of points.
	coords = sten->coords;
	for(i = -nsts+1; i<nsts; i++)
	{
	    for(j=0; j<3; j++)
		pt[i+nsts-1][j] = coords[j] + i*dn*nor[j];
	}
	//set point and tri flags.
	for(i=0; i<2*nsts-2; i++)
	{
	    if(set_blocks_of_line(imin[i], imax[i], 
				    pt[i], pt[i+1], intfc))
	    {
		//printf("#reset enter.\n");
		//print_int_vector("imin", imin[i], 3, "\n");
	        //print_int_vector("imax", imax[i], 3, "\n");

		reset_point_tri_flag(NULL,NULL, NULL,NULL, 
			    imin[i],imax[i],surf,intfc);
	    }
	}

	pnp = 0;
	tnt = 0;
	for(i=0; i<2*nsts-2; i++)
	    reset_point_tri_flag(tset,&tnt, pset,&pnp, 
			    imin[i],imax[i],surf,intfc);

	//printf("#reset %d %d %d\n", nsts, pnp, tnt);

	//finding all tris around point p and deacting all of them.
	tri = Tri_of_hse(sten->hse);
	nt = set_tri_list_around_point(sten->p,tri,&ptris,intfc);
	for(i=0; i<nt; i++)
	{
	    k = Tri_order(ptris[i]);
	    tset[k].t = NULL;
	    Tri_order(ptris[i]) = -1;
	}

	t[0][0] = -nor[1];
	t[0][1] = nor[0];
	t[0][2] = 0.0;
	Cross3d(nor, t[0], t[1]);
	for(i=0; i<2; i++)
	{
	    //calc the coords of p in the projected plan.
	    for(j=0; j<pnp; j++)
	    {
		if(pset[j].p == NULL)
		    continue;

		//distance from point pset[j].p to plane (coords, t[i])
		difference(Coords(pset[j].p), coords, v, 3);
		pset[j].val = Dot3d(v,t[i]);
	    }
	    for(j=0; j<tnt; j++)
	    {
		if(tset[j].t == NULL)
		    continue;

		//if the 3 points of tset[j].t are on the same side of the
		//plane, deacting the tri.
		p = Point_of_tri(tset[j].t);
		n0 = pset[Point_order(p[0])].val;
		n1 = pset[Point_order(p[1])].val;
		n2 = pset[Point_order(p[2])].val;
		if(min3(n0,n1,n2) > 0.0 || max3(n0,n1,n2) < 0.0)
		{
		     Tri_order(tset[j].t) = -1;
		     tset[j].t = NULL;
		}
	    }
	}

	trin = NULL;
	trip = NULL;
	distl = -HUGE_VAL;
	distr = HUGE_VAL;
	tol = 1.0e-6*dn;
	
	for(j=0; j<tnt; j++)
	{
	    tri1 = tset[j].t;
	    if(tri1 == NULL)
		continue;

	    //printf("#tri1 %d\n", tri1);
	    if(!line_tri_crossing(crx,tri1,pt[0],pt[2*nsts-2],tol))
		continue;

	    /*
	    printf("#crx found. %15.9e  %15.9e\n", 
			    distance_between_positions(crx,coords,3),
			    distance_between_positions(pt[0],pt[2*nsts-2],3));
	    print_tri(tri1, intfc);
	    print_general_vector("crx", crx, 3, "\n");
	    print_general_vector("pt0", pt[0], 3, "\n");
	    print_general_vector("pt1", pt[2*nsts-2], 3, "\n");
	    */

	    difference(crx,coords,v,3);
	    dist = Dot3d(v,nor);
	    if(dist < 0.0 && dist > distl)
	    {
		distl = dist;
		trin = tri1;
		*crxl = crxl_store;
		ft_assign(*crxl, crx, 3*FLOAT);
	    }
	    if(dist > 0.0 && dist < distr)
	    {
		distr = dist;
		trip = tri1;
		*crxr = crxr_store;
		ft_assign(*crxr, crx, 3*FLOAT);
	    }
	}

	if(trin != NULL || trip != NULL)
	{
	    //print_general_vector("crds", coords, 3, "\n");
	    //printf("#crxlr  %15.9e %15.9e %15.9e\n",dn, distl, distr);
	    
	    //tecplot_interface_in_ball("intfc_csten", intfc);
	    //tecplot_sten_points("intfc_sten", sten);
	    //add_to_debug("intfc_crx_sten");
	}
	
	if(trin != NULL)
	{
	    interpolate_tri_state(trin, hs, *crxl, stl,tmpst, fr);
	    //verbose_print_state("crx stl", stl);
	}
	if(trip != NULL)
	{
	    interpolate_tri_state(trip, hs, *crxr, tmpst, str, fr);
	    //verbose_print_state("crx str", str);
	}

	//*crxl = NULL;
	//*crxr = NULL;
	
	DEBUG_LEAVE(interface_crosses_stencil3d)
}


