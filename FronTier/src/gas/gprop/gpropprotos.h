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
*				gpropprotos.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#if !defined(_GPROPPROTOS_H)
#define _GPROPPROTOS_H

	/* gprop EXPORTED Function Declarations */
/*	gprop/gipolar.c */
IMPORT	int	i_polar_2(Locstate,Locstate,Locstate,Locstate,ANGLE_DIRECTION,
			  float*,float,float,float*,float*,float);
IMPORT	int	intersection_of_two_shock_polars(Locstate,Locstate,float*,
						 float*,float*,float*,
						 bool,bool,
						 RIEMANN_SOLVER_WAVE_TYPE*,
						 RIEMANN_SOLVER_WAVE_TYPE*);

/*	gprop/gpolar.c */
IMPORT	float	max_behind_shock_pr(float,Locstate);
IMPORT	int	prandtl_meyer_wave(Locstate,float,bool,float*,Locstate,
				   float*,float*,float*);
IMPORT	int	s_polar_2(Locstate,int,int,float,float*,Locstate,float*);
IMPORT	int	s_polar_3(Locstate,int,float,int,int,float*,Locstate,
			  float*,float*);
IMPORT	int	s_polar_4(int,float,float*,float*,Locstate,Locstate,int);
IMPORT	int	s_polar_5(Locstate,int,float,int,int,float*,Locstate,
			  float*,float*);

/*	gprop/gprop.c */
IMPORT	USWSSten2d	*AllocUSWSSten2d(Front*,int,int);
IMPORT	WSSten	*g_AllocWSSten(int,int,Front*);
IMPORT	bool	passive_point_propagate(int,POINT*,POINT*,float*,INTERFACE*);
IMPORT	float	set_pjump_at_wave(POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,
				  Front*,float*);
IMPORT	void	g_pt_prop_by_w_speed(Front*,POINT*,POINT*,
				     HYPER_SURF_ELEMENT*,HYPER_SURF*,
				     float,float*);
IMPORT	void	include_source(float*,Locstate,float,float,float*,float*,float*,
			       RECT_GRID*,float,int);
IMPORT  void    include_source_1(float*,Locstate,float,float,float*,float*,float*,
                               RECT_GRID*,float*,int);
IMPORT	void	print_point_propagate_data(POINT*,HYPER_SURF_ELEMENT*,
					   HYPER_SURF*,int);
IMPORT	void	reflect_wssten(WSSten*,SIDE,Front*);
IMPORT	void	set_point_propagate(Front*,bool);
IMPORT	void	states_near_location(WSSten*,float*,float*,COMPONENT,COMPONENT,
				     Locstate,Locstate);
//#bjet2
IMPORT	void	g_point_propagate_along_wall(Front*,POINTER,POINT*,BOND*,CURVE*,
                                  HYPER_SURF_ELEMENT*,HYPER_SURF*,POINT*,float,float*);
IMPORT	void	pseudo_g_point_propagate_along_wall(Front*,POINTER,POINT*,BOND*,CURVE*,
                                  HYPER_SURF_ELEMENT*,HYPER_SURF*,POINT*,float,float*);


#if defined(TWOD)
IMPORT	void	g_unsplit_w_speed2d(USWSSten2d*,Locstate,Locstate,float*);
IMPORT	void	set_USWSSten2d_off_front_states(USWSSten2d*,HYPER_SURF*,
						Front*,Wave*);
IMPORT	void	set_uswssten_geometry_and_center(USWSSten2d*,POINT*,
						 HYPER_SURF_ELEMENT*,
						 HYPER_SURF*,Front*,float);
IMPORT	void	g_curve_propagate2d(Front*,POINTER,CURVE*,CURVE*,float);
#endif /* defined(TWOD) */

/*	gprop/griemann.c */
IMPORT	bool	oned_state_in_rarefaction_fan(float,float,Locstate,Locstate,
					      Locstate,int,float*,int);
IMPORT	bool	find_mid_state(Locstate,Locstate,float,float*,float*,float*,
			       float*,float*,float*,RIEMANN_SOLVER_WAVE_TYPE*,
			       RIEMANN_SOLVER_WAVE_TYPE*);
IMPORT	int	onedrsoln(float,Locstate,Locstate,Locstate,float*,int);
IMPORT	void	riemann_solution(float,const float*,Locstate,
                                 Locstate,Locstate,int);
#if defined(DEBUG_GRIEMANN)
IMPORT	void	set_debug_find_mid_state(bool);
IMPORT	void	set_debug_riem_sol(bool);
#endif /* defined(DEBUG_GRIEMANN) */
#if defined(COMBUSTION_CODE)

/*	griemcombst.c */
IMPORT	float	CJ_det(Locstate,int,Locstate,int);
IMPORT	bool	combust_find_mid_state(Locstate,Locstate,float*,float*,float*,
				       float*,float*,float*,
				       RIEMANN_SOLVER_WAVE_TYPE*,
				       RIEMANN_SOLVER_WAVE_TYPE*);
IMPORT	int	combust_onedrsoln(float,Locstate,Locstate,Locstate,float*,int);
#endif /* defined(COMBUSTION_CODE) */

/*	gprop/gwspeed.c */
#if defined(DEBUG_W_SPEED)
IMPORT	void	g_set_debug_w_speed(bool);
#endif /* defined(DEBUG_W_SPEED) */

IMPORT	void	g_npt_w_speed(WSSten*,Locstate,Locstate,float*);
IMPORT	void	g_w_speed(float*,Locstate,Locstate,Locstate,Locstate,
			  float*,float,float*,int,Front*);
IMPORT	void	print_g_npt_w_speed_opts(NptWSpeedOpts*);
IMPORT	void	set_npt_wspeed_options(NptWSpeedOpts*);
IMPORT	void	state_behind_sound_wave(Locstate,Locstate,float*,float*,float,
					float,float,float,int,int,
					RIEMANN_SOLVER_WAVE_TYPE,int);

#endif /* !defined(_GPROPPROTOS_H) */
