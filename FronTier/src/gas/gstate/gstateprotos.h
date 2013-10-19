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
*				gstateprotos.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	This file contains the structures for gas dynamics states.
*/

#if !defined(_GSTATEPROTOS_H)
#define _GSTATEPROTOS_H

	/* gstate EXPORTED Function Declarations*/
	/* gintrp.c*/
IMPORT	CONSISTENT_PARAMS consistent_params_in_lin_comb(Locstate,Locstate);
IMPORT	bool	g_bad_state_data(const char*,Front*,Wave*);
IMPORT	bool	g_invalid_state(const char*,Locstate,bool);
IMPORT	bool	is_bad_state(Locstate,bool,const char*);
IMPORT	float	pos_radius(float,RECT_GRID*);
IMPORT	void	g_lin_comb_states(float,float,float*,Locstate,float*,
				  Locstate,RECT_GRID*,Locstate);
IMPORT	void	g_linear_combination_of_states(const float*,const Locstate*,
					       const int,Locstate);
IMPORT	void	gt_lin_comb_states(float,float,float*,Locstate,float*,
				   Locstate,RECT_GRID*,Locstate);
IMPORT	void    inconsistent_params_in_lin_comb_message(float,float,float*,
						        Locstate,float*,
							Locstate,const char*,
							bool);
IMPORT	void    obstacle_state_found_in_lin_comb_message(float,float,float*,
						         Locstate,float*,
	                                                 Locstate,const char*,
	                                                 bool);
IMPORT	void	set_interpolation_weights(INTERPOLATION_WEIGHT);
#if defined(TWOD) || defined(THREED)
IMPORT	CONSISTENT_PARAMS consistent_params_in_tri_lin_comb(Locstate,
							    Locstate,
							    Locstate);
IMPORT	bool	g_tri_interpolator(float*,LINEAR_ELEMENT*,TRI_SOLN*,Locstate);
IMPORT	bool	g_tri_lin_comb_states(float,float,float,float*,Locstate,float*,
				      Locstate,float*,Locstate,RECT_GRID*,
				      Locstate);
IMPORT	bool	gt_tri_interpolator(float*,LINEAR_ELEMENT*,TRI_SOLN*,Locstate);
IMPORT	bool	gt_tri_lin_comb_states(float,float,float,float*,Locstate,
				       float*,Locstate,float*,
				       Locstate,RECT_GRID*,Locstate);
IMPORT void inconsistent_params_in_tri_lin_comb_mesg(float,float,float,float*,
	                                             Locstate,float*,Locstate,
						     float*,Locstate,
						     const char*,bool);
IMPORT bool g_least_sqr_interpolator(LEAST_SQR_CLUSTER*,TRI_SOLN*,
					float*,Locstate);
#endif /* defined(TWOD) || defined(THREED) */
#if defined(TWOD)
IMPORT	void	g_grad_quad_interpolator(float*,float*,BILINEAR_ELEMENT*,
					 TRI_SOLN*,Locstate*);
IMPORT	bool	g_grad_tri_interpolator(float**,LINEAR_ELEMENT*,TRI_SOLN*,
					Locstate*);
IMPORT	void	g_quad_interpolator(float*,BILINEAR_ELEMENT*,TRI_SOLN*,
				    Locstate);
IMPORT	void	gt_quad_interpolator(float*,BILINEAR_ELEMENT*,
				     TRI_SOLN*,Locstate);
#endif /* defined(TWOD) */
#if defined(THREED)
IMPORT	bool	g_tetra_interpolator(float*,LINEAR_ELEMENT*,
				     TRI_SOLN*,Locstate);
IMPORT	bool	gt_tetra_interpolator(float*,LINEAR_ELEMENT*,
				      TRI_SOLN*,Locstate);
IMPORT	void	g_cube_interpolator(float*,BILINEAR_ELEMENT*,TRI_SOLN*,
				    Locstate);
IMPORT	void	gt_cube_interpolator(float*,BILINEAR_ELEMENT*,
				     TRI_SOLN*,Locstate);
#endif /* defined(THREED) */

	/* gbstate.c*/
IMPORT	void	constant_pressure_flow_through_boundary_state(float*,
				HYPER_SURF*,Front*,POINTER,Locstate);
IMPORT	void	time_dep_pressure_flow_through_boundary_state(float*,
				HYPER_SURF*,Front*,POINTER,Locstate);
IMPORT	void	flow_through_boundary_state(float*,HYPER_SURF*,Front*,
					    POINTER,Locstate);
IMPORT	void	g_fixed_boundary_state(float*,HYPER_SURF*,Front*,
				       POINTER,Locstate);
IMPORT	int	g_neumann_bdry_state(float*,COMPONENT,POINT*,HYPER_SURF*,
				     Front*,POINTER,Locstate);
IMPORT	void    g_time_dependent_boundary_state(float*,HYPER_SURF*,Front*,
                                                POINTER,Locstate);
IMPORT	void 	rk_new_open_bdry( double *crds, Locstate sim_state,
				Locstate bdry_state,
				float *dir, // out of domain direction
				/*out*/ Locstate sten_state);
IMPORT  void    rk_new_open_bdry_a( double *crds, Locstate sim_state,
                                Locstate bdry_state,
                                float *dir, // out of domain direction
                                /*out*/ Locstate sten_state);
	/* grstate.c*/
IMPORT	bool axisymmetric_random_region_about_origin(RANDOM_STATE*,
							INTERFACE*);
IMPORT	int	in_correlation_ellipse(int,int,int,RANDOM_STATE*);
IMPORT	void	generate_random_region(float,Locstate***,RANDOM_STATE*,
				       INTERFACE*);
IMPORT	void	random_velocity_inlet(float*,HYPER_SURF*,Front*,
				      POINTER,Locstate);
IMPORT	void	relax_random_inlet_level(Locstate***,Locstate***,
					 RANDOM_STATE*);

	/* gstate.c*/
IMPORT	Locstate	g_alloc_intfc_state(size_t);
IMPORT	OUTPUT_VALUE	*g_solution(OUTPUT_SOLN*,float*,int*);
IMPORT	float	*MomentumVector(Locstate,float*);
IMPORT	float	*VelocityVector(Locstate,float*);
IMPORT	float	RadialComponentOfVelocity(Locstate,float*,int);
IMPORT	float	energy(Locstate);
IMPORT	float	kinetic_energy(Locstate);
IMPORT	float	mach_number(Locstate,float*);
IMPORT	float	mach_number_squared(Locstate,float*,float*);
IMPORT	float	max_speed(Locstate);
IMPORT	float	mom(int,Locstate);
IMPORT	float	react(Locstate);
IMPORT	float	reaction_rate(Locstate);
IMPORT	float	vel(int,Locstate);
IMPORT	void	add_velocity_to_state(Locstate,float*);
IMPORT	void	g_alloc_state(Locstate*,size_t);
IMPORT	void	g_bundle_states(int*,int*,Wave*,byte*);
IMPORT	void	g_clear_state(Locstate,size_t);
IMPORT	void	g_initialize_max_front_speed(Front*);
IMPORT	void	g_initialize_max_wave_speed(Wave*);
IMPORT	void	g_intfc_solution(OUTPUT_SOLN*,POINT*,HYPER_SURF_ELEMENT*,
				 HYPER_SURF*,OUTPUT_VALUE*,OUTPUT_VALUE*);
IMPORT	void	g_principal_tangent(POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,
				    float*,float*);
IMPORT	void	g_reflect_and_stratify_state(Locstate,INTERFACE*,
					     float*,float*,float*);
IMPORT	void	g_reflect_state(Locstate,INTERFACE*,float*,float*,float*);
IMPORT	void	g_set_state(Locstate,int,Locstate);
IMPORT	void	g_set_max_front_speed(int,float,Locstate,float*,Front*);
IMPORT	void	g_set_max_wave_speed(int,float,Locstate,float*,Wave*);
IMPORT	void	g_transform_state(Locstate,AFLIN*);
IMPORT	void	g_unbundle_states(int*,int*,Wave*,byte*);
IMPORT	void	zero_normal_velocity(Locstate,float*,int);
IMPORT	void	zero_state_velocity(Locstate,int);
IMPORT	void	alpha_state_velocity(float,Locstate,int);
IMPORT	void	g_check_front_state_consistency(Front*);
#if defined(TWOD)
IMPORT	float	TangCmptOfXYVelocity(Locstate,float*,int);
#endif /* defined(TWOD) */
#if defined(THREED)
IMPORT  void    check_for_consistent_tri_states(INTERFACE*);
#endif /* defined(THREED) */
#if defined(COMBUSTION_CODE)
IMPORT float    flame_velocity(Locstate,Locstate);
#endif /* defined(COMBUSTION_CODE) */

	/* gstglobs.c*/
IMPORT	GEOMETRY	Geometry(float*);
IMPORT	GEOMETRY	coord_system(void);
IMPORT	Locstate	return_obst_state(void);
IMPORT	bool	is_gravity(void);
IMPORT	bool	source_terms_exist(void);
IMPORT	const float	*gravity(const float*,const float);
IMPORT	float	rotational_symmetry(void);
IMPORT	int	g_composition_type(void);
IMPORT	int	g_nfloats(void);
IMPORT	size_t	g_sizest(void);
IMPORT	void	eval_gravity(const float*,const float,float*);
IMPORT	void	default_artificial_viscosity(AVISC*);
IMPORT	void	g_set_sizeof_state(CHART*,size_t,int);
IMPORT	void	g_obstacle_state(Locstate,size_t);
IMPORT	void	set_composition_type(int);
IMPORT	void	set_coord_sys(int,int);
IMPORT	void	set_default_artificial_viscosity(AVISC*);
IMPORT	void	set_gravity(GRAVITY*);
IMPORT  void	use_artificial_dissipation(AVISC*);
IMPORT  float	astrophys_grav_constant();

	/* gstate/gparams.c*/
IMPORT	Gas_param	*Params_of_index(int);
IMPORT	Gas_param	*alloc_Gas_param(size_t);
IMPORT	bool	valid_params(Locstate);
IMPORT	int	index_of_Gas_param(Gas_param*);
IMPORT	int	return_params_list(Gas_param***);
IMPORT  void	g_set_params(Locstate,Locstate);
IMPORT	void	set_params_list(INTERFACE*);
IMPORT	uint64_t gas_param_number(Gas_param*);

        /* gstate/goverinerp.c */
#if defined(USE_OVERTURE)
IMPORT  void    g_overture_assign_wave_params(Locstate,Locstate);
IMPORT  void    g_overture_assign_wave_st_type(Locstate,Locstate);
IMPORT  void    g_overture_to_ft_st(Locstate,POINTER,int,int*);
IMPORT  void    g_ft_to_overture_st(Locstate,POINTER,int,int*);

IMPORT  void    g_trans_wv_st_to_overfunc(Wave*,Front*,POINTER,int,int);
IMPORT  void    g_fill_root_extr_overture_cell_st_from_wv(Wave*,Front*,
                                     POINTER,int,int);
IMPORT  bool    g_overture_init_interpolation_coarse_to_fine(Wave**,Front**);
IMPORT  bool    g_overture_interpolation_fine_to_coarse(Wv_on_pc**,
                   Wave**,Front**,Wave***,Front***,int,int);
IMPORT  bool    g_overture_undistribute_interpolation_fine_to_coarse(
                     Wave**,Front**);
IMPORT  bool    g_scatter_patch_states(Overparam*,Wv_on_pc**,Wave**,
                   Front**,int*);
IMPORT  bool    g_scatter_patch_states_in_sweep_dir(Overparam*,Wv_on_pc**,
                  struct _Wave**,Front**,int*,int);
IMPORT  bool    g_overture_fill_patch_amr_buffer2d(Wave**,Front**);
IMPORT  bool    g_overture_fill_patch_amr_buffer_pt2d(int*,int,Wave**,Front**);
IMPORT  bool    g_overture_injection_after_repatch(Wave**,Front**,Wave**,Front**);
#endif /* defined(USE_OVERTURE) */

#endif /* !defined(_GSTATEPROTOS_H) */
