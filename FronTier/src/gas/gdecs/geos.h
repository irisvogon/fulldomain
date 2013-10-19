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
*				geos.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#if !defined(_GEOS_H)
#define _GEOS_H

#include <gdecs/gstate.h>

struct _EOS {
	/* PRIMARY THERMODYNAMIC FUNCTIONS */
	float	(*_internal_energy)(Locstate);
	float	(*_pressure)(Locstate);
	float   (*_density)(Locstate);
	float	(*_sound_speed_squared)(Locstate);
	float	(*_acoustic_impedance_squared)(Locstate);
	float	(*_specific_internal_energy)(Locstate);
        float   (*_MGamma)(Locstate);

	/* SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS */
	float	(*_specific_enthalpy)(Locstate);
	float	(*_temperature)(Locstate);
	float	(*_entropy)(Locstate);
	float	(*_adiabatic_gamma)(Locstate);
	float	(*_gruneisen_gamma)(Locstate);
	float	(*_fundamental_derivative)(Locstate);
	float	(*_C_V)(Locstate);
	float	(*_C_P)(Locstate);
	float	(*_K_S)(Locstate);
	float	(*_K_T)(Locstate);
        float   (*_specific_enthalpy_species)(Locstate,int);
        float   (*_dynamic_thermal_conductivity)(Locstate,float);
        float   (*_dynamic_viscosity)(Locstate,float);
	void	(*_dynamic_viscosity_thermalconduct)(Locstate,float,float*,float*);
	/* MATERIAL PROPERTY FUNCTIONS */
	float   (*_bulk_viscosity)(Locstate);
	float   (*_shear_viscosity)(Locstate);
	float   (*_heat_coeff)(Locstate);

        /* PHASE TRANSITION FUNCTIONS */
        /* Jun 7 2004: Zhiliang : for phase boundary */
        /** IN pressure and Temperature plane **/
        /** On the saturated vapor and liquid curve **/
        float   (*_saturated_pres)(float);
        float   (*_saturated_temperature)(float);
        /* Jun 10 2004: Zhiliang : for dynamics of phase boundary */
        float   (*_vaporization_heat)(Locstate,float);
        float   (*_thermal_conduct)(Locstate);
        float   (*_density_via_PT)(Locstate,float,float);

        /* Jun 16 2004: Myoung-Nyoun: heat of vaporization for n-heptane */
        float   (*_heat_of_vaporization_via_T)(float);
        float   (*_heat_of_vaporization_via_P)(float);
        float   (*_temperature_via_rhoS)(Locstate,float,float);

	/* VECTORIZED THERMODYNAMIC FUNCTIONS */
	void	(*_single_eos_load_pressure_and_sound_speed2)(struct _Vec_Gas*,
							      int,int);
	void	(*_single_eos_load_pressure_and_gammas)(struct _Vec_Gas*,int,
							int);
	void	(*_single_eos_load_pressure)(struct _Vec_Gas*,int,int);
	void	(*_single_eos_load_sound_speed2)(struct _Vec_Gas*,int,int);

	/* RIEMANN SOLUTIONS UTILITY FUNCTIONS */
	/* Purely Thermodynamic Hugoniot Functions */
	float	(*_dens_Hugoniot)(float,Locstate);
	void	(*_state_w_pr_on_Hugoniot)(Locstate,float,Locstate,int);
	bool	(*_state_w_mf_sqr_on_Hugoniot)(Locstate,float,Locstate,int);

	/* Velocity Related Hugoniot Functions */
	float	(*_pr_normal_vel_wave_curve)(float,Locstate);

	/* Purely Thermodynamic Adiabatic Wave Curve Functions */
	float	(*_dens_rarefaction)(float,Locstate);
	float	(*_pressure_rarefaction)(float,Locstate);
	void	(*_state_on_adiabat_with_pr)(Locstate,float,Locstate,int);
	void	(*_state_on_adiabat_with_dens)(Locstate,float,Locstate,int);

	/* General Wave Curve Functions */
	float	(*_mass_flux)(float,Locstate);
	float	(*_mass_flux_squared)(float,Locstate);

	/* Functions for the Evaluation of Riemann Solutions */
	float	(*_oned_fan_state)(float,Locstate,Locstate,Locstate,int,
				   bool*);

	/* Functions to Compute Riemann Solutions */
	float	(*_riemann_wave_curve)(Locstate,float);
	void	(*_set_state_for_find_mid_state)(Locstate,Locstate);
	float	(*_eps_for_Godunov)(Locstate,float,float);
	void	(*_initialize_riemann_solver)(Locstate,Locstate,float*,float*,
					      float,float*,float*,
					      bool(*)(Locstate,Locstate,
						      float,float*,float*,
						      float*,float*,float*,
						      float*,
						      RIEMANN_SOLVER_WAVE_TYPE*,
						      RIEMANN_SOLVER_WAVE_TYPE*));


	/* TWO DIMENSIONAL RIEMANN SOLUTION UTILTITY FUNCTIONS */
	bool	(*_steady_state_wave_curve)(float,float,float*,Locstate);
	float	(*_pressure_at_sonic_point)(float,Locstate);
	bool	(*_pr_at_max_turn_angle)(float*,float,Locstate);
	float	(*_state_in_prandtl_meyer_wave)(float,float,Locstate,float,
						Locstate,Locstate,int);
#if defined(COMBUSTION_CODE)
	/* DETONATION SPECIFIC UTILITY FUNCTIONS */
	float	(*_CJ_state)(Locstate,int,Locstate,int,int);
	void	(*_progress_state)(float,Locstate,Locstate,float);
	void	(*_fprint_combustion_params)(FILE*,struct _Gas_param*);
	/* flame */
        float   (*_state_w_pr_and_heat_on_Hugoniot)(float,Locstate,Locstate,int);
        float   (*_CJ_deflagration_state)(Locstate,int,Locstate,int,int);
        float   (*_pre_deflagration_density)(float,Locstate);
        float   (*_post_deflagration_pressure)(float,Locstate,float,float);
        float   (*_post_deflagration_vel)(float,Locstate,float,float,float,int);
	/* flame */
#endif /* defined(COMBUSTION_CODE) */

	/* METHOD OF CHARACTERISTIC FUNCTIONS FOR W_SPEED */
	void	(*_neumann_riem_inv_moc)(float*,Locstate,float,float,Locstate,
					 SIDE,Locstate,float,float*,Front*);
	void	(*_shock_ahead_state_riem_inv_moc)(float*,Locstate,Locstate,
						   Locstate,Locstate,Locstate,
						   float,float,float,float,
						   float*,float*,int,float,
						   Front*);
	bool	(*_shock_moc_plus_rh)(float*,Locstate,Locstate,Locstate,
				      Locstate,float,float*,float*,int,Front*);

	/* INITIALIZATION UTILITY FUNCTIONS */
	void	(*_prompt_for_state)(Locstate,int,struct _Gas_param*,
				     const char*);
	void	(*_prompt_for_thermodynamics)(Locstate,struct _Gas_param*,
					      const char*);
	void	(*_fprint_EOS_params)(FILE*,struct _Gas_param*);
	void	(*_read_print_EOS_params)(INIT_DATA*,const IO_TYPE*,
	                                  struct _Gas_param*);
	struct _EOS*	(*_free_EOS_params)(struct _EOS*);
	void	(*_prompt_for_EOS_params)(INIT_DATA*,struct _Gas_param*,
					  const char*,const char*);
	/* Problem Type Specific Initialization Functions */
	float	(*_RT_RS_f)(float,Locstate,float,float,float);
	void	(*_RT_single_mode_perturbation_state)(Locstate,float*,float,
						      Locstate,float,float,
						      struct _MODE*,float);
	void	(*_KH_single_mode_state)(Locstate,float*,float,
					 Locstate,float,float,float,
					 struct _MODE*);
	void	(*_compute_isothermal_stratified_state)(Locstate,float,float,
							Locstate);
	void	(*_compute_isentropic_stratified_state)(Locstate,float,float,
							Locstate);
	void	(*_compute_constant_density_stratified_state)(Locstate,float,
							      float,Locstate);
	float	(*_pressure_from_entropy)(Locstate,float);

	bool	_multiphase_eos;
        bool    _compute_ns_terms;

	/* Equation of state domain functions */
	float   (*_Min_energy)(Locstate);
	float   (*_Min_pressure)(Locstate);
	float   (*_Vacuum_dens)(Locstate);
	float   (*_Raref_press)(Locstate);
#if defined(COMBUSTION_CODE)
	float   (*_Tol_alpha)(Locstate);
	float   (*_Tol_pressure)(Locstate);
#endif /* defined(COMBUSTION_CODE) */
};
typedef struct _EOS EOS;

#define Eos(state)			Params(state)->eos
#define	acoustic_impedance(state)	sqrt(acoustic_impedance_squared(state))
#define	sound_speed(state)		sqrt(sound_speed_squared(state))
#define	multiphase_eos(eos)		((EOS *)eos)->_multiphase_eos

/**** PHASE TRANSITION FUNCTION ****/
/* Jun 8 2004: Myoung-Nyoun : for phase boundary */
#define saturated_pres(state,temperature)                               \
        (*Eos(state)->_saturated_pres)(temperature)
#define saturated_temperature(state,pres)                               \
        (*Eos(state)->_saturated_temperature)(pres)
/* Jun 10 2004: Zhiliang : for dynamics of phase boundary */
#define vaporization_heat(state,pres)                                   \
        (*Eos(state)->_vaporization_heat)((state),(pres))
#define thermal_conduct(state)                                          \
        (*Eos(state)->_thermal_conduct)((state))
#define density_via_PT(state,pres,temp)                                 \
        (*Eos(state)->_density_via_PT)((state),(pres),(temp))
#define temperature_via_rhoS(state,entropy,rho)                         \
        (*Eos(state)->_temperature_via_rhoS)(state,(entropy),(rho))

/* Jun 16 2004: Myoung-Nyoun: heat of vaporization for n-heptane */
#define heat_of_vaporization_via_T(state,temp)                          \
        (*Eos(state)->_heat_of_vaporization_via_T)((temp))
#define heat_of_vaporization_via_P(state,press)                         \
        (*Eos(state)->_heat_of_vaporization_via_P)((press))

	/* PRIMARY THERMODYNAMIC FUNCTIONS */
#define	internal_energy(state)						\
	(*Eos(state)->_internal_energy)(state)
#define	pressure(state)							\
	(*Eos(state)->_pressure)(state)
#define	density(state)							\
	(*Eos(state)->_density)(state)
#define	sound_speed_squared(state)					\
	(*Eos(state)->_sound_speed_squared)(state)
#define	acoustic_impedance_squared(state)				\
	(*Eos(state)->_acoustic_impedance_squared)(state)
#define	specific_internal_energy(state)					\
	(*Eos(state)->_specific_internal_energy)(state)
#define MGamma(state)                                                   \
        (*Eos(state)->_MGamma)(state)

	/* SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS */
#define	specific_enthalpy(state)					\
	(*Eos(state)->_specific_enthalpy)(state)
#define	temperature(state)						\
	(*Eos(state)->_temperature)(state)
#define	entropy(state)							\
	(*Eos(state)->_entropy)(state)
#define	adiabatic_gamma(state)						\
	(*Eos(state)->_adiabatic_gamma)(state)
#define	gruneisen_gamma(state)						\
	(*Eos(state)->_gruneisen_gamma)(state)
#define	fundamental_derivative(state)					\
	(*Eos(state)->_fundamental_derivative)(state)
#define	C_V(state)							\
	(*Eos(state)->_C_V)(state)
#define	C_P(state)							\
	(*Eos(state)->_C_P)(state)
#define	K_S(state)							\
	(*Eos(state)->_K_S)(state)
#define	K_T(state)							\
	(*Eos(state)->_K_T)(state)
#define specific_enthalpy_species(state,i)                              \
        (*Eos(state)->_specific_enthalpy_species)(state,i)             
#define dynamic_viscosity(state,T_ref)                              	\
        (*Eos(state)->_dynamic_viscosity)(state,T_ref)             
#define dynamic_thermal_conductivity(state,T_ref)                       \
        (*Eos(state)->_dynamic_thermal_conductivity)(state,T_ref)     
#define dynamic_viscosity_thermalconduct(state,T_ref,mu,kappa)         \
        (*Eos(state)->_dynamic_viscosity_thermalconduct)(state,T_ref,mu,kappa)

	/* MATERIAL PROPERTY FUNCTIONS */
#define	bulk_viscosity(state)						\
	(*Eos(state)->_bulk_viscosity)(state)
#define	shear_viscosity(state)						\
	(*Eos(state)->_shear_viscosity)(state)
#define heat_coeff(state)                                               \
        (*Eos(state)->_heat_coeff)(state)

	/* RIEMANN SOLUTIONS UTILITY FUNCTIONS */
	/* Purely Thermodynamic Hugoniot Functions */
#define	dens_Hugoniot(p1,state0)					\
	(*Eos(state0)->_dens_Hugoniot)(p1,state0)
#define	state_w_pr_on_Hugoniot(state0,p1,state1,stype1)			\
	(*Eos(state0)->_state_w_pr_on_Hugoniot)(state0,p1,state1,stype1)
#define	state_w_mf_sqr_on_Hugoniot(state0,m2,state1,stype1)		\
	(*Eos(state0)->_state_w_mf_sqr_on_Hugoniot)(state0,m2,state1,stype1)

	/* Velocity Related Hugoniot Functions */
#define	pr_normal_vel_wave_curve(du,state0)				\
	(*Eos(state0)->_pr_normal_vel_wave_curve)(du,state0)

	/* Purely Thermodynamic Adiabatic Wave Curve Functions */
#define	dens_rarefaction(p1,state0)					\
	(*Eos(state0)->_dens_rarefaction)(p1,state0)
#define	pressure_rarefaction(rho1,state0)				\
	(*Eos(state0)->_pressure_rarefaction)(rho1,state0)
#define	state_on_adiabat_with_pr(state0,p1,state1,stype1)		\
	(*Eos(state0)->_state_on_adiabat_with_pr)(state0,p1,state1,stype1)
#define	state_on_adiabat_with_dens(state0,rho1,state1,stype1)		\
	(*Eos(state0)->_state_on_adiabat_with_dens)(state0,rho1,state1,stype1)

	/* General Wave Curve Functions */
#define	mass_flux(p,state0)						\
	(*Eos(state0)->_mass_flux)(p,state0)
#define	mass_flux_squared(p,state0)					\
	(*Eos(state0)->_mass_flux_squared)(p,state0)

	/* Functions for the Evaluation of Riemann Solutions */
#define	oned_fan_state(w,sta,stb,stm,stype_m,vacuum)			\
	(*Eos(sta)->_oned_fan_state)(w,sta,stb,stm,stype_m,vacuum)

	/* Functions to Compute Riemann Solutions */
#define	riemann_wave_curve(state0,pstar)				\
	(*Eos(state0)->_riemann_wave_curve)(state0,pstar)
#define	set_state_for_find_mid_state(Tst,st)				\
	(*Eos(st)->_set_state_for_find_mid_state)(Tst,st)
#define	eps_for_Godunov(state,pstar,r_eps)				\
	(*Eos(state)->_eps_for_Godunov)(state,pstar,r_eps)
#define	initialize_riemann_solver(sl,sr,ps,pm,eps,eps_u,eps_p,f)	\
	(*Eos(sl)->_initialize_riemann_solver)(sl,sr,ps,pm,eps,eps_u,eps_p,f)


	/* TWO DIMENSIONAL RIEMANN SOLUTION UTILTITY FUNCTIONS */
#define	steady_state_wave_curve(p1,M0sq,theta,state0)			\
	(*Eos(state0)->_steady_state_wave_curve)(p1,M0sq,theta,state0)
#define	pressure_at_sonic_point(M0sq,state0)				\
	(*Eos(state0)->_pressure_at_sonic_point)(M0sq,state0)
#define	pr_at_max_turn_angle(prm,M0sq,state0)				\
	(*Eos(state0)->_pr_at_max_turn_angle)(prm,M0sq,state0)
#define	state_in_prandtl_meyer_wave(w,A_a,sta,A_b,stb,stm,stype_m)	\
	(*Eos(sta)->_state_in_prandtl_meyer_wave)(w,A_a,sta,A_b,stb,stm,stype_m)

#if defined(COMBUSTION_CODE)
	/* DETONATION SPECIFIC UTILITY FUNCTIONS */
#define	CJ_state(CJ,st_type_CJ,start,l_or_r,avail)			\
	(*Eos(start)->_CJ_state)(CJ,st_type_CJ,start,l_or_r,avail)
#define	progress_state(prog,init,ans,max_vol)				\
	(*Eos(init)->_progress_state)(prog,init,ans,max_vol)
#define	fprint_combustion_params(file,params)				\
	(*(params)->eos->_fprint_combustion_params)(file,params)
#endif /* defined(COMBUSTION_CODE) */

	/* METHOD OF CHARACTERISTIC FUNCTIONS FOR W_SPEED */
#define	neumann_riem_inv_moc(pt,st0,u0,c0,st1,side,ans,dn,nor,fr)	\
	(*Eos(st0)->_neumann_riem_inv_moc)(pt,st0,u0,c0,st1,side,ans,dn,nor,fr)
#define	shock_ahead_state_riem_inv_moc(p,s0,s1,s2,s3,a,dn,f1,f2,f3,n,W,as,dt,f)\
	(*Eos(s0)->_shock_ahead_state_riem_inv_moc)(p,s0,s1,s2,s3,a,dn, \
						    f1,f2,f3,n,W,as,dt,f)
#define	shock_moc_plus_rh(pt,sta,stm,stb,ans,dn,n,W,wt,fr) \
	(*Eos(sta)->_shock_moc_plus_rh)(pt,sta,stm,stb,ans,dn,n,W,wt,fr)

	/* INITIALIZATION UTILITY FUNCTIONS */
#define	prompt_for_state(state,stype,params,mesg)			\
	(*(params)->eos->_prompt_for_state)(state,stype,params,mesg)
#define	prompt_for_thermodynamics(state,params,mesg)			\
	(*(params)->eos->_prompt_for_thermodynamics)(state,params,mesg)
#define	fprint_EOS_params(file,params)					\
	(*(params)->eos->_fprint_EOS_params)(file,params)
#define	read_print_EOS_params(init,io_type,params)			\
	(*(params)->eos->_read_print_EOS_params)(init,io_type,params)
#define	free_EOS_params(eos)						\
	(*(eos)->_free_EOS_params)(eos)
#define	prompt_for_EOS_params(init,params,message1,message2)		\
	(*(params)->eos->_prompt_for_EOS_params)(init,params,message1,message2)
#define	RT_RS_f(s,amb_st,dz,k_sqr,g_z)					\
	(*Eos(amb_st)->_RT_RS_f)(s,amb_st,dz,k_sqr,g_z)
#define	RT_single_mode_perturbation_state(ans,crds,t,st,zi,zb,m,gz)	\
	(*Eos(st)->_RT_single_mode_perturbation_state)(ans,crds,t,st,zi,zb,m,gz)
#define	KH_single_mode_state(ans,crds,t,st,sv,zi,zb,m)			\
	(*Eos(st)->_KH_single_mode_state)(ans,crds,t,st,sv,zi,zb,m)
#define	compute_isothermal_stratified_state(ans,dz,gz,ref_st)		\
	(*Eos(ref_st)->_compute_isothermal_stratified_state)(ans,dz,gz,ref_st)
#define	compute_isentropic_stratified_state(ans,dz,gz,ref_st)		\
	(*Eos(ref_st)->_compute_isentropic_stratified_state)(ans,dz,gz,ref_st)
#define	compute_constant_density_stratified_state(ans,dz,gz,ref_st)	\
	(*Eos(ref_st)->_compute_constant_density_stratified_state)(ans,dz,gz,ref_st)
#define specific_internal_energy_h_species(state)                       \
         (*Eos(state)->_specific_internal_energy_h_species)(state)
#define specific_internal_energy_l_species(state)                       \
         (*Eos(state)->_specific_internal_energy_l_species)(state)

	/* Equation of state domain functions */
#define	Min_energy(state)	(*Eos(state)->_Min_energy)(state)
#define	Min_pressure(state)	(*Eos(state)->_Min_pressure)(state)
#define	Vacuum_dens(state)	(*Eos(state)->_Vacuum_dens)(state)
#define	Raref_press(state)	(*Eos(state)->_Raref_press)(state)
#if defined(COMBUSTION_CODE)
#define	Tol_alpha(state)	(*Eos(state)->_Tol_alpha)(state)
#define	Tol_pressure(state)	(*Eos(state)->_Tol_pressure)(state)
#endif /* defined(COMBUSTION_CODE) */

#define	pressure_from_entropy(state,beta)	(*Eos(state)->_pressure_from_entropy)(state,beta)
		
		/* Equation of state types */

enum {
	UNKNOWN_EOS             = -3,
	OBSTACLE_EOS		=  1,
	FIRST_PHYSICAL_EOS_TYPE	= 10
};

#endif /* !defined(_GEOS_H) */
