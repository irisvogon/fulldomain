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
*				mpoly-eos.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Template file for the implementation of new equation of state models
*
*/

#if defined(MULTI_COMPONENT)
#include <geos/mpoly.h>

	/* LOCAL Function Prototypes */
	/* PRIMARY THERMODYNAMIC FUNCTIONS */
LOCAL	float	MPOLY_internal_energy(Locstate);
LOCAL	float	MPOLY_pressure(Locstate);
LOCAL	float	MPOLY_sound_speed_squared(Locstate);
LOCAL	float	MPOLY_acoustic_impedance_squared(Locstate);
LOCAL	float	MPOLY_specific_internal_energy(Locstate);

	/* SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS */
LOCAL	float	MPOLY_specific_enthalpy(Locstate);
LOCAL	float	MPOLY_temperature(Locstate);
LOCAL	float	MPOLY_entropy(Locstate);
LOCAL	float	MPOLY_adiabatic_gamma(Locstate);
LOCAL	float	MPOLY_gruneisen_gamma(Locstate);
LOCAL	float	MPOLY_fundamental_derivative(Locstate);
LOCAL	float	MPOLY_C_V(Locstate);
LOCAL	float	MPOLY_C_P(Locstate);
LOCAL	float	MPOLY_K_S(Locstate);
LOCAL	float	MPOLY_K_T(Locstate);
LOCAL   float   MPOLY_specific_enthalpy_species(Locstate,int);
LOCAL   float   MPOLY_dynamic_viscosity(Locstate,float);
LOCAL   float   MPOLY_dynamic_thermal_conductivity(Locstate,float);
LOCAL   void    MPOLY_dynamic_viscosity_thermalconduct(Locstate,float,float*,float*);
LOCAL   float   Collision_intergral_2(float);
LOCAL   float   Collision_intergral_1(float);
LOCAL   float   C_P_species(Locstate,float,int);

	/* VECTORIZED THERMODYNAMIC FUNCTIONS */
LOCAL	void	MPOLY_single_eos_load_pressure_and_sound_speed2(Vec_Gas*,
							        int,int);
LOCAL	void	MPOLY_single_eos_load_pressure_and_gammas(Vec_Gas*,int,int);
LOCAL	void	MPOLY_single_eos_load_pressure(Vec_Gas*,int,int);
LOCAL	void	MPOLY_single_eos_load_sound_speed2(Vec_Gas*,int,int);

	/* RIEMANN SOLUTIONS UTILITY FUNCTIONS */
	/* Purely Thermodynamic Hugoniot Functions */
LOCAL	float	MPOLY_dens_Hugoniot(float,Locstate);
LOCAL	void	MPOLY_state_w_pr_on_Hugoniot(Locstate,float,Locstate,int);
LOCAL	bool	MPOLY_state_w_mf_sqr_on_Hugoniot(Locstate,float,Locstate,int);

	/* Velocity Related Hugoniot Functions */
LOCAL	float	MPOLY_pr_normal_vel_wave_curve(float,Locstate);

	/* Purely Thermodynamic Adiabatic Wave Curve Functions */
LOCAL	float	MPOLY_dens_rarefaction(float,Locstate);
LOCAL	float	MPOLY_pressure_rarefaction(float,Locstate);
LOCAL	void	MPOLY_state_on_adiabat_with_pr(Locstate,float,Locstate,int);
LOCAL	void	MPOLY_state_on_adiabat_with_dens(Locstate,float,Locstate,int);

	/* General Wave Curve Functions */
LOCAL	float	MPOLY_mass_flux(float,Locstate);
LOCAL	float	MPOLY_mass_flux_squared(float,Locstate);

	/* Functions for the Evaluation of Riemann Solutions */
LOCAL	float	MPOLY_oned_fan_state(float,Locstate,Locstate,Locstate,int,
				     bool*);

	/* Functions to Compute Riemann Solutions */
LOCAL	float	MPOLY_riemann_wave_curve(Locstate,float);
LOCAL	void	MPOLY_set_state_for_find_mid_state(Locstate,Locstate);
LOCAL	float	MPOLY_eps_for_Godunov(Locstate,float,float);
LOCAL	void	MPOLY_initialize_riemann_solver(Locstate,Locstate,float*,
						float*,float,float*,float*,
						bool(*)(Locstate,Locstate,
							float,float*,float*,
							float*,float*,float*,
							float*,
							RIEMANN_SOLVER_WAVE_TYPE*,
							RIEMANN_SOLVER_WAVE_TYPE*));


	/* TWO DIMENSIONAL RIEMANN SOLUTION UTILTITY FUNCTIONS */
LOCAL	bool	MPOLY_steady_state_wave_curve(float,float,float*,Locstate);
LOCAL	float	MPOLY_pressure_at_sonic_point(float,Locstate);
LOCAL	bool	MPOLY_pr_at_max_turn_angle(float*,float,Locstate);
LOCAL	float	MPOLY_state_in_prandtl_meyer_wave(float,float,Locstate,float,
						  Locstate,Locstate,int);

#if defined(COMBUSTION_CODE)
	/* DETONATION SPECIFIC UTILITY FUNCTIONS */
LOCAL	float	MPOLY_CJ_state(Locstate,int,Locstate,int,int);
LOCAL	void	MPOLY_progress_state(float,Locstate,Locstate,float);
LOCAL	void	MPOLY_fprint_combustion_params(FILE*,Gas_param*);
#endif /* defined(COMBUSTION_CODE) */

	/* METHOD OF CHARACTERISTIC FUNCTIONS FOR W_SPEED */
LOCAL	void	MPOLY_neumann_riem_inv_moc(float*,Locstate,float,float,Locstate,
	                                   SIDE,Locstate,float,float*,Front*);
LOCAL	void	MPOLY_shock_ahead_state_riem_inv_moc(float*,Locstate,Locstate,
	                                             Locstate,Locstate,
	                                             Locstate,float,float,
	                                             float,float,float*,float*,
						     int,float,Front*);
LOCAL	bool	MPOLY_shock_moc_plus_rh(float*,Locstate,Locstate,Locstate,
					Locstate,float,float*,float*,int,
					Front*);

	/* INITIALIZATION UTILITY FUNCTIONS */
LOCAL	void	MPOLY_prompt_for_state(Locstate,int,Gas_param*,const char*);
LOCAL	void	MPOLY_prompt_for_thermodynamics(Locstate,Gas_param*,
						const char*);
LOCAL	void	MPOLY_fprint_EOS_params(FILE*,Gas_param*);
LOCAL	void	MPOLY_read_print_EOS_params(INIT_DATA*,const IO_TYPE*,
                                            Gas_param*);
LOCAL	EOS*	MPOLY_free_EOS_params(EOS*);
LOCAL	void	MPOLY_prompt_for_EOS_params(INIT_DATA*,Gas_param*,
					    const char*,const char*);

	/* Problem Type Specific Initialization Functions */
LOCAL	float	MPOLY_RT_RS_f(float,Locstate,float,float,float);
LOCAL	void	MPOLY_RT_single_mode_perturbation_state(Locstate,float*,float,
							Locstate,float,float,
							MODE*,float);
LOCAL	void	MPOLY_KH_single_mode_state(Locstate,float*,float,Locstate,
					   float,float,float,MODE*);
LOCAL	void	MPOLY_compute_isothermal_stratified_state(Locstate,float,float,
							  Locstate);
LOCAL	void	MPOLY_compute_isentropic_stratified_state(Locstate,float,float,
							  Locstate);
LOCAL   float   MPOLY_MGamma(Locstate);

LOCAL	void	set_eos_function_hooks(EOS*);

LOCAL	float	*Vec_Gas_Gamma(Vec_Gas*,int,int);
LOCAL	float	Gamma(Locstate);
LOCAL	float	Mu(Locstate);
LOCAL	float	Coef1(Locstate);
LOCAL	float	Coef2(Locstate);
LOCAL	float	Coef3(Locstate);
LOCAL	float	Coef4(Locstate);
LOCAL	float	Coef5(Locstate);
LOCAL	float	Coef6(Locstate);
LOCAL	float	Coef7(Locstate);
LOCAL	float	Coef8(Locstate);
LOCAL	float	molecular_weight(Locstate);
LOCAL	float	R(Locstate);


typedef struct {
	float ca, cb, dS, alpha, beta, gam, c2, c4, psi, eta;
	float ua, rad, delta;
	bool cyl_coord;
} MPOLY_S_MPRH;

	/* Polytropic specific utility functions */
LOCAL	bool	mpoly_fS(float,float*,POINTER);

EXPORT	EOS	*set_MPOLY_eos(
	EOS	*eos)
{
	if (eos == NULL)
		scalar(&eos,sizeof(MPOLY_EOS));
	(void) set_GENERIC_eos(eos);
	set_eos_function_hooks(eos);
	return eos;
}

LOCAL	void	set_eos_function_hooks(
	EOS *eos)
{
	/* PRIMARY THERMODYNAMIC FUNCTIONS */
	eos->_internal_energy = MPOLY_internal_energy;
	eos->_pressure = MPOLY_pressure;
	eos->_sound_speed_squared = MPOLY_sound_speed_squared;
	eos->_acoustic_impedance_squared = MPOLY_acoustic_impedance_squared;
	eos->_specific_internal_energy = MPOLY_specific_internal_energy;

	/* SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS */
	eos->_specific_enthalpy = MPOLY_specific_enthalpy;
	eos->_temperature = MPOLY_temperature;
	eos->_entropy = MPOLY_entropy;
	eos->_adiabatic_gamma = MPOLY_adiabatic_gamma;
	eos->_gruneisen_gamma = MPOLY_gruneisen_gamma;
	eos->_fundamental_derivative = MPOLY_fundamental_derivative;
	eos->_C_V = MPOLY_C_V;
	eos->_C_P = MPOLY_C_P;
	eos->_K_S = MPOLY_K_S;
	eos->_K_T = MPOLY_K_T;

	/* VECTORIZED THERMODYNAMIC FUNCTIONS */
	eos->_single_eos_load_pressure_and_sound_speed2 =
		MPOLY_single_eos_load_pressure_and_sound_speed2;
	eos->_single_eos_load_pressure_and_gammas =
		MPOLY_single_eos_load_pressure_and_gammas;
	eos->_single_eos_load_pressure = MPOLY_single_eos_load_pressure;
	eos->_single_eos_load_sound_speed2 = MPOLY_single_eos_load_sound_speed2;

	/* RIEMANN SOLUTIONS UTILITY FUNCTIONS */
	/* Purely Thermodynamic Hugoniot Functions */
	eos->_dens_Hugoniot = MPOLY_dens_Hugoniot;
	eos->_state_w_pr_on_Hugoniot = MPOLY_state_w_pr_on_Hugoniot;
	eos->_state_w_mf_sqr_on_Hugoniot = MPOLY_state_w_mf_sqr_on_Hugoniot;

	/* Velocity Related Hugoniot Functions */
	eos->_pr_normal_vel_wave_curve = MPOLY_pr_normal_vel_wave_curve;

	/* Purely Thermodynamic Adiabatic Wave Curve Functions */
	eos->_dens_rarefaction = MPOLY_dens_rarefaction;
	eos->_pressure_rarefaction = MPOLY_pressure_rarefaction;
	eos->_state_on_adiabat_with_pr = MPOLY_state_on_adiabat_with_pr;
	eos->_state_on_adiabat_with_dens = MPOLY_state_on_adiabat_with_dens;

	/* General Wave Curve Functions */
	eos->_mass_flux = MPOLY_mass_flux;
	eos->_mass_flux_squared = MPOLY_mass_flux_squared;

	/* Functions for the Evaluation of Riemann Solutions */
	eos->_oned_fan_state = MPOLY_oned_fan_state;

	/* Functions to Compute Riemann Solutions */
	eos->_riemann_wave_curve = MPOLY_riemann_wave_curve;
	eos->_set_state_for_find_mid_state = MPOLY_set_state_for_find_mid_state;
	eos->_eps_for_Godunov = MPOLY_eps_for_Godunov;
	eos->_initialize_riemann_solver = MPOLY_initialize_riemann_solver;

	/* TWO DIMENSIONAL RIEMANN SOLUTION UTILTITY FUNCTIONS */
	eos->_steady_state_wave_curve = MPOLY_steady_state_wave_curve;
	eos->_pressure_at_sonic_point = MPOLY_pressure_at_sonic_point;
	eos->_pr_at_max_turn_angle = MPOLY_pr_at_max_turn_angle;
	eos->_state_in_prandtl_meyer_wave = MPOLY_state_in_prandtl_meyer_wave;

#if defined(COMBUSTION_CODE)
	/* DETONATION SPECIFIC UTILITY FUNCTIONS */
	eos->_CJ_state = MPOLY_CJ_state;
	eos->_progress_state = MPOLY_progress_state;
	eos->_fprint_combustion_params = MPOLY_fprint_combustion_params;
#endif /* defined(COMBUSTION_CODE) */

	/* METHOD OF CHARACTERISTIC FUNCTIONS FOR W_SPEED */
	eos->_neumann_riem_inv_moc = MPOLY_neumann_riem_inv_moc;
	eos->_shock_ahead_state_riem_inv_moc =
					MPOLY_shock_ahead_state_riem_inv_moc;
	eos->_shock_moc_plus_rh = MPOLY_shock_moc_plus_rh;

	/* INITIALIZATION UTILITY FUNCTIONS */
	eos->_prompt_for_state = MPOLY_prompt_for_state;
	eos->_prompt_for_thermodynamics = MPOLY_prompt_for_thermodynamics;
	eos->_fprint_EOS_params = MPOLY_fprint_EOS_params;
	eos->_read_print_EOS_params = MPOLY_read_print_EOS_params;
	eos->_free_EOS_params = MPOLY_free_EOS_params;
	eos->_prompt_for_EOS_params = MPOLY_prompt_for_EOS_params;

	/* Problem Type Specific Initialization Functions */
	eos->_RT_RS_f = MPOLY_RT_RS_f;
	eos->_RT_single_mode_perturbation_state =
				MPOLY_RT_single_mode_perturbation_state;
	eos->_KH_single_mode_state = MPOLY_KH_single_mode_state;
	eos->_compute_isothermal_stratified_state = MPOLY_compute_isothermal_stratified_state;
	eos->_compute_isentropic_stratified_state = MPOLY_compute_isentropic_stratified_state;
        eos->_specific_enthalpy_species = MPOLY_specific_enthalpy_species;
        eos->_dynamic_viscosity = MPOLY_dynamic_viscosity;
        eos->_dynamic_thermal_conductivity = MPOLY_dynamic_thermal_conductivity;
        eos->_dynamic_viscosity_thermalconduct = MPOLY_dynamic_viscosity_thermalconduct;
        eos->_MGamma = MPOLY_MGamma;
}


/***************PRIMARY THERMODYNAMIC FUNCTIONS ****************************/

/*
*			MPOLY_internal_energy():
*
*	Returns the internal energy per unit volume of a state.
*/

LOCAL	float	MPOLY_internal_energy(
	Locstate state)
{
	float gam;

	switch (state_type(state)) 
	{
	case GAS_STATE:
	    return	Energy(state) - kinetic_energy(state);

	case EGAS_STATE:
	    return	Energy(state)*Dens(state);

	case VGAS_STATE:
	    return Dens(state) * Int_en(state);
		
	case TGAS_STATE:
	    return	Press(state) / (Gamma(state) - 1.0);

	case FGAS_STATE:
	    gam = Gamma(state);
	//    return Dens(state)*Temperature(state)/(gam - 1.0);
	    return R(state)*Dens(state)*Temperature(state)/(gam - 1.0);

	default:
	    screen("ERROR: in MPOLY_internal_energy(), no such state type\n");
	    clean_up(ERROR);
	}
	return ERROR_FLOAT;
}		/*end MPOLY_internal_energy*/


/*
*			MPOLY_pressure():
*
*	Returns the thermodynamic pressure of a state.
*
*				     dE  |
*			     P = -  ---- |
*		                     dV  |S
*
*	Where E = specific internal energy,  V = specific volume,  and
*	S = specific entropy.
*/

LOCAL	float	MPOLY_pressure(
	Locstate state)
{
	float pr, rho;

	if (is_obstacle_state(state))
	    return HUGE_VAL;
	rho = Dens(state);
	switch (state_type(state)) 
	{
	case GAS_STATE:
	    pr = (Gamma(state) - 1.0) * (Energy(state) - kinetic_energy(state));
	    break;

	case EGAS_STATE:
	    pr = (Gamma(state) - 1.0) * Energy(state) * rho;
	    break;

	case FGAS_STATE:
	    pr = R(state)*Temperature(state)*Dens(state);
	    break;

	case TGAS_STATE:
	case VGAS_STATE:
	    pr = Press(state);
	    break;

	default:
	    screen("ERROR in MPOLY_pressure(), no such state type\n");
	    clean_up(ERROR);
	}
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	pr = max(pr,Min_pressure(state));
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	return pr;
}		/*end MPOLY_pressure*/


/*
*			MPOLY_sound_speed_squared():
*
*	Returns the square of the local sound speed of the state.
*
*                        2   dP  |
*			c = ---- |
*                           drho |S
*/

LOCAL	float	MPOLY_sound_speed_squared(
	Locstate state)
{
	float c2;
	if (state_type(state) == VGAS_STATE)
	{
	    float c = Sound_speed(state);
	    return c*c;
	}
	c2 =  Gamma(state)*pressure(state)/Dens(state);
	return c2;
}		/*end MPOLY_sound_speed_squared*/


/*
*		MPOLY_acoustic_impedance_squared():
*
*	Returns the square of the local acoustic impedance of the state.
*
*                        2     dP  |
*			i = - ---- |
*                              dV  |S
*/

LOCAL	float	MPOLY_acoustic_impedance_squared(
	Locstate state)
{
	float i2;
	if (state_type(state) == VGAS_STATE)
	{
	    float i = Dens(state)*Sound_speed(state);
	    return i*i;
	}
	i2 = Gamma(state)*pressure(state)*Dens(state);
	return i2;
}		/*end MPOLY_acoustic_impedance_squared*/

/*
*			MPOLY_specific_internal_energy():
*
*	Returns the specific internal energy = internal energy per unit
*	mass of the state.
*/

LOCAL	float	MPOLY_specific_internal_energy(
	Locstate state)
{
	float rho = Dens(state);

	switch (state_type(state))
	{
	case GAS_STATE:
	    return	(Energy(state) - kinetic_energy(state))/rho;

	case EGAS_STATE:
	    return	Energy(state);

	case TGAS_STATE:
	    return	Coef6(state)* Press(state) / rho;

	case FGAS_STATE:
	    return R(state)*Temperature(state)*Coef6(state);
	
	case VGAS_STATE:
	    return Int_en(state);

	default:
	    screen("ERROR in MPOLY_specific_internal_energy(), "
	           "no such state type\n");
	    clean_up(ERROR);
	    break;
	}
	return ERROR_FLOAT;

}		/*end MPOLY_specific_internal_energy*/


/***************END PRIMARY THERMODYNAMIC FUNCTIONS ************************/
/***************SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS ***********/

/*
*			MPOLY_specific_enthalpy():
*
*	This function computes the specific enthalpy of the given state.
*
*			H = E + P*V
*
*	E = specific internal energy, P = pressure, V = specific volume.
*
*/

LOCAL	float	MPOLY_specific_enthalpy(
	Locstate state)
{
#if defined(VERBOSE_PLUS_GAS)
	if (state_type(state) == VGAS_STATE)
	    return Enthalpy(state);
#endif /* defined(VERBOSE_PLUS_GAS) */

	return Coef7(state) * pressure(state) / Dens(state);
}		/*end MPOLY_specific_enthalpy*/


/*
*			MPOLY_temperature():
*
*	Returns the thermodynamic temperature of a state.
*
*                            dE |
*			T = --- |
*                            dS |V
*/

LOCAL	float	MPOLY_temperature(
	Locstate state)
{
	if (state_type(state) == FGAS_STATE) return Temperature(state);

	return pressure(state)/(Dens(state)*R(state));
}		/*end MPOLY_temperature*/

/*
*			MPOLY_entropy():
*
*	Returns the specific entropy of a state.
*/

LOCAL	float	MPOLY_entropy(
	Locstate state)
{
	if (state_type(state) == VGAS_STATE)
	    return Entropy(state);
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	if (Dens(state) < Vacuum_dens(state))
	    return 0.0;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	return (log(pressure(state)) - Gamma(state)*log(Dens(state))) *
			Coef6(state)*R(state);
}		/*end MPOLY_entropy*/

/*
*			MPOLY_adiabatic_gamma():
*
*	Returns the dimensionless sound speed
*
*		gamma = - d(log P)/d(log V) | .
*					     S
*	As usual P = thermodynamic pressure,  V = specific volume
*	and S = specific entropy.
*/

LOCAL	float	MPOLY_adiabatic_gamma(
	Locstate state)
{
	return Gamma(state);
}		/*end MPOLY_adiabatic_gamma*/


/*
*			MPOLY_gruneisen_gamma():
*
*	Returns the dimensionless Gruneisen exponent
*
*
*                                                 dP/dE |
*		GAMMA = - d(log T)/d(log V) |  =  -----  V
*                                            S     rho
*
*	As usual P = thermodynamic pressure,  V = specific volume
*	rho = density, E = specific internal energy,
*	and  S = specific entropy.
*
*
*/

LOCAL	float	MPOLY_gruneisen_gamma(
	Locstate state)
{
	return Gamma(state) - 1.0;
}		/*end MPOLY_gruneisen_gamma*/

/*
*			MPOLY_fundamental_derivative():
*
*	Returns the fundamental derivative of gas dynamics for the state.
*	This quantity is defined by the formula
*
*			    2      2
*		           d P / dV  |
*                                    |S
*             G = -0.5 V -----------------
*                          dP / dV |
*                                  |S
*
*	Where P is the thermodynamic pressure,  V is the specific volume
*	and S is the specific entropy.  Both derivatives are taken at
*	constant S.
*/

LOCAL	float	MPOLY_fundamental_derivative(
	Locstate state)
{
	return	0.5*(Gamma(state) + 1.0);
}		/*end MPOLY_fundamental_derivative*/

/*
*			MPOLY_C_V():
*
*	Specific heat at constant volume.
*
*                        dS  |
*		C_V = T ---- |
*                        dT  | V
*/

LOCAL	float	MPOLY_C_V(
	Locstate state)
{
	return R(state)/(Gamma(state) - 1.0);
}	/* end MPOLY_C_V */

/*
*			MPOLY_C_P():
*
*	Specific heat at constant pressure.
*
*
*                        dS  |
*		C_P = T ---- |
*                        dT  | P
*/

LOCAL	float	MPOLY_C_P(
	Locstate state)
{
	return R(state)*Gamma(state)/(Gamma(state) - 1.0);
}	/* end MPOLY_C_P */

/*
*			K_S():
*
*	Isentropic compressibility.
*
*                        1   dV  |
*		K_S = - --- ---- |
*                        V   dP  | S
*/

LOCAL	float	MPOLY_K_S(
	Locstate state)
{
	return 1.0/(Gamma(state)*pressure(state));
}	/* end MPOLY_K_S */

/*
*			MPOLY_K_T():
*
*	Isothermal compressibility.
*
*                        1   dV  |
*		K_T = - --- ---- |
*                        V   dP  | T
*/

LOCAL	float	MPOLY_K_T(
	Locstate state)
{
	return 1.0/pressure(state);
}	/* end MPOLY_K_T */



/***************END SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS *******/
/***************VECTORIZED THERMODYNAMIC FUNCTIONS *************************/

/*
*		MPOLY_single_eos_load_pressure_and_sound_speed2():
*
*	Loads a vector of pressures and sound speeds into the
*	appropriate fields of the Vec_Gas structure.
*
*	NOTE :
*	Only callable via the function wrapper load_pressure_and_sound_speed.
*	Assumes that the specific internal energy field is set.
*	This function could be written in terms of the locstate
*	thermodynamic functions,  but is provided in primitive
*	form for increased efficiency of execution time code.
*/

LOCAL	void	MPOLY_single_eos_load_pressure_and_sound_speed2(
	Vec_Gas *vst,
	int     offset,
	int     vsize)
{
	float	*rho = vst->rho + offset;
	float	*p = vst->p + offset, *c2 = vst->c2 + offset;
	int     k;
	float	*gm = Vec_Gas_Gamma(vst,offset,vsize);

	if (Vec_Gas_field_set(vst,re))
	{
	    float *re = vst->re + offset;
	    for (k = 0; k < vsize; ++k)
	        p[k] = (gm[k]-1.0)*re[k];
	}
	else
	{
	    float *e = vst->e + offset;
	    for (k = 0; k < vsize; ++k)
	        p[k] = (gm[k]-1.0)*rho[k]*e[k];
	}
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	limit_pressure(p,vst->min_pressure + offset,vsize);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	for (k = 0; k < vsize; ++k)
	    c2[k] = gm[k]*p[k]/rho[k];
	if (vst->FD != NULL)
	{
	    float *FD = vst->FD + offset;
	    for (k = 0; k < vsize; ++k)
	        FD[k] = 0.5*(gm[k] + 1.0);
	}
}		/*end MPOLY_single_eos_load_pressure_and_sound_speed2*/


/*
*		MPOLY_single_eos_load_pressure_and_gammas():
*
*	Loads the pressure, adiabatic exponent, and Gruneisen
*	coefficient uni_arrays of the Vec_Gas state vst.
*	This function assumes that the specific internal energy
*	uni_array vst->e is already loaded.
*
*	NOTE :
*	Only callable via the function wrapper load_pressure_and_sound_speed.
*	Assumes that the specific internal energy field is set.
*	This function could be written in terms of the locstate
*	thermodynamic functions,  but is provided in primitive
*	form for increased efficiency of execution time code.
*/

LOCAL	void	MPOLY_single_eos_load_pressure_and_gammas(
	Vec_Gas *vst,
	int     offset,
	int     vsize)
{
	float *rho = vst->rho + offset;
	float *p = vst->p + offset;
	float *c2 = vst->c2 + offset, *GAM = vst->GAM + offset;
	float *gm = Vec_Gas_Gamma(vst,offset,vsize);
	int   k;

	if (Vec_Gas_field_set(vst,re))
	{
	    float *re = vst->re + offset;
	    for (k = 0; k < vsize; ++k)
	    {
	        GAM[k] = gm[k] - 1.0;
	        p[k] = (gm[k]-1.0)*re[k];
	    }
	}
	else
	{
	    float *e = vst->e + offset;
	    for (k = 0; k < vsize; ++k)
	    {
	        GAM[k] = gm[k] - 1.0;
	        p[k] = (gm[k]-1.0)*rho[k]*e[k];
	    }
	}

	if (vst->FD != NULL)
	{
	    float *FD = vst->FD + offset;
	    for (k = 0; k < vsize; ++k)
	        FD[k] = 0.5*(gm[k] + 1.0);
	}
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	limit_pressure(p,vst->min_pressure + offset,vsize);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	for (k = 0; k < vsize; ++k)
	    c2[k] = gm[k]*p[k]/rho[k];
}		/*end MPOLY_single_eos_load_pressure_and_gammas*/

/*
*			MPOLY_single_eos_load_pressure():
*
*	Loads a vector of pressures into the appropriate field of the 
*	Vec_Gas structure.
*
*	NOTE :
*	Only callable via the function wrapper load_pressure_and_sound_speed.
*	Assumes that the specific internal energy field is set.
*	This function could be written in terms of the locstate
*	thermodynamic functions,  but is provided in primitive
*	form for increased efficiency of execution time code.
*/

LOCAL	void	MPOLY_single_eos_load_pressure(
	Vec_Gas *vst,
	int offset,
	int vsize)
{
	float *p = vst->p + offset;
	float *rho = vst->rho + offset;
	float *gm = Vec_Gas_Gamma(vst,offset,vsize);
	int   k;

	if (Vec_Gas_field_set(vst,re))
	{
	    float *re = vst->re + offset;
	    for (k = 0; k < vsize; ++k)
	        p[k] = (gm[k]-1.0)*re[k];
	}
	else
	{
	    float *e = vst->e + offset;
	    for (k = 0; k < vsize; ++k)
	        p[k] = (gm[k]-1.0)*rho[k]*e[k];
	}
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	limit_pressure(p,vst->min_pressure + offset,vsize);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
}		/*end MPOLY_single_eos_load_pressure*/

/*
*			MPOLY_single_eos_load_sound_speed2():
*
*	Loads a vector of sound speeds into the appropriate field of the 
*	Vec_Gas structure.
*
*	NOTE :
*	Only callable via the function wrapper load_pressure_and_sound_speed.
*	Assumes that the specific internal energy field is set.
*	This function could be written in terms of the locstate
*	thermodynamic functions,  but is provided in primitive
*	form for increased efficiency of execution time code.
*/

LOCAL	void	MPOLY_single_eos_load_sound_speed2(
	Vec_Gas *vst,
	int offset,
	int vsize)
{
	float *e = vst->e + offset;
	float *c2 = vst->c2 + offset;
	float *gm = Vec_Gas_Gamma(vst,offset,vsize);
	int   k;

	for (k = 0; k < vsize; ++k)
	    c2[k] = gm[k]*(gm[k]-1.0)*e[k];
	if (vst->FD != NULL)
	{
	    float *FD = vst->FD + offset;
	    for (k = 0; k < vsize; ++k)
	        FD[k] = 0.5*(gm[k] + 1.0);
	}
}		/*end MPOLY_single_eos_load_sound_speed2*/


/***************END VECTORIZED THERMODYNAMIC FUNCTIONS *********************/

/***************RIEMANN SOLUTIONS UTILITY FUNCTIONS ************************/

/***************Purely Thermodynamic Hugoniot Functions*********************/

/*
*			MPOLY_dens_Hugoniot():
*
*	Given the state st0 on one side of an oblique shock and the pressure
*	p1 on the other side, this function returns the density rho1 of the
*	state with pressure p1.  Rho1 is found by solving the Hugoniot relation
*
*		(p1 + p0)*(1/rho0 - 1/rho1) = 2*(e1 - e0)
*
*	where e0 and e1 are the specific internal energies of the two
*	respective states.  For a given equation of state the specific
*	internal energy can be expressed as a function of the
*	pressure and density.  Thus the above equation can be solved to
*	give rho1 as a function of st0 and p1.
*
*
*	Reference: Courant and Friedrichs page 302 ff.
*/


LOCAL	float	MPOLY_dens_Hugoniot(
	float p1,
	Locstate st0)
{
	float	p0, c4;

	p0 = pressure(st0);
	c4 = Coef4(st0);
	return Dens(st0)*(p1 + p0*c4)/ (p0 + p1*c4);
}		/*end MPOLY_dens_Hugoniot*/


/*
*			MPOLY_state_w_pr_on_Hugoniot():
*
*	Given the state st0 on one side of an oblique shock and the pressure
*	p1 on the other side, this function returns the thermodynamic variables
*	of the state with pressure p1 (density and internal energy for a
*	GAS_STATE, pressure and density for a TGAS_STATE).  Rho1 is found by
*	solving the Hugoniot relation
*
*		e1 - e0 - 0.5*(p1 + p0)*(V0 - V1) = 0
*
*	where e0 and e1 are the specific internal energies of the two
*	respective states.  For a given equation of state the specific
*	internal energy can be expressed as a function of the
*	pressure and density.  Thus the above equation can be solved to
*	give rho1 and e1 as a function of st0 and p1.  The internal
*	energy is then given by E1 = r1 * e1.
*
*	IMPORTANT NOTE:
*		If stype1 == GAS_STATE the energy in st1 is
*		the internal energy.  The kinetic energy must
*		be added separately.  The reason for this is
*		that this function is a purely theromdynamic
*		function and is independent of the state
*		velocities.
*
*	Reference: Courant and Friedrichs page 302 ff.
*/

LOCAL	void	MPOLY_state_w_pr_on_Hugoniot(
	Locstate st0,
	float p1,
	Locstate st1,
	int stype1)
{
	float	p0;
	float	rho0;
	float	c4;

	p0 = pressure(st0);
	rho0 = Dens(st0);
	c4 = Coef4(st0);
	zero_state_velocity(st1,Params(st0)->dim);
	Set_params(st1,st0);
	set_type_of_state(st1,stype1);
	Dens(st1) = rho0*(p1 + p0*c4) / (p0 + p1*c4);
	switch(stype1)
	{
	case TGAS_STATE:
		Press(st1) = p1;
		break;
	case GAS_STATE:
		Energy(st1) = p1*Coef6(st0);
		break;
	case EGAS_STATE:
		Energy(st1) = p1*Coef6(st0)/Dens(st1);
		break;
	case FGAS_STATE:
		Temperature(st1) = (p1/Dens(st1))/R(st1);
		break;
	case VGAS_STATE:
		Press(st1) = p1;
		set_type_of_state(st1,TGAS_STATE);
		set_state(st1,VGAS_STATE,st1);
		break;
	default:
		screen("ERROR in state_w_pr_on_Hugoniot(), ");
		screen("Unknown state type %d\n",stype1);
		clean_up(ERROR);
	}
}		/*end MPOLY_state_w_pr_on_Hugoniot*/

/*
*			MPOLY_state_w_mf_sqr_on_Hugoniot():
*
*	Given the state st0 on one side of an oblique shock and the square
*	of the mass flux across the shock, this function returns the
*	thermodynamic variables of the state on the opposite side of the shock.
*
*	By definition the square of the mass flux across a shock is given by
*
*			mf_sqr = (p1 - p0) / (V0 - V1)
*
*	where pi and Vi denote the pressure and specific volume on
*	of the two states on either side of the shock.
*
*	IMPORTANT NOTE:
*		If stype1 == GAS_STATE the energy in st1 is
*		the internal energy.  The kinetic energy must
*		be added separately.  The reason for this is
*		that this function is a purely theromdynamic
*		function and is independent of the state
*		velocities.
*
*/

LOCAL	bool	MPOLY_state_w_mf_sqr_on_Hugoniot(
	Locstate st0,
	float m2,
	Locstate st1,
	int stype1)
{
	float	p0, p1;
	float 	r0;
	float	c4;
	int	k;

	p0 = pressure(st0);
	r0 = Dens(st0);
	c4 = Coef4(st0);
	zero_state_velocity(st1,Params(st0)->dim);
	Set_params(st1,st0);
	set_type_of_state(st1,stype1);
	p1 = Coef5(st0)*m2/r0 - c4*p0;
	Dens(st1) = Dens(st0)*(p1 + p0*c4)/(p0 + p1*c4);
        for(k = 0; k < Params(st0)->n_comps; k++)
            pdens(st1)[k] = (pdens(st0)[k]/Dens(st0))*Dens(st1);
	switch(stype1)
	{
	case TGAS_STATE:
		Press(st1) = p1;
		break;
	case GAS_STATE:
		Energy(st1) = p1 * Coef6(st1);
		break;
	case EGAS_STATE:
		Energy(st1) = p1 * Coef6(st1)/Dens(st1);
		break;
	case FGAS_STATE:
		Temperature(st1) = (p1/Dens(st1))/R(st1);
		break;
	case VGAS_STATE:
		Press(st1) = p1;
		set_type_of_state(st1,TGAS_STATE);
		set_state(st1,VGAS_STATE,st1);
		break;
	default:
		screen("ERROR in state_w_mf_sqr_on_Hugoniot()\n");
		screen("Unknown state type %d\n",stype1);
		clean_up(ERROR);
	}
	return FUNCTION_SUCCEEDED;
}		/*end MPOLY_state_w_mf_sqr_on_Hugoniot*/


/***************End Purely Thermodynamic Hugoniot Functions*****************/
/***************Velocity Related Hugoniot Functions*************************/

/*
*			MPOLY_pr_normal_vel_wave_curve():
*
*	Computes the pressure on the forward Riemann wave curve given the
*	velocity difference across the wave.
*
*	If du > 0 returns the solution to the system:
*
*                 2
*		du   = (p - p0)*(V0 - V)
*		de   = 0.5*(p + p0)*(V0 - V)
*
*	if du < 0 returns the solution to the system:
*
*		     /p    dP   |
*		du = \    ----  |
*		     /p0  rho c |S
*/

LOCAL	float	MPOLY_pr_normal_vel_wave_curve(
	float du,	/* normal velocity change across shock = (u1 - u0)*/
	Locstate st0)
{
	float b, c, disc, p0;

	p0 = pressure(st0);

	if (du > 0.0)
	{
		c = Gamma(st0)*Dens(st0)*sqr(du) / p0;
		b = c / (1.0 + Coef4(st0));
		disc = (b*b + 4.0*c);
		return p0 * (1.0 + 0.5*(b + sqrt(disc)));
	}
	else if (du < 0.0)
	{
		float y, p1;

		y = 1.0 + Coef2(st0)*du/sound_speed(st0);
		if (y <= 0.0)
			return Min_pressure(st0);
		p1 = pow(y,1.0/Coef3(st0));
		if (p1 < Min_pressure(st0))
			p1 = Min_pressure(st0);
		return p1;
	}
	else
		return p0;
}		/*end MPOLY_pr_normal_vel_wave_curve*/


/***************End Velocity Related Hugoniot Functions*********************/
/***************Purely Thermodynamic Adiabatic Wave Curve Functions*********/

/*	
*			MPOLY_dens_rarefaction():
*
*	Given the state st0 and the pressure on the other side of
*	a simple wave in steady irrotational flow, this
* 	function returns the density on the other side.
*
*	The answer is give by the solution of the ordinary differential
*	equation
*
*		dh/dP = V,  h(p0) = h0;
*
*	where h is the specific enthalpy,  and the derivatives are taken
*	at constant entropy.
*/

LOCAL	float	MPOLY_dens_rarefaction(
	float p1,
	Locstate st0)
{
	return Dens(st0) * pow(p1/pressure(st0),1.0/Gamma(st0));
}		/*end MPOLY_dens_rarefaction*/

/*	
*			MPOLY_pressure_rarefaction():
*
*	Given the state st0 and the density on the other side of
*	a simple wave in steady irrotational flow, this
* 	function returns the pressure on the other side.
*
*	The answer is give by the solution of the ordinary differential
*	equation
*
*		de/dV = -P,  e(V0) = e0;
*
*	where e is the specific internal energy,  and the derivatives are taken
*	at constant entropy.
*/

LOCAL	float	MPOLY_pressure_rarefaction(
	float rho1,
	Locstate st0)
{
	return pressure(st0) * pow(rho1/Dens(st0),Gamma(st0));
}		/*end MPOLY_pressure_rarefaction*/


/*	
*			MPOLY_state_on_adiabat_with_pr():
*
*	Given the state st0 and the pressure on the other side of
*	a simple wave in steady irrotational flow, this function returns
*	the thermodynamic variable on the other side.
*
*	IMPORTANT NOTE:
*		If stype1 == GAS_STATE the energy in st1 is
*		the internal energy.  The kinetic energy must
*		be added separately.  The reason for this is
*		that this function is a purely theromdynamic
*		function and is independent of the state
*		velocities.
*
*/

LOCAL	void	MPOLY_state_on_adiabat_with_pr(
	Locstate st0,
	float p1,
	Locstate st1,
	int stype1)
{
	int	k;
	zero_state_velocity(st1,Params(st0)->dim);
	Set_params(st1,st0);
	set_type_of_state(st1,stype1);
	Dens(st1) = Dens(st0)*pow(p1/pressure(st0),1.0/Gamma(st0));
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            for(k = 0; k < Params(st0)->n_comps; k++)
            {
                pdens(st1)[k] = (pdens(st0)[k]/Dens(st0))*Dens(st1);
                /* New 051005 */
                if(fabs(pdens(st1)[k]) < 10.0*MACH_EPS && pdens(st1)[k] < 0.0)
                    pdens(st1)[k] = 0.0;
                /* End of New 051005 */
            }
        }
	switch(stype1)
	{
	case TGAS_STATE:
		Press(st1) = p1;
		break;
	case GAS_STATE:
		Energy(st1) = p1 * Coef6(st1);
		break;
	case EGAS_STATE:
		Energy(st1) = p1*Coef6(st1)/Dens(st1);
		break;
	case FGAS_STATE:
		Temperature(st1) = R(st1)*p1/Dens(st1);
		break;
	case VGAS_STATE:
		Press(st1) = p1;
		set_type_of_state(st1,TGAS_STATE);
		set_state(st1,VGAS_STATE,st1);
		break;
	default:
		screen("ERROR in state_on_adiabat_with_pr()\n");
		screen("Unknown state type %d\n",stype1);
		clean_up(ERROR);
	}
}		/*end MPOLY_state_on_adiabat_with_pr*/

/*	
*			MPOLY_state_on_adiabat_with_dens():
*
*	Given the state st0 and the density on the other side of
*	a simple wave in steady irrotational flow, this	function returns
*	the pressure and internal energy on the other side.
*
*	IMPORTANT NOTES:
*		1.  If stype1 == GAS_STATE the energy in st1 is
*		the internal energy.  The kinetic energy must
*		be added separately.  The reason for this is
*		that this function is a purely theromdynamic
*		function and is independent of the state
*		velocities.
*
*		2.  Dens(st1) cannot be set to rho1 before the evaluation of
*		the pressure of st0.  This allows this function to work
*		even in the case were st0 = st1 (ie they both point to the
*		same area in storage).
*/

LOCAL	void	MPOLY_state_on_adiabat_with_dens(
	Locstate st0,
	float rho1,
	Locstate st1,
	int stype1)
{
	float p1; 	/* pressure of the answer state */
	int	i, nc = Num_gas_components(st0);

	Set_params(st1,st0);
	for (i = 0; i < nc; ++i)
		pdens(st1)[i] = rho1*pdens(st0)[i]/Dens(st0);
	zero_state_velocity(st1,Params(st0)->dim);
	set_type_of_state(st1,stype1);
	p1 =  pressure(st0) * pow(rho1/Dens(st0),Gamma(st0));
	Dens(st1) = rho1;
	switch(stype1)
	{
	case TGAS_STATE:
		Press(st1) = p1;
		break;
	case GAS_STATE:
		Energy(st1) = p1 * Coef6(st0);
		break;
	case EGAS_STATE:
		Energy(st1) = p1*Coef6(st0)/rho1;
		break;
	case FGAS_STATE:
		Temperature(st1) = p1*R(st1)/rho1;
		break;
	case VGAS_STATE:
		Press(st1) = p1;
		set_type_of_state(st1,TGAS_STATE);
		set_state(st1,VGAS_STATE,st1);
		break;
	default:
		screen("ERROR in state_on_adiabat_with_dens()\n");
		screen("Unknown state type %d\n",stype1);
		clean_up(ERROR);
	}
}		/*end MPOLY_state_on_adiabat_with_dens*/

/***************End Purely Thermodynamic Adiabatic Wave Curve Functions*****/
/***************General Wave Curve Functions********************************/

/*
*			MPOLY_mass_flux():
*
*	Returns the mass flux across a wave.
*
*				 
*		     | (P - P0) |
*		m  = | -------  |
*		     | (U - U0) |
*
*	Where 
*		P0 = pressure ahead of the shock
*		U0 = velocity ahead of the shock
*		P = pressure behind the shock
*		U = velocity behind the shock
*
*/

LOCAL	float	MPOLY_mass_flux(
	float p,
	Locstate st0)
{
	float p0, rho0;
	float xi, m, i0;

	p0 = pressure(st0);
	rho0 = Dens(st0);
	if (p < p0)
	{
		i0 = acoustic_impedance(st0);
		xi = p/p0;
		if ( (1.0 - xi) < EPS) return i0;
		m = Coef3(st0)*(1.0-xi)/(1.0-pow(xi,Coef3(st0)));
		return i0*m;
	}
	else
		return sqrt(rho0*(Coef1(st0)*p + Coef2(st0)*p0));
}		/*end MPOLY_mass_flux*/

/*
*			MPOLY_mass_flux_squared():
*
*	Returns the square of the mass flux across a wave.
*
*				 2
*		 2   | (P - P0) |
*		m  = | -------  |
*		     | (U - U0) |
*
*	Where 
*		P0 = pressure ahead of the shock
*		U0 = velocity ahead of the shock
*		P = pressure behind the shock
*		U = velocity behind the shock
*
*/

LOCAL	float	MPOLY_mass_flux_squared(
	float p,
	Locstate st0)
{
	float p0, rho0;
	float xi, m, i02;

	p0 = pressure(st0);
	rho0 = Dens(st0);
	if (p < p0)
	{
		i02 = acoustic_impedance_squared(st0);
		xi = p/p0;
		if ( (1.0 - xi) < EPS) return i02;
		m = Coef3(st0)*(1.0-xi)/(1.0-pow(xi,Coef3(st0)));
		return i02*m*m;
	}
	else
		return rho0*(Coef1(st0)*p + Coef2(st0)*p0);
}		/*end MPOLY_mass_flux_squared*/


/***************End General Wave Curve Functions****************************/
/***************Functions for the Evaluation of Riemann Solutions***********/

/*
*				MPOLY_oned_fan_state():
*
*	This is a utility function provided for the evaluation of states
*	in a simple wave.   Given sta, it solves for stm using the
*	equation:
*
*	                     / p_m        |            / c_m        |
*	                    /             |           /             |
*	                    \       dP    |           \        dc   |
*	    w = c_m - c_a +  \    -----   |         =  \     ------ |
*	                      \   rho c   |             \     mu^2  |
*	                      /           |             /           |
*	                     /p_a         | S = S_a    / c_a        | S = S_a
*
*	here c is the sound speed,  rho the density,  S the specific entropy,
*	p the pressure,  and mu^2 = (G - 1)/G,  where G is the fundamental
*	derivative of gas dynamics.  The returned st1 contains only
*	the thermodyanics of the state in the rarefaction fan.  In particular
*	st1 can be used to evaluate the pressure, density, and sound speed
*	of the state inside the rarefaction fan.
*	
*	Input data:
*		w = value of w as defined above
*		sta = state ahead of fan
*		stb = state behind fan
*
*	Output data:
*		stm = state inside fan
*		vacuum = 1 if stm is a vacuum,  0 otherwise
*
*	Returns the sound speed of the answer state stm.
*/

/*ARGSUSED*/
LOCAL	float	MPOLY_oned_fan_state(
	float    w,
	Locstate sta,
	Locstate stb,
	Locstate stm,
	int      stype_m,
	bool  *vacuum)
{
	float	c_a, c_m, p_a;
	float	c2, c3, c4;
	int	k;

	zero_state_velocity(stm,Params(sta)->dim);
	*vacuum = NO;

	p_a = pressure(sta);
	set_type_of_state(stm,TGAS_STATE);
	c_a = sound_speed(sta);
	c2 = 1.0 / Coef2(sta);
	c3 = 1.0 / Coef3(sta);
	c4 = Coef4(sta);
	c_m = c_a + c4*w;
	if (c_m <= 0.0)
	{
	    /* rarefaction to vacuum */
	    state_on_adiabat_with_pr(sta,Min_pressure(sta),stm,TGAS_STATE);
	    c_m = 0.0;
	    *vacuum = YES;
	}
	else
	{
	    Set_params(stm,sta);
	    Dens(stm) = Dens(sta)*pow(c_m/c_a,c2);
            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
                for(k = 0; k < Params(sta)->n_comps; k++)
                {
                    pdens(stm)[k] = (pdens(sta)[k]/Dens(sta))*Dens(stm);
                    /* New 051005 */
                    if(fabs(pdens(stm)[k]) < 10.0*MACH_EPS && pdens(stm)[k] < 0.0)
                        pdens(stm)[k] = 0.0;
                    /* End of New 051005 */
                }
            }
	    Press(stm) = p_a*pow(c_m/c_a,c3);
	}

#if !defined(UNRESTRICTED_THERMODYNAMICS)
	if (Press(stm) < Min_pressure(sta))
	{
	    state_on_adiabat_with_pr(sta,Min_pressure(sta),stm,TGAS_STATE);
	    c_m = 0.0;
	    *vacuum = YES;
	}
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	set_state(stm,stype_m,stm);
	return c_m;
}		/* end oned_fan_state*/
/***************End Functions for the Evaluation of Riemann Solutions********/



/***************Functions to Compute Riemann Solutions**********************/


/*
*			MPOLY_riemann_wave_curve():
*
*	Evalutes the forward wave family wave curve defined by
*
*		 _
*		|
*		|
*		|                                1/2
*               |   [ (Pstar  -  P0) * ( V0 - V) ]     if Pstar > P0
*		|
*		|
*	        / 
*	       /
*              \
*		\		
*		|
*               |        / Pstar     |
*               |       /            |
*               |       \      dP    |
*               |        \   ------  |		       if Pstar < P0
*               |         \   rho c  |
*               |         /          |
*               |        / P0        | S
*               |_
*
*/

LOCAL	float	MPOLY_riemann_wave_curve(
	Locstate st0,
	float pstar)
{
	float rho0 = Dens(st0), p0 = pressure(st0);
	float c1, c2, c3;

#if !defined(UNRESTRICTED_THERMODYNAMICS)
	if (pstar < Min_pressure(st0))
		pstar = Min_pressure(st0);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	c1 = Coef1(st0);
	c2 = Coef2(st0);
	c3 = Coef3(st0);

	return (pstar < p0) ?
		sound_speed(st0)*(pow(pstar/p0,c3) - 1.0)/ c2 :
		(pstar-p0)/sqrt(rho0*(c1*pstar+c2*p0));
}		/*end MPOLY_riemann_wave_curve*/


/*
*		MPOLY_set_state_for_find_mid_state():
*
*	Copies the Gas state st into the thermodynamic
*	version Tst, for some EOSs a VGas state is set.
*
*	Technical function added for enhanced performance.
*/

LOCAL	void	MPOLY_set_state_for_find_mid_state(
	Locstate Tst,
	Locstate st)
{
	set_state(Tst,TGAS_STATE,st);
}		/*end MPOLY_set_state_for_find_mid_state*/

/*
*			MPOLY_eps_for_Godunov():
*
*	Returns a tolerance to be used to determine convergence of the
*	of the Riemann solver.
*
*	Technical function added for enhanced performance.
*/

/*ARGSUSED*/
LOCAL	float	MPOLY_eps_for_Godunov(
	Locstate state,
	float pstar,
	float r_eps)
{
	return r_eps;
}		/*end MPOLY_eps_for_Godunov*/

/*
*			MPOLY_initialize_riemann_solver()
*
*	Computes the epsilons and the initial guess for the pressure
*	in the secant iteration of find_mid_state.
*
*	Technical function added for enhanced performance.
*/

/*ARGSUSED*/
LOCAL	void	MPOLY_initialize_riemann_solver(
	Locstate Tsl,
	Locstate Tsr,
	float *pstar,
	float *p_min,
	float eps,
	float *eps_u,
	float *eps_p,
	bool (*fd_md_st)(Locstate,Locstate,float,float*,float*,
			 float*,float*,float*,float*,
			 RIEMANN_SOLVER_WAVE_TYPE*,RIEMANN_SOLVER_WAVE_TYPE*))
{
	float pl, pr;
	float cl, cr, ul_tdl, ur_tdl, z;
	float c2l, c2r, c3l, c3r;
	float vl = vel(0,Tsl), vr = vel(0,Tsr);
	float dutdl;

	*eps_u = *eps_p = eps;
	pl = pressure(Tsl), pr = pressure(Tsr);
	if (Eos(Tsl) != Eos(Tsr))
	{
#if defined(UNRESTRICTED_THERMODYNAMICS)
	    *p_min = -HUGE_VAL;
#else /* defined(UNRESTRICTED_THERMODYNAMICS) */
	    *p_min = max(Min_pressure(Tsl),Min_pressure(Tsr));
#endif /* defined(UNRESTRICTED_THERMODYNAMICS) */
	    *pstar = 0.5*(pl + pr);
	    *pstar = max(*pstar,*p_min);
	    return;
	}

	c2l = Coef2(Tsl);
	c3l = Coef3(Tsl);
	c2r = Coef2(Tsr);
	c3r = Coef3(Tsr);
	cl = sound_speed(Tsl);
	cr = sound_speed(Tsr);
	ul_tdl = vl + cl/c2l;
	ur_tdl = vr - cr/c2r;
	dutdl = ul_tdl - ur_tdl;
	if (pl >= pr)
	{
	    z = (c2l*cr/(c2r*cl))*pow(pl/pr,c3l);
	    *pstar = (dutdl > 0.0) ? pl*pow(c2l*dutdl/((1.0 + z)*cl),1.0/c3l) :
				     0.5*min(pl,pr);
	}
	else
	{
	    z = (c2r*cl/(c2l*cr))*pow(pr/pl,c3r);
	    *pstar = (dutdl > 0.0) ? pr*pow(c2r*dutdl/((1.0 + z)*cr),1.0/c3r) :
				     0.5*min(pl,pr);
	}
#if defined(UNRESTRICTED_THERMODYNAMICS)
	*p_min = -HUGE_VAL;
#else /* defined(UNRESTRICTED_THERMODYNAMICS) */
	*pstar = max(*pstar,Min_pressure(Tsl));
	*pstar = max(*pstar,Min_pressure(Tsr));
	*p_min = Min_pressure(Tsl);
#endif /* defined(UNRESTRICTED_THERMODYNAMICS) */
	*pstar = max(*pstar,*p_min);
}		/*end MPOLY_initialize_riemann_solver*/


/***************End Functions to Compute Riemann Solutions******************/
/***************END RIEMANN SOLUTIONS UTILITY FUNCTIONS ********************/

/***************TWO DIMENSIONAL RIEMANN SOLUTION UTILTITY FUNCTIONS*********/

/*
*			MPOLY_steady_state_wave_curve():
*
*	Calculates the steady state wave curve through the state
*	st0 with steady state flow speed q0,  parameterized
*	by pressure.  In general the value returned is given by
*
*                       __                          --                  -- --
*                       |                           |      2             |  |
*                       |     p1/p0 - 1             |    M0              |  |
*                       |                           |                    |  |
*       theta =  arctan |-------------------  * sqrt| ---------    -    1|  |
*                       |         2                 |      2             |  |
*                       |gamma0*M0  - (p1/p0 - 1)   |    Mn              |  |
*			|                           |                    |  |
*                       --                          --                  -- --
*							        for p1 > p0
*
*	where gamma0 = adiabatic_gamma(st0), Mn = m/(rho0*c0) =
*	shock normal Mach number, and m = mass_flux across the shock.
*	and
*
*
*			/ p1
*		       /                 2
*	theta =        \     sqrt(1 - (M) ) *  dp		for p1 < p0
*			\		     ------
*			 \                        2
*			 /		     M * c * rho
*		        / p0
*
*	
*	Returns FUNCTION_SUCCEEDED on success,  FUNCTION_FAILED on failure.
*/

LOCAL	bool	MPOLY_steady_state_wave_curve(
	float p1,
	float M0sq,
	float *theta,
	Locstate st0)
{
	float	p0 = pressure(st0);
	float	dp, tmp;
	float	c4, cf8;
	float	Mnsq, cotb;
	float	A0, A1;
	float	tan_theta;
	float	gam0;

	if (p1 >= p0)		/* shock */
	{
		dp = (p1 - p0) / p0;
		gam0 = Gamma(st0);
		tan_theta = dp / (gam0*M0sq - dp);
		c4 = Coef4(st0);
		cf8 = Coef8(st0);
		Mnsq = (p1/p0 + c4) / cf8;
		cotb = M0sq/Mnsq - 1.0;
		cotb = max(0.0,cotb);
		cotb = sqrt(cotb);
		tan_theta *= cotb;
		*theta = (tan_theta > 0.0) ? atan(tan_theta) : 0.0;
		return FUNCTION_SUCCEEDED;
	}
	else		/* rarefaction */
	{
		float mu;
		float GC3 = Coef3(st0);
	
		c4 = Coef4(st0);
		mu = Mu(st0);
		if (M0sq < SONIC_MINUS_SQR)
		{
			(void) printf("WARNING in ");
			(void) printf("MPOLY_steady_state_wave_curve(), ");
			(void) printf("Subsonic state in ");
			(void) printf("rarefaction\n");
			(void) printf("Mach number = %g, (squared %g)\n",
					sqrt(M0sq),M0sq);
			return FUNCTION_FAILED;
		}
		else
		{
			M0sq = max(1.0,M0sq);
			A0 = asin(sqrt(1.0/M0sq));
			A1 = atan(mu*pow(p1/p0,GC3)/
					sqrt(1.0 + c4*(M0sq - 1.0) -
						pow(p1/p0,2.0*GC3)));
			tmp = angle(tan(A1),mu);
			*theta = -(A1 - A0 + (tmp - atan(mu/tan(A0)))/mu);
			return FUNCTION_SUCCEEDED;
		}
	}
}		/*end MPOLY_steady_state_wave_curve*/


/*
*			MPOLY_pressure_at_sonic_point():
*
*	Returns the pressure at the sonic point of the shock polar
*	through the state st0 with steady state Mach number M0.
*/

LOCAL	float	MPOLY_pressure_at_sonic_point(
	float M0sq,
	Locstate st0)
{
	float x;

	x = 0.5*(M0sq - 1.0);
	return pressure(st0) * (x + sqrt(1.0 + 2.0*Coef4(st0)*x + x*x));
}		/*end MPOLY_pressure_at_sonic_point*/


/*
*			MPOLY_pr_at_max_turn_angle():
*
*	Given st0 and the Mach number (squared) of st0 in the frame
*	of a shock, this function calculates the pressure at the point of
*	maximum turning angle on the turning angle pressure shock polar
*	through st0.
*
*	Returns FUNCTION_SUCCEEDED if sucessful,  FUNCTION_FAILED otherwise.
*/

LOCAL	bool	MPOLY_pr_at_max_turn_angle(
	float *prm,
	float M0sq,	/* Mach number of st0 in the frame of the shock */
	Locstate st0)
{
	float xi;
	float c4;

	if (M0sq < SONIC_MINUS_SQR)
	{
		(void) printf("WARNING in MPOLY_pr_at_max_turn_angle(), ");
		(void) printf("subsonic ahead state\n");
		return FUNCTION_FAILED;
	}

	M0sq = max(1.0,M0sq);
	c4 = Coef4(st0);

	xi = 0.5 * (M0sq - 4.0) + 
	sqrt(0.25 * sqr(M0sq-4.0) + 2.0*(1.0+c4)*(M0sq-1.0) );
	*prm = pressure(st0)*(xi+1.0);
	return FUNCTION_SUCCEEDED;
}		/*end MPOLY_pr_at_max_turn_angle*/

/*
*		MPOLY_state_in_prandtl_meyer_wave():
*
*	This is a utility function provided for the evaluation of states
*	in a Prandtl-Meyer wave.   Given sta, it solves for stm using the
*	equation:
*
*	                / p_m        |             /A_m                |
*	               /             |            /            2       |
*	               \   cos(a) dP |            \       csc A  dA    |
*	w = A_m - A_a + \  --------  |             \  ------------     |
*	                 \  rho c q  |             /          2    2   |
*	                 /           |B=B_a       /    (1 + mu  cot A) |B=B_a
*	                /p_a         |S=S_a      /A_a                  |S=S_a
*
*	The integrals are evaluted at constant entropy and constant
*	B = 0.5*q*q + i, where i is the specific enthalpy.  Here
*	c is the sound speed, q is the flow speed, sin(A) = c/q is the
*	              2
*	Mach angle, mu  = (G - 1)/G, and G is the fundamental derivative
*	G = 1 + rho c dc/dp|S. Note that mu may be complex for some
*	equations of state.
*
*	The returned st1 contains only
*	the thermodyanics of the state in the rarefaction fan.  In particular
*	st1 can be used to evaluate the pressure, density, and sound speed
*	of the state inside the rarefaction fan.
*	
*	Input data:
*		w = value of w as defined above
*		sta = state ahead of fan
*		A_a = Positive Mach angle of sta, sin(A_a) = c_a/q_a
*		A_b = Positive Mach angle of stb, sin(A_b) = c_b/q_b
*		stype_m = state type of stm
*
*	Output data:
*		stm = state inside fan
*
*	Returns the Mach angle of stm.
*/


/*ARGSUSED*/
LOCAL	float	MPOLY_state_in_prandtl_meyer_wave(
	float w,
	float A_a,
	Locstate sta,
	float A_b,
	Locstate stb,
	Locstate stm,
	int stype_m)
{
	float	c_a, c_m, p_a, A_m;
	float	cos_muw, sin_muw, cot_A_a;
	float	c2, c3, mu;

	zero_state_velocity(stm,Params(sta)->dim);

	set_type_of_state(stm,TGAS_STATE);
	p_a = pressure(sta);
	c_a = sound_speed(sta);
	c2 = 1.0 / Coef2(sta);
	c3 = 1.0 / Coef3(sta);
	mu = Mu(sta);
	cos_muw = cos(mu*w);
	sin_muw = sin(mu*w);
	cot_A_a = 1.0/tan(A_a);
	c_m = cos_muw + mu*cot_A_a*sin_muw;
	A_m = mu*c_m/(mu*cos_muw*cot_A_a - sin_muw);
	if (c_m <= 0.0)
	{
		/* rarefaction to vacuum */
		state_on_adiabat_with_pr(sta,Min_pressure(sta),stm,TGAS_STATE);
		c_m = 0.0;
		A_m = 0.0;
	}
	else
	{
		Set_params(stm,sta);
		Dens(stm) = Dens(sta)*pow(c_m/c_a,c2);
		Press(stm) = p_a*pow(c_m/c_a,c3);
	}

#if !defined(UNRESTRICTED_THERMODYNAMICS)
	if (Press(stm) < Min_pressure(sta))
	{
		state_on_adiabat_with_pr(sta,Min_pressure(sta),stm,TGAS_STATE);
		A_m = 0.0;
	}
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

	set_state(stm,stype_m,stm);
	return A_m;
}		/* end MPOLY_state_in_prandtl_meyer_wave*/

/***************END TWO DIMENSIONAL RIEMANN SOLUTION UTILTITY FUNCTIONS*****/

#if defined(COMBUSTION_CODE)
/***************DETONATION SPECIFIC UTILITY FUNCTIONS*********************/

/*
*			MPOLY_CJ_state():
*
* 	This routine finds the state behind a CJ_detonation.
*	The inputs are the initial state "start"
*	and the side (l_or_r, -1 or 1) we are on.
*/

/*ARGSUSED*/
LOCAL	float	MPOLY_CJ_state(
	Locstate CJ,
	int st_type_CJ,
	Locstate start,
	int l_or_r,
	int avail)
{
	screen("ERROR in MPOLY_CJ_state(), nonreactive EOS\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}	/* end MPOLY_CJ_state*/


/*
*	 		MPOLY_progress_state(): 
*
*	Finds the gas state as a function of reaction progress
*	in the steady frame.  
*/
	
/*ARGSUSED*/
LOCAL	void	MPOLY_progress_state(
	float prog,		/* reaction progress */
	Locstate init,		/* init = state behind front */
	Locstate ans,		/* TGas states, init = state behind front */
	float max_vol)		/* maximum allowed volume of reacted state */
{
	screen("ERROR in MPOLY_progress_state(), nonreactive EOS\n");
	clean_up(ERROR);
}	/* end MPOLY_progress_state*/

/*
*			MPOLY_fprint_combustion_params():
*
*	Prints combustion related parameters.
*/

/*ARGSUSED*/
LOCAL	void	MPOLY_fprint_combustion_params(
	FILE *file,
	Gas_param *params)
{
	(void) fprintf(file,"\tcomposition_type = %d %s\n",
		PURE_NON_REACTIVE,"PURE_NON_REACTIVE");
}
/***************END DETONATION SPECIFIC UTILITY FUNCTIONS*****************/
#endif /* defined(COMBUSTION_CODE) */


/***************METHOD OF CHARACTERISTIC FUNCTIONS FOR W_SPEED**************/

/*
*			MPOLY_neumann_riem_inv_moc():
*
*	Uses the characteristic equations for the k-Riemann invariants
*	to update the states along a Neumann boundary.
*
*	This function integrates the characteristic form of Euler's
*	equations:
*
*	         1     dP       dU             acU
*	       -----  ----  +  ----  =  g  -  -----       (+)
*	       rho c   dl       dl              r
*	                 +        +
*		     
*	         1     dP       dU             acU
*	       -----  ----  -  ----  = -g  -  -----       (-)
*	       rho c   dl       dl              r
*	                 -        -
*
*	               dS
*	              ---- = 0                            (0S)
*	               dl
*	                 0
*		     
*	               dV
*	              ---- = 0                            (0V)
*	               dl
*	                 0
*
*	Here:
*		rho = density
*		P   = pressure
*		S   = specific entropy
*		c   = sound speed
*		U   = component of velocity in the normal direction
*		V   = component of velocity orthogonal to the normal
*		g   = gravitational acceleration
*		a   = geometric coefficient = 0 for rectangular geometry
*		                            = 1 for cylindrical symmetry
*		                            = 2 for spherical symmetry
*
*	Basic geometry:
*
*		side = POSITIVE_SIDE  if the flow is to the right of the wall
*		                      (-) family is directed towards the wall
*
*			/|ans      flow region
*			/|   \
*		    wall/|    \
*			/|     \(-) characteristic
*			/|      \
*			/|_______\________________
*		       pt[0]   pt[0]+dn
*		        st0     st1
*
*		side = NEGATIVE_SIDE if the flow is to the left of the wall
*		                     (+) family is directed towards the wall
*
*			        flow region    ans|\
*			                      /   |\
*			                     /    |\
*		         (+) characteristic /     |\
*			                   /      |\
*			  ________________/_______|\
*				     pt[0]+dn   pt[0]
*		                        st1      st0
*		                flow region
*
*	Basic algorithm:
*	In this function we use equations (0S), (0V), and the slip
*	boundary condition at the wall to compute the updated entropy and
*	velocity at the new state.  Note that we allow a possibly nonzero wall
*	velocity u0.  Thus the entropy at the new state equals the entropy
*	at state st0,  the tangential velocity of the new state equals that
*	of st0,  and the normal component of velocity of the new state
*	equals u0.  The pressure of the new state is found by integrating
*	either equation (+) or (-) depending on the characteristic family
*	directed towards the wall.  The input variable side determines
*	the side of the computational region with respect to the wall.
*
*	NOTE:  some applications may include an artificial heat
*	conduction.  This can be implemented in a variety of ways.
*	One method is to allow an entropy change between st0 and ans
*	that is proportional to the quantity (T1 - T0)/T0  where
*	T0 and T1 are the termperatures of states st0 and st1 respectively.
*
*	Input variables:
*		pt	coordinates of the wall
*		st0	state at the wall at start of time step
*		u0	wall normal velocity
*		c0	sound speed at wall
*		st1	state at position pt + nor*dn at start of time step
*		side	side of flow relative to the wall
*		dn	grid spacing in wall normal direction
*		nor	wall normal
*		front   Front structure
*	Output variables:
*		ans	time updated state at the wall
*		
*/

LOCAL	void	MPOLY_neumann_riem_inv_moc(
	float     *pt,
	Locstate  st0,
	float     u0,
	float     c0,
	Locstate  st1,
	SIDE      side,
	Locstate  ans,
	float     dn,
	float     *nor,
	Front     *front)
{
	RECT_GRID *gr = front->rect_grid;
	float     r0, c0sqr, c1sqr, u1, c1;
	float     v0[3], v1[3];
	float     pr;
	float     time = front->time;
	float     dt = front->dt;
	float     entropy_corr, entropy_coef, entropy_coef1, T1, T0;
	float     vnor;
	float     Sans, S1, dS, dS1;
	float     cans;
	float     sgn;
	float     gam, GAM, Rst;
	float     num, den;
	float     heat_cond;
	float     cf2, cf3, cf6;
	float     a, x;
	float     gans[3], g1[3], g;
	float     ptans[3], pt1[3];
	float     cv;
	int       i, dim;
	float     alpha = nor[0]*rotational_symmetry();

	heat_cond = Params(st0)->avisc.heat_cond;
	dim = gr->dim;
	for (i = 0; i < dim; ++i)
	{
	    v0[i] = vel(i,st0);
	    v1[i] = vel(i,st1);
	    pt1[i] = pt[i]+dn*nor[i];
	    ptans[i] = pt[i] + u0*nor[i]*dt;
	}
	eval_gravity(pt1,time,g1);
	eval_gravity(ptans,time+dt,gans);
	for (g = 0.0, i = 0; i < dim; ++i)
	    g += 0.5*(g1[i]+gans[i])*nor[i];

	u1 = scalar_product(v1,nor,dim);
	sgn = (side == POSITIVE_SIDE) ? -1.0 : 1.0;
	r0 = Dens(st0);
	T0 = temperature(st0);
	c0sqr = sqr(c0);
	T1 = temperature(st1);
	c1sqr = sound_speed_squared(st1);
	c1 = sqrt(c1sqr);
	Set_params(ans,st0);

	/*
	 * Units of heat_cond = units of entropy = 
	 * units of Cv = internal energy / temperature
	 */
	entropy_corr = heat_cond*(T1/T0 - 1.0);
	Sans = entropy(st0) + entropy_corr;
	S1 = entropy(st1);
	gam = Gamma(st0);
	GAM = gam - 1.0;
	cf2 = Coef2(st0);
	cf3 = Coef3(st0);
	cf6 = Coef6(st0);
	Rst = R(st0);
	cv = Rst/GAM;

	dS = dS1 = Sans - S1;
	x = exp(entropy_corr/cv);
	a = x*c0sqr/pow(r0,GAM);
	entropy_coef  = 0.5 * cf3 * dS/Rst;
	entropy_coef1 = 0.5 * cf3 * dS1/Rst;
	num = (1.0 + entropy_coef1)*c1 + sgn*cf2*(u1-u0) + sgn*cf2*g*dt;
	den = 1.0 - entropy_coef;

	if (is_rotational_symmetry() && alpha != 0.0)
	{
	    float rmin, rn, rd;

	    rmin = fabs(pos_radius(0.0,gr));
	    rn = pos_radius(pt1[0],gr);
	    rd = pos_radius(ptans[i],gr);
	    if ((fabs(rn) > rmin) && (fabs(rd) > rmin))
	    {
	        num -= 0.5*alpha*cf2*c1*u1*dt / rn;
	        den += 0.5*alpha*cf2*u0*dt / rd;
	    }
	}
	cans = num/den;
	Dens(ans) = pow(cans*cans/a,cf6);
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            for(i = 0; i < Params(st0)->n_comps; i++)
                pdens(ans)[i] = (pdens(st0)[i]/Dens(st0))*Dens(ans);
        }
	pr = sqr(cans)*Dens(ans)/gam;
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	pr = max(pr,Min_pressure(ans));
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	Press(ans) = pr;
	vnor = scalar_product(nor,v0,dim);
	for (i = 0; i < dim; ++i)
	    Vel(ans)[i] = u0*nor[i] + v0[i] - vnor * nor[i];
	set_type_of_state(ans,TGAS_STATE);
}		/*end MPOLY_neumann_riem_inv_moc*/

/*
*		MPOLY_shock_ahead_state_riem_inv_moc():
*
*	Uses the characteristic equations for the k-Riemann invariants
*	to update the state ahead of a shock wave.
*
*	This function integrates the characteristic form of Euler's
*	equations:
*
*	         1     dP       dU             acU
*	       -----  ----  +  ----  =  g  -  -----       (+)
*	       rho c   dl       dl              r
*	                 +        +
*		     
*	         1     dP       dU             acU
*	       -----  ----  -  ----  = -g  -  -----       (-)
*	       rho c   dl       dl              r
*	                 -        -
*
*	               dS
*	              ---- = 0                            (0S)
*	               dl
*	                 0
*		     
*	               dV
*	              ---- = 0                            (0V)
*	               dl
*	                 0
*
*	Here:
*		rho = density
*		P   = pressure
*		S   = specific entropy
*		c   = sound speed
*		U   = component of velocity in the normal direction
*		V   = component of velocity orthogonal to the normal
*		g   = gravitational acceleration
*		a   = geometric coefficient = 0 for rectangular geometry
*		                            = 1 for cylindrical symmetry
*		                            = 2 for spherical symmetry
*
*	Basic geometry:
*
*		side = POSITIVE_SIDE
*		                    /
*		                   /ans(position = pt0 + W*dt)
*		                  /  + 0 -   
*		                 /    +  0  -   
*		                /      +   0   -   
*		          shock/        +    0    -   
*		         front/          +     0     -   
*		             /            +      0      -   
*			____/st0__________st3____st2____st1____
*		           pt0
*
*		side = NEGATIVE_SIDE
*	                                  \
*			                ans\(position = pt0 + W*dt)
*			             + 0 -  \
*			          +  0  -    \shock
*			       +   0   -      \front
*			    +    0    -        \
*		         +     0     -          \
*                     +      0      -            \
*		____st3____st2____st1__________st0\______
*	                                         pt0
*
*		+ = forward characteristic (velocity = U + c)
*		0 = middle characteristic (velocity = U)
*		- = backward characteristic (velocity = U - c)
*
*	Basic algorithm:
*	The entropy and tangential component of velocity of the state
*	ans are found by integrating the (0S) and (0V) characteristic
*	equations.  State st2 is the state at the foot of these
*	characteristics,  so the entropy and tangential velocity of
*	ans is the same as that of st2.  The pressure and normal
*	component of velocity of ans are then obtained by integrating
*	the characteristic equations (+) and (-).
*
*	Input variables:
*		st0	state at foot of shock at start of time step
*		st1	state at foot of (-) or (U-c) charateristic 
*		st2	state at foot of (0) or (U) characteristic
*		st3	state at foot of (+) or (U+c) characteristic
*		pt0	coordinates of shock at start of time step
*		side	ahead side of shock
*		add_source	include gravitational and geometric sources
*		dn	grid spacing in wall normal direction
*	        f1      location of st1 = pt0 + f1*dn*nor
*	        f2      location of st2 = pt0 + f2*dn*nor
*	        f3      location of st3 = pt0 + f3*dn*nor
*		nor	wall normal
*		W	predicted shock front velocity
*		front   Front structure
*	Ouput variables:
*		ans	time updated state ahead of the shock
*
*/

/*ARGSUSED*/
LOCAL	void	MPOLY_shock_ahead_state_riem_inv_moc(
	float     *pt0,
	Locstate  st0,
	Locstate  st1,
	Locstate  st2,
	Locstate  st3,
	Locstate  ans,
	float     dn,
	float     f1,
	float     f2,
	float     f3,
	float     *nor,
	float     *W,
	int       add_source,
	float     dt,
	Front     *front)
{
	RECT_GRID *gr = front->rect_grid;
	float     time = front->time;
	float     u1, u2, u3;
	float     S1, S2, S3;
	float     u;
	float     c, c1, c3, c1sqr, c2sqr, c3sqr;
	float     pr, p2, r2;
	float     alpha1, alpha3;
	float     gam, GAM, cv, cf6;
	float     ptans[3], pt1[3], pt3[3];
	float     gans[3], g1[3], g3[3], g1_bar, g3_bar;
	float     A1, A3;
	float     ua1, ua3;
	int       i, dim;
	float     radans = 0.0;
	float     alpha = nor[0]*rotational_symmetry();

	dim = gr->dim;
	u1 = 0.0;	u2 = 0.0;	u3 = 0.0;
	for (i = 0; i < dim; ++i)
	{
	    u1 += nor[i] * Vel(st1)[i];
	    u2 += nor[i] * Vel(st2)[i];
	    u3 += nor[i] * Vel(st3)[i];
	    pt1[i] = pt0[i] + f1*dn*nor[i];
	    pt3[i] = pt0[i] + f3*dn*nor[i];
	    ptans[i] = pt0[i] + W[i]*dt;
	}
	if (add_source)
	{
	    eval_gravity(pt1,time,g1);
	    eval_gravity(pt3,time,g3);
	    eval_gravity(ptans,time+dt,gans);
	    for (g1_bar = 0.0, g3_bar = 0.0, i = 0; i < dim; ++i)
	    {
	        g1_bar += 0.5*(gans[i] + g1[i])*nor[i];
	        g3_bar += 0.5*(gans[i] + g3[i])*nor[i];
	    }
	}
	else
	{
	    g1_bar = 0.0;
	    g3_bar = 0.0;
	}

	c1sqr = sound_speed_squared(st1);
	S1 = entropy(st1);
	c1 = sqrt(c1sqr);

	r2 = Dens(st2);
	p2 = pressure(st2);
	c2sqr = sound_speed_squared(st2);
	S2 = entropy(st2);

	c3sqr = sound_speed_squared(st3);
	S3 = entropy(st3);
	c3 = sqrt(c3sqr);

	Set_params(ans,st0);
	gam = Gamma(st0);
	GAM = gam - 1.0;
	cv = R(st0)/GAM;
	cf6 = Coef6(st0);

	alpha1 = (S2 - S1)/(4.0*gam*cv);
	alpha3 = (S2 - S3)/(4.0*gam*cv);
	A1 = 0.5*GAM*(u1 + g1_bar*dt) - c1*(1.0 + alpha1);
	A3 = 0.5*GAM*(u3 + g3_bar*dt) + c3*(1.0 + alpha3);

	if (is_rotational_symmetry() && add_source && alpha != 0.0)
	{
	    float rmin = fabs(pos_radius(0.0,gr));
	    float ra = pos_radius(ptans[0],gr);
	    float r1 = pos_radius(pt1[0],gr);
	    float r3 = pos_radius(pt3[0],gr);

	    radans = (fabs(ra) > rmin) ? 0.5*alpha*dt/ra : 0.0;

	    if (fabs(r1) > rmin)
	        A1 += 0.25*GAM*alpha*c1*u1*dt/r1;
	    if (fabs(r3) > rmin)
	        A3 -= 0.25*GAM*alpha*c3*u3*dt/r3;
	}

	if (is_rotational_symmetry() && radans != 0.0)
	{
	    float A, B, C, D;

	    A = radans*(alpha1 - alpha3);
	    B = 2.0 - (alpha1 + alpha3) +  radans*(A1 + A3);
	    C = A3 - A1;
	    D = 1.0 - 4.0*A*C/(B*B);
	    if (D < 0.0)
	    {
		screen("ERROR in MPOLY_shock_ahead_state_riem_inv_moc(), "
		       "negative discriminate\n");
		clean_up(ERROR);
	    }
	    c = (C/B)/(0.5*(1.0 + sqrt(D)));
	}
	else
	    c = (A3 - A1)/(2.0 - (alpha1 + alpha3));

	Dens(ans) = r2 * pow(c*c/c2sqr,cf6);
            /* Mass fraction is obtained from st2, scaling is performed */
            /* Shall we use char. to update partial density ?? */
            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
                int    k;
                float  sum; 
                if(Params(st0)->n_comps != 1 && Params(st2)->n_comps != 1)
                {
                    if(Params(st0)->n_comps != Params(st2)->n_comps)
                    {
                        (void)printf("ERROR, MPOLY_shock_ahead_state_riem_inv_moc,"
                          " different # of components\n");
                        clean_up(ERROR); 
                    }
                    sum = 0.0;
                    for(k = 0; k < Params(st2)->n_comps; k++)
                    {
                        pdens(ans)[k] = (pdens(st2)[k]/Dens(st2))*Dens(ans);
                        sum += pdens(ans)[k];
                    }
                    for(k = 0; k < Params(st2)->n_comps; k++)
                        pdens(ans)[k] *= Dens(ans)/sum;
                }
            }
	pr = p2 * pow(Dens(ans)/r2,gam);
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	pr = max(pr,Min_pressure(ans));
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	Press(ans) = pr;
	set_type_of_state(ans,TGAS_STATE);
	ua1 = (2.0/GAM)*(A1 + c*(1.0 - alpha1));
	ua3 = (2.0/GAM)*(A3 - c*(1.0 - alpha3));
	if (is_rotational_symmetry() && radans != 0.0)
	{
	    ua1 /= 1.0 - radans*c;
	    ua3 /= 1.0 + radans*c;
	}
	u = 0.5*(ua1 + ua3);
	for (i = 0; i < dim; ++i)
	    Vel(ans)[i] = Vel(st2)[i] + (u - u2)*nor[i];
}		/*end MPOLY_shock_ahead_state_riem_inv_moc*/


/*
*			MPOLY_shock_moc_plus_rh():
*
*	Given the updated state st0 on the ahead shock, and the state
*	st1 at the foot of the characteristic through the behind state
*	this function uses the method of characteristics and the
*	Rankine-Hugoniot conditions to solve for the updated state
*	behind the shock.
*
*	This function integrates the characteristic form of Euler's
*	equations:
*
*	         1     dP       dU             acU
*	       -----  ----  +  ----  =  g  -  -----       (+)
*	       rho c   dl       dl              r
*	                 +        +
*		     
*	         1     dP       dU             acU
*	       -----  ----  -  ----  = -g  -  -----       (-)
*	       rho c   dl       dl              r
*	                 -        -
*
*	               dS
*	              ---- = 0                            (0S)
*	               dl
*	                 0
*		     
*	               dV
*	              ---- = 0                            (0V)
*	               dl
*	                 0
*
*	together with the Hugoniot conditions across a shock
*
*		      rho  *(U  - s) = rho  * (U  - s) = m
*                        0    0           1     1
*
*                             2                       2
*		rho  *(U  - s)  + P  = rho1 * (U  - s)  + P1
*		   0    0          0      1     1
*
*	                                1
*		             E  - E  = --- (P + P ) * (V  - V )
*		              1    0    2    1   0      0    1 
*
*				  V  = V
*	                           0    1
*
*			
*
*	Here:
*		rho = density
*		P   = pressure
*		S   = specific entropy
*		c   = sound speed
*		s   = shock normal velocity
*		U   = component of velocity in the normal direction
*		V   = component of velocity orthogonal to the normal
*		g   = gravitational acceleration
*		a   = geometric coefficient = 0 for rectangular geometry
*		                            = 1 for cylindrical symmetry
*		                            = 2 for spherical symmetry
*
*		The subscripts on the Hugoniot equations represent
*		the ahead (0) and behind (1) shock states respectively.
*
*	Basic geometry:
*
*		side = POSITIVE_SIDE
*		                    /
*		                ans/sta(position = pt + W*dt)
*		                + /  + 0 -   
*		             +   /    +  0  -   
*		          +     /      +   0   -   
*		       +       /shock   +    0    -   
*		    +         /front     +     0     -   
*	         +           /            +      0      -   
*	    __stb________stm/______________+_______0_______-___
*		           pt
*
*		side = NEGATIVE_SIDE
*	                                  \
*			                sta\ans(position = pt + W*dt)
*			             + 0 -  \ -
*			          +  0  -    \    -    
*			       +   0   -      \       - 
*			    +    0    -   shock\         -
*		         +     0     -     front\            -
*                     +      0      -            \              -
*		___+_______0_______-______________\stm___________stb__
*	                                          pt
*
*		+ = forward characteristic (velocity = U + c)
*		0 = middle characteristic (velocity = U)
*		- = backward characteristic (velocity = U - c)
*
*	Basic algorithm:
*	The state behind the shock and the shock velocity is determined from
*	the state ahead of the shock and one additional piece of information
*	which is obtained by integrating the incoming behind shock
*	characteristic.  Basically the discretized integration of the
*	behind shock incoming characteristic and the Rankine-Hugoniot 
*	equations across the shock provide a complete set of equations
*	to determine the time updated state behind the shock,  and the
*	time updated shock velocity.  It is common in practice to take
*	the net shock velocity and the average of the wave velocity
*	computed from the Riemann solution at the start of the time step
*	and the valued computed from the above set of equations for the
*	velocity at the end of the time step.
*
*	Input variables:
*		sta	state ahead of time updated shock
*		stm	state at foot of shock at start of time step
*		stb	state behind shock at foot of incoming chacteristic
*		pt	coordinates of shock at start of time step
*		dn	grid spacing in wall normal direction
*		nor	wall normal
*		W	first prediction of shock front velocity
*		front   Front structure
*	Output variables:
*		ans	times updated state behind shock
*		W	updated shock speed
*
*		The answer state ans is returned in Gas format.
*/


LOCAL	bool	MPOLY_shock_moc_plus_rh(
	float     *pt,
	Locstate  sta,
	Locstate  stm,
	Locstate  stb,
	Locstate  ans,
	float     dn,
	float     *nor,
	float     *W,
	int       w_type,
	Front     *front)
{
	RECT_GRID    *gr = front->rect_grid;
	float        time = front->time;
	float        dt = front->dt;
	float        g, ga[3], gb[3];
	float	     hld_va[MAXD], hld_vb[MAXD];
	float	     ua, ub;
	float	     Ws, u_ans;
	float	     nr = nor[0];
	float        ca;
	float        pa, ra, pb, pm, rb, xa, pi;
	float        xl, xr;
	float        sgn;
	float        p, r, u;
	float        f0;
	float        epsilon, delta;
	const float  eps = 10.0*MACH_EPS;/*TOLERANCE*/
	float        c4, gam, eta, psi;
	float        pta[3], ptb[3];
	MPOLY_S_MPRH Smprh;
	float        d;
	float        cf2, cf3;
	bool         status = FUNCTION_SUCCEEDED;
	int	     i, dim;
	static const int    mnth = 100; /*TOLERANCE*/
	float alpha = rotational_symmetry();

	sgn = (is_forward_wave(w_type)) ? 1.0 : -1.0;
	dim = gr->dim;

	if (is_rarefaction_wave(w_type))
	{
	    set_state(ans,GAS_STATE,sta);
	    for (ua = 0.0, i = 0; i < dim; ++i)
	    	ua += nor[i] * Vel(sta)[i];
	    ca = sound_speed(sta);
	    Ws = ua + sgn*ca;
	    for (i = 0; i < dim; ++i)
	    	W[i] = 0.5 * (W[i] + Ws * nor[i]);
	    return status;
	}

	ua = ub = 0.0;
	for (i = 0; i < dim; ++i)
	{
	    hld_va[i] = Vel(sta)[i];
	    hld_vb[i] = Vel(stb)[i];
	    ua += nor[i] * Vel(sta)[i];
	    ub += nor[i] * Vel(stb)[i];
	    pta[i] = pt[i] + W[i]*dt;
	    ptb[i] = pt[i] - dn*nor[i];
	}
	eval_gravity(pta,time + dt,ga);
	eval_gravity(ptb,time,gb);
	g = 0.5*(scalar_product(ga,nor,dim) + scalar_product(gb,nor,dim));

	Vel(sta)[0] = ua;
	Vel(stb)[0] = ub;


	ra = Dens(sta);		pa = Press(sta);	ua = Vel(sta)[0];
	rb = Dens(stb);		pb = Press(stb);	ub = Vel(stb)[0];
	pm = Press(stm);

	gam = Gamma(sta);
	cf2 = Coef2(sta);
	cf3 = Coef3(sta);
	c4 = Coef4(sta);

	Smprh.beta = 0.25/gam;
	Smprh.gam = gam;
	Smprh.c2 = cf2;
	Smprh.c4 = c4;
	Smprh.ca = ca = sound_speed(sta);
	Smprh.cb = sound_speed(stb);
	Smprh.psi = psi = cf3;
	Smprh.eta = eta = 1.0 - psi;	 /* (gam + 1)/(2*gam) */
	Smprh.alpha = ca*psi;
	Smprh.dS = Smprh.beta * (log(pa/pb) - gam*log(ra/rb));

	f0 = ub - ua;
	if (is_gravity() == YES)
	    f0 += g*dt;
	f0 *= sgn;

	if (is_rotational_symmetry() && alpha > 0.0)
	{
	    float rmin, rada, radb;

	    Smprh.ua = ua;
	    Smprh.delta = sgn/cf2;
	    rmin = fabs(pos_radius(0.0,gr));

	    Ws = scalar_product(W,nor,dim);
	    rada = pos_radius(pta[0],gr);
	    radb = pos_radius(ptb[0],gr);
	    if ((fabs(rada) > rmin) && (fabs(radb) > rmin))
	    {
		float cb;
	        Smprh.cyl_coord = YES;
	        cb = sound_speed(stb);
	        f0 -= 0.5*nr*cb*ub*dt*alpha/radb;
	        Smprh.rad = 0.25*nr*cf2*dt*alpha/rada;
	    }
	    else
	    {
	        Smprh.cyl_coord = NO;
		Smprh.rad = 0.0;
	    }
	}
	else
	    Smprh.cyl_coord = NO;

	xa = pm/pa;
	xa = max(1.0,xa);

	f0 *= cf2;

	delta = max(xa*EPS, eps);
	epsilon = max(f0*EPS, eps);
	xl = max(1.0 - delta,0.5*xa);/*TOLERANCE*/
	xr = 1.5*xa;/*TOLERANCE*/

	if ((find_root(mpoly_fS,(POINTER) &Smprh,f0,&pi,xl,xr,
				 epsilon, delta) == FUNCTION_FAILED)
	    &&
	    (search_harder_for_root(mpoly_fS,(POINTER) &Smprh,f0,&pi,
				    xl,xr,&xl,&xr,1.0,HUGE_VAL,
				    mnth,epsilon,delta) == FUNCTION_FAILED))
	{
	    float rpi;
	    /* Try looking for root < 1 */
	    if (find_root(mpoly_fS,(POINTER) &Smprh,f0,&rpi,0.0,1.0,
			  epsilon, delta) == FUNCTION_FAILED)
	    {
	        (void) printf("WARNING in MPOLY_shock_moc_plus_rh(), "
	                      "for %s shock solution\n",
	                      (is_backward_wave(w_type)) ?
			      "backward" : "forward");
	        status = FUNCTION_FAILED;
	    }
	    else
	    {
	        float fpi, frpi;
	        (void) printf("WARNING in MPOLY_shock_moc_plus_rh(), "
                              "shock reduced to zero strength\n");
                if (mpoly_fS(pi,&fpi,(POINTER)&Smprh) &&
                    mpoly_fS(rpi,&frpi,(POINTER)&Smprh))
                {
                    if (fabs(frpi - f0) < fabs(fpi - f0))
                        pi = 1.0;
                }
                else
		{
		    status = FUNCTION_FAILED;
	            pi = 1.0;
		}
	    }
	}
	if (pi <= (1.0 - delta)) 
	{
	    (void) printf("WARNING in MPOLY_shock_moc_plus_rh(), "
	                  "Non-physical solution for R-H conditions\n");
	    status = FUNCTION_FAILED;
	}
	if (pi < 1.0)
	    pi = 1.0;
	p = pa * pi;
	r = ra * (pi + c4)/(1.0 + c4*pi);
	d = sgn * ca / sqrt(psi + eta*pi);
	u = ua + d * (pi - 1.0) / gam;

	Dens(ans) = r;
        /* Mass fraction of stb is preserved, scaling is performed */
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            int    k;
            float  sum;
            if(Params(stb)->n_comps != 1)
            {
                sum = 0.0;
                for(k = 0; k < Params(stb)->n_comps; k++)
                    sum += pdens(stb)[k];
                for(k = 0; k < Params(stb)->n_comps; k++)
                    pdens(ans)[k] = Dens(ans)*pdens(stb)[k]/sum;
            }
        }
	Vel(ans)[0] = u;
	Press(ans) = p;
	Set_params(ans,sta);
	set_type_of_state(ans,TGAS_STATE);
	Ws = ua + d * (pi + c4) / (1.0 + c4);

	for (i = 0; i < dim; ++i)
	{
	    Vel(sta)[i] = hld_va[i];
	    Vel(stb)[i] = hld_vb[i];
	}

	u_ans = Vel(ans)[0];
	for (i = 0; i < dim; ++i)
	    Vel(ans)[i] = Vel(sta)[i] + (u_ans - ua) * nor[i];
	set_type_of_state(ans,TGAS_STATE);

	set_state(ans,GAS_STATE,ans);

	for (i = 0; i < dim; ++i)
	    W[i] = 0.5 * (W[i] + Ws * nor[i]);
	return status;
}		/*end MPOLY_shock_moc_plus_rh*/


/***************END METHOD OF CHARACTERISTIC FUNCTIONS FOR W_SPEED**********/

/***************INITIALIZATION UTILITY FUNCTIONS****************************/

/*
*			MPOLY_prompt_for_state():
*
*	Prompts for a hydrodynamical state.  The form of
*	the states depends of the Eos. 	The type of the state
*	is returned.
*/

LOCAL	void	MPOLY_prompt_for_state(
	Locstate   state,
	int        stype,
	Gas_param  *params,
	const char *mesg)
{
	int i, dim;
	static  char velmesg[3][11] = {"x velocity","y velocity","z velocity"};

	if (params == NULL)
	{
		g_obstacle_state(state,g_sizest());
		return;
	}
	dim = params->dim;
	set_type_of_state(state,TGAS_STATE);
	Params(state) = params;
	screen("Enter the density, pressure");
	for (i = 0; i < dim; ++i)
	{
		screen(", ");
		if (i == (dim - 1)) screen("and ");
		screen("%s",velmesg[i]);
	}
	screen("%s: ",mesg);
	(void) Scanf("%f %f",&Dens(state),&Press(state));
	for (i = 0; i < dim; ++i)
		(void) Scanf("%f",&Vel(state)[i]);
	(void) getc(stdin); /*read trailing newline*/

	set_state(state,stype,state);
}		/* end MPOLY_prompt_for_state */

/*
*			MPOLY_prompt_for_thermodynamics():
*
*	Prompts for the thermodynamic variables in a state.  Returns
*	a state with the appropriate thermodynamic state and zero velocity.
*	The return status gives the state type representation of the state.
*/

LOCAL	void	MPOLY_prompt_for_thermodynamics(
	Locstate   state,
	Gas_param  *params,
	const char *mesg)
{
	if (params == NULL)
	{
		g_obstacle_state(state,g_sizest());
		return;
	}
	set_type_of_state(state,TGAS_STATE);
	zero_state_velocity(state,MAXD);
	Params(state) = params;
	screen("Enter the density and pressure");
	screen("%s: ",mesg);
	(void) Scanf("%f %f\n",&Dens(state),&Press(state));
}		/* end MPOLY_prompt_for_thermodynamics */

/*
*			MPOLY_fprint_EOS_params():
*
*	Prints the parameters that define the given equation of state.
*	NOTE:  This is not really an initialization function,  but it is
*	convenient to locate it next the the corresponding read function.
*/

LOCAL	void	MPOLY_fprint_EOS_params(
	FILE *file,
	Gas_param *params)
{
	float	*gam = ((MPOLY_EOS *) params->eos)->_gamma;
	float	*M = ((MPOLY_EOS *) params->eos)->_M;
	float	R = ((MPOLY_EOS *) params->eos)->R;
	int	i, nc = params->n_comps;

	(void) fprintf(file,"\tEquation of state = %d MULTI_COMP_POLYTROPIC\n",
		MULTI_COMP_POLYTROPIC);
	(void) fprintf(file,"\tnumber of components = %d\n",params->n_comps);
	(void) fprintf(file,"\tgamma = ");
	if (is_binary_output() == YES)
	{
	    (void) fprintf(file,"\f%c",nc);
	    (void) fwrite((const void *) gam,FLOAT,nc,file);
	}
	else
	{
	    for (i = 0; i < nc; ++i)	
	    	(void) fprintf(file,"%"FFMT"%s",gam[i],(i == (nc-1))?"\n":" ");
	}
	(void) fprintf(file,"\tM = ");
	if (is_binary_output() == YES)
	{
	    (void) fprintf(file,"\f%c",nc);
	    (void) fwrite((const void *) M,FLOAT,nc,file);
	}
	else
	{
	    for (i = 0; i < nc; ++i)	
	    	(void) fprintf(file,"%"FFMT"%s",M[i],(i == (nc-1))?"\n":" ");
	}
	(void) fprintf(file,", R = ");
	if (is_binary_output() == YES)
	{
	    (void) fprintf(file,"\f%c",1);
	    (void) fwrite((const void *) &R,FLOAT,1,file);
	}
	else
	    (void) fprintf(file,"%"FFMT,R);
	(void) fprintf(file,"\n");
}		/*end MPOLY_fprint_EOS_params */

/*
*			MPOLY_read_print_EOS_params():
*
*	Reads the equation of state data as printed by MPOLY_fprint_EOS_params.
*	This is restart function.
*/

/*ARGSUSED*/
LOCAL	void	MPOLY_read_print_EOS_params(
	INIT_DATA     *init,
	const IO_TYPE *io_type,
	Gas_param     *params)
{
	FILE      *file = io_type->file;
	MPOLY_EOS *mpeos = (MPOLY_EOS *)params->eos;
	float	  *gam =  mpeos->_gamma;
	float	  *M = mpeos->_M;
	int	  c, i, nc;

	(void) fgetstring(file,"number of components = ");
	(void) fscanf(file,"%d",&nc);
	params->n_comps = nc;
	(void) fgetstring(file,"gamma = ");
	if ((c = getc(file)) != '\f') /*NOBINARY*/
	{
	    (void) ungetc(c,file);
	    for (i = 0; i < nc; ++i)
	    	(void) fscan_float(file,gam+i);
	}
	else
	{
	    (void) getc(file);
	    (void) read_binary_real_array(gam,nc,io_type);
	}
	(void) fgetstring(file,"M = ");
	if ((c = getc(file)) != '\f') /*NOBINARY*/
	{
	    (void) ungetc(c,file);
	    for (i = 0; i < nc; ++i)
	    	(void) fscan_float(file,M+i);

	}
	else
	{
	    (void) getc(file);
	    (void) read_binary_real_array(M,nc,io_type);
	}
	mpeos->R = fread_float("R = ",io_type);
}		/*end MPOLY_read_print_EOS_params*/

/*
*			MPOLY_free_EOS_params():
*
*	Frees the storage allocated for an equation of state parameter
*	function.
*/

LOCAL	EOS*	MPOLY_free_EOS_params(
	EOS *eos)
{
	free(eos);
	return NULL;
}		/*end MPOLY_free_EOS_params*/

/*
*			MPOLY_prompt_for_EOS_params():
*
*	Prompts for equation of state parameters.
*/

/*ARGSUSED*/
LOCAL	void	MPOLY_prompt_for_EOS_params(
	INIT_DATA  *init,
	Gas_param  *params,
	const char *message1,
	const char *message2)
{
	MPOLY_EOS	*mpeos = (MPOLY_EOS*) params->eos;
	float		*gam =  mpeos->_gamma;
	float		*M = mpeos->_M;
	char		s[120];
	int		n, nc = 2;
	size_t		len1, len2, len;


	len1 = (message1 == NULL) ? 0 : strlen(message1);
	len2 = (message2 == NULL) ? 0 : strlen(message2);
	len = len1 + len2;

	screen("Enter the number of components (< %d)",MAX_NUM_GAS_COMPS);
	screen("%s",(len > 25) ? "\n\t" : " ");
	screen("for the%s gas%s (default = %d): ",message1,message2,nc);
	(void) Gets(s);
	if (s[0] != '\0')
		(void) sscanf(s,"%d\n",&nc);
	if (nc <= 0 || nc > MAX_NUM_GAS_COMPS)
	{
	    screen("ERROR in MPOLY_prompt_for_EOS_params(), "
	           "invalid number of gas components %d\n"
	           "Value must be in the range 1 <= n <= %d\n",
		   nc,MAX_NUM_GAS_COMPS);
	    clean_up(ERROR);
	}
        else
            params->n_comps = nc;
	for (n = 0; n < nc; ++n)
	{
	    screen("Enter the %d%s component's "
	           "polytropic exponent (gamma) %s"
	           "for the%s gas%s: ",n+1,ordinal_suffix(n+1),
		   (len > 25) ? "\n\t" : " ",message1,message2);
	    (void) Scanf("%f\n",gam+n);
	    screen("Enter the %d%s component's "
	           "molecular weight (M)%s for the%s gas%s: ",
		   n+1,ordinal_suffix(n+1),(len > 25) ? "\n\t" : " ",
		   message1,message2);
	    (void) Scanf("%f\n",M+n);
	}
	mpeos->R = 1;
	screen("Enter the ideal gas constant (R, PV = (R/m)T, ");
	screen("(m = molecular weight),\ndefault for R = %g)\n",mpeos->R);
	screen("\tfor the%s gas%s: ",message1,message2);
	(void) Gets(s);
	if (s[0] != '\0')
		(void) sscan_float(s,&mpeos->R);
}		/*end MPOLY_prompt_for_EOS_params*/




/***************Problem Type Specific Initialization Functions**************/

/*
*			MPOLY_RT_RS_f():
*
*	Support function for the computation of a solution to the linearized
*	Rayleigh-Taylor flow.
*
*	NEEDED:  More complete documentation
*/

LOCAL	float	MPOLY_RT_RS_f(
	float		s,
	Locstate	amb_st,
	float		dz,
	float		k_sqr,
	float		g_z)
{
	float gam;
	float h_sqr;
	float csqr = sound_speed_squared(amb_st);
	float rho = Dens(amb_st);
	float alpha1, alpha2;
	float D, D1, D2, N;
	float arg1, arg2;
	float beta;

	gam = Gamma(amb_st);
	beta = gam * g_z / csqr;
	h_sqr = 0.25*sqr(beta) + s/csqr + k_sqr
			+ (gam-1.0)*sqr(g_z)*k_sqr/(s*csqr);
	if (h_sqr < 0.0)
	{
	    screen("ERROR in RT_RS_f(), h_sqr = %g < 0\n",h_sqr);
	    screen("s = %g, beta = %g, csqr = %g\n",s,beta,csqr);
	    screen("k_sqr = %g, gam = %g, g_z = %g\n",k_sqr,gam,g_z);
	    clean_up(ERROR);
	}
	alpha1 = -sqrt(h_sqr) - 0.5*beta;
	alpha2 =  sqrt(h_sqr) - 0.5*beta;
	if (alpha1*dz >= alpha2*dz)
	{
	    arg1 = alpha1*dz;
	    arg2 = alpha2*dz;
	    D = 1.0 - exp(arg2 - arg1);
	}
	else
	{
	    arg1 = alpha2*dz;
	    arg2 = alpha1*dz;
	    D = exp(arg2 - arg1) - 1.0;
	}
	N = (gam-1.0)*sqr(g_z) + s*csqr;
	D1 = (gam-1.0)*g_z + alpha1*csqr;
	D2 = (gam-1.0)*g_z + alpha2*csqr;
	return  rho*N*(exp(alpha1*dz-arg1)/D1 - exp(alpha2*dz-arg1)/D2)/D - g_z*rho;
}		/*end MPOLY_RT_RS_f*/


/*
*			MPOLY_RT_single_mode_perturbation_state():
*
*	Computes the perturbation term for the solution to the linearized
*	Euler equations in the Rayleigh-Taylor problem.  See the appendix to
*
*	``The Dynamics of Bubble Growth for
*				Rayleigh-Taylor Unstable Interfaces''
*	C. L. Gardner, J. Glimm, O. McBryan, R. Menikoff, D. H. Sharp,
*	and Q. Zhang, Phys. Fluids 31 (3), 447-465 (1988).
*
*	for an explanation of the formulas.
*
*       To avoid overflows in some exponential terms, the solution is divided
*       into the light and heavy fluid cases, as determined by the sign of
*       zh = z_bdry - z_intfc.  For the case where zh > 0 (light fluid), the
*       equation for the pressure perturbation dP is multiplied by
*       exp(alpham*zh)/exp(alpham*zh), and for the case where zh < 0 (heavy
*       fluid), it's multiplied by exp(alphap*zh)/exp(alphap*zh).  All other
*       formulas are exactly as they appear in Gardner et al, with the
*       exception of a sign difference in the density perturbation which
*       corrects a typo in the paper.
*
*       Note that ans is only the perturbation to the ambient state.
*/

/*ARGSUSED*/
LOCAL	void	MPOLY_RT_single_mode_perturbation_state(
	Locstate	ans,
	float		*coords,
	float		t,
	Locstate	amb_st,
	float		z_intfc,
	float		z_bdry,
	MODE		*mode,
	float		g_z)
{
	int 		j, dim = Params(amb_st)->dim;

	float 		gam = Gamma(amb_st);
	float 		csqr = sound_speed_squared(amb_st);
	float 		beta = gam * g_z / csqr;

	float 		A = mode->amplitude;
	float 		sigma = mode->growth_rate;
	float 		timefac = exp(sigma*t);
	float 		*k = mode->wave_number;
	float 		phi = scalar_product(k,coords,dim-1) - mode->phase;

	float		h, alphap, alpham;
	float		z, zh, N, D, Dp, Dm;
	float           prefac, expp, expm;
	float           dP, dPdz, drho, dv_z;

	h = 0.25*sqr(beta) + sqr(sigma)/csqr +
	  scalar_product(k,k,dim-1)*(1.0+(gam-1.0)*sqr(g_z)/(sqr(sigma)*csqr));
	h = sqrt(h);

	alphap = -0.5*beta + h;
	alpham = -0.5*beta - h;
	
	z = coords[dim-1];
	zh = z_bdry - z_intfc;
	N = (gam-1.0) * sqr(g_z) + sqr(sigma) * csqr;
	Dp = (gam-1.0) * g_z + alphap * csqr;
	Dm = (gam-1.0) * g_z + alpham * csqr;

	/* To avoid overflows in the exponential terms, we have to
	   separate the light (above the interface) and heavy (below
	   the interface) fluid cases, as determined by the sign of
	   zh. */
	
	if (zh > 0.0)      /* light fluid */
	{
	    D = 1.0 - exp(-2.0*h*zh);
	    expm = exp((alpham+beta) * (z-z_intfc)) / Dm;
	    expp = exp((alpham+beta) * zh
			   + (alphap+beta) * (z-z_bdry)) / Dp;
	}
	else               /* heavy fluid */
	{
	    D = exp(2.0*h*zh) - 1.0;
	    expm = exp((alpham+beta) * (z-z_bdry)
			   + (alphap+beta) * zh) / Dm;
	    expp = exp((alphap+beta) * (z-z_intfc)) / Dp;
	}
	prefac = N * Dens(amb_st) * A / D;

	dP = -prefac * (expm - expp);
	dPdz = -prefac * ((alpham+beta) * expm - (alphap+beta) * expp);
	drho = (sqr(sigma) * dP + (gam-1.0) * g_z * dPdz) / N;
	dv_z = (g_z * drho - dPdz) / (sigma * Dens(amb_st));

	set_type_of_state(ans,TGAS_STATE);
	Set_params(ans,amb_st);

	Press(ans) = dP*timefac*sin(phi);
	Dens(ans) = drho*timefac*sin(phi);
	Vel(ans)[dim-1] = dv_z*timefac*sin(phi);

	for (j = 0; j < dim-1; ++j)
	    Vel(ans)[j] = -k[j]*dP*cos(phi)/(sigma*Dens(amb_st));

	if (debugging("RT_RS_all")) 
	{
	    (void) printf("\nIn POLY_RT_single_mode_perturbation_state(),\n");
	    (void) printf("A = %g; sigma = %g; timefac = %g; "
			  "k[0] = %g; xfac = %g\n",
			  A,sigma,timefac,k[0],sin(phi));
	    (void) printf("gam = %g; csqr = %g; beta = %g\n",gam,csqr,beta);
	    (void) printf("h = %g; alphap = %g; alpham = %g\n",h,alphap,alpham);
	    (void) printf("zh = %g; N = %g; Dp = %g; Dm = %g\n",zh,N,Dp,Dm);
	    (void) printf("D = %g; expm = %g; expp = %g; prefac = %g\n",
			   D,expm,expp,prefac);
	    (void) printf("dP = %g; dPdz = %g; drho = %g; dv_z = %g\n",
			  dP,dPdz,drho,dv_z);
	    (void) printf("x = %g; z = %g; t = %g; z_intfc = %g; z_bdry = %g\n",
		          coords[0],coords[1],t,z_intfc,z_bdry);
	    (void) printf("P = %g; d = %g; Vx = %g; Vz = %g\n",
		       Press(amb_st)+Press(ans),Dens(amb_st)+Dens(ans),
		       Vel(amb_st)[0]+Vel(ans)[0],
		       Vel(amb_st)[dim-1]+Vel(ans)[dim-1]);
	}	
}		/*end MPOLY_RT_single_mode_perturbation_state*/



/*
*		KH_single_mode_state():
*
*	Computes the state at location coords and time t for the solution of
*	the linearized Euler equations for a single mode Kelvin-Helmholtz
*	perturbation.
*	
*	See I. G. Currie, Fundamental Mechanics of Fluids, Chapter 6.,
*	or see Lamb's Hydrodynamics for the incompressible analysis;
*	for the compressible analysis, use Crocco's equation in place
*	of Bernoulli's equation and the wave equation in place of
*	Laplace's equation.
*/

LOCAL	void	MPOLY_KH_single_mode_state(
	Locstate	ans,
	float		*coords,
	float		t,
	Locstate	amb_st,
	float		stream_velocity,
	float		z_intfc,
	float		z_bdry,
	MODE		*mode)
{
	float gam;
	float amb_press = pressure(amb_st);

	Set_params(ans,amb_st);
	set_type_of_state(ans,TGAS_STATE);
	gam = Gamma(amb_st);
			/*TODO 3D*/
	FORTRAN_NAME(khstate)(coords,coords+1,&t,
		&Dens(amb_st),&z_intfc,
		&z_bdry,&mode->amplitude,mode->wave_number,
		&mode->growth_rate,&mode->propagation_speed,
		&mode->phase,Vel(ans),Vel(ans)+1,
		&Dens(ans),&Press(ans),&amb_press,
		&stream_velocity,&gam);
}		/*end MPOLY_KH_single_mode_state*/


/*
*		MPOLY_compute_isothermal_stratified_state():
*
*	Solves for the state at height dz above the reference state
*	ref_st in an isothermal one dimensional steady flow.
*
*	The solution is computed by solving the differential
*	equation:
*
*		P_z = rho gz,	P(0) = P_R, rho(0) = rho_r, T = T_r.
*/

LOCAL	void	MPOLY_compute_isothermal_stratified_state(
	Locstate	ans,
	float		dz,	/* distance from reference position */
	float		gz,	/* gravity */
	Locstate	ref_st)
{
	float gam, tmp;
	float rho, pr;
	int dim = Params(ref_st)->dim;

	gam = Gamma(ref_st);
	Set_params(ans,ref_st);
	tmp = exp((gam*gz/sound_speed_squared(ref_st))*dz);
	rho = Dens(ref_st);
	pr = pressure(ref_st);
	Dens(ans) = rho*tmp;
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
	    int i;
            for(i = 0; i < Params(ref_st)->n_comps; i++)
            {
                pdens(ans)[i] = (pdens(ref_st)[i]/rho)*Dens(ans);
            }
        }
	Press(ans) = pr*tmp;
	set_type_of_state(ans,TGAS_STATE);
	zero_state_velocity(ans,dim);
}		/*end MPOLY_compute_isothermal_stratified_state */


/*
*		MPOLY_compute_isentropic_stratified_state():
*
*	Solves for the state at height dz above the reference state
*	ref_st in an isentropic one dimensional steady flow.
*
*	The solution is computed by solving the differential
*	equation:
*
*		P_z = rho gz,	P(0) = P_R, rho(0) = rho_r, S = S_r.
*/

LOCAL	void	MPOLY_compute_isentropic_stratified_state(
	Locstate	ans,
	float		dz,	/* distance from reference position */
	float		gz,	/* gravity */
	Locstate	ref_st)
{
	float rho_r, p_r, gam;
	float cr2;
	float csq_ratio;

	Set_params(ans,ref_st);
	gam = Gamma(ref_st);
	cr2 = sound_speed_squared(ref_st);
	csq_ratio = 1.0 + (gam - 1.0)*gz*dz/cr2;
	if (csq_ratio <= 0.0)
	{
	    Dens(ans) = Vacuum_dens(ref_st);
	    Press(ans) = Min_pressure(ref_st);
	}
	else
	{
	    rho_r = Dens(ref_st);
	    p_r = pressure(ref_st);
	    Dens(ans) = rho_r*pow(csq_ratio,1.0/(gam-1.0));
	    Press(ans) = p_r*pow(csq_ratio,gam/(gam-1.0));
	}
	set_type_of_state(ans,TGAS_STATE);
	zero_state_velocity(ans,Params(ref_st)->dim);
}	/*end MPOLY_compute_isentropic_stratified_state*/


/***************End Problem Type Specific Initialization Functions**********/
/***************END INITIALIZATION UTILITY FUNCTIONS************************/


LOCAL	bool	mpoly_fS(float x, float *fans, POINTER prm)
{
	float c, rr, v;
	float ca, cb, dS, alpha, beta, c4, eta, psi;

	ca = ((MPOLY_S_MPRH *) prm)->ca;
	cb = ((MPOLY_S_MPRH *) prm)->cb;
	dS = ((MPOLY_S_MPRH *) prm)->dS;
	alpha = ((MPOLY_S_MPRH *) prm)->alpha;
	beta = ((MPOLY_S_MPRH *) prm)->beta;
	eta = ((MPOLY_S_MPRH *) prm)->eta;
	psi = ((MPOLY_S_MPRH *) prm)->psi;
	c4 = ((MPOLY_S_MPRH *) prm)->c4;

	if (x >= 1.0)
	{
	    rr = (x + c4)/(1.0 + c4*x);
	    c = ca * sqrt(x/rr);
	    v = alpha * (x - 1.0)/sqrt(psi + eta*x);
	}
	else
	{
	    float c2, gam;
	    gam = ((MPOLY_S_MPRH *) prm)->gam;
	    c2 = ((MPOLY_S_MPRH *) prm)->c2;
	    rr = pow(x,1.0/gam);
	    c = ca * pow(x,(gam-1.0)/gam);
	    v = (c - ca)/c2;
	}

	*fans = c - cb + v - (c + cb)*(dS + beta*log(x) - 0.25*log(rr));
	if (is_rotational_symmetry() && ((MPOLY_S_MPRH *) prm)->cyl_coord == YES)
	{
	    float ua, rad, delta;

	    ua = ((MPOLY_S_MPRH *) prm)->ua;
	    rad = ((MPOLY_S_MPRH *) prm)->rad;
	    delta = ((MPOLY_S_MPRH *) prm)->delta;
		
	    *fans = *fans + rad*c*(ua + delta*v);
	}
	return FUNCTION_SUCCEEDED;
}		/*end mpoly_fS*/

LOCAL	float	Gamma(
	Locstate	state)
{
	int	i, nc = Num_gas_components(state);
	float	w;
	float	cv, *gam = gamma(state);
	float	*M = M(state);
	float	mole_dens;

	mole_dens = 0.0;
	cv = 0.0;
	for (i = 0; i < nc; ++i)
	{
	    w = pdens(state)[i]/(M[i]*Dens(state));
	    mole_dens += w;
	    cv += w/(gam[i] - 1.0);
	}
	return 1.0 + mole_dens/cv;
}	/*end Gamma*/

LOCAL   float   MPOLY_MGamma(
        Locstate        state)
{
        int     i, nc = Num_gas_components(state);
        float   w;
        float   cv, *gam = gamma(state);
        float   *M = M(state);
        float   mole_dens;

        mole_dens = 0.0;
        cv = 0.0;
        for (i = 0; i < nc; ++i)
        {
            w = pdens(state)[i]/(M[i]*Dens(state));
            mole_dens += w;
            cv += w/(gam[i] - 1.0);
        }
        return 1.0 + mole_dens/cv;
}       /*end MPOLY_gamma*/

LOCAL	float	Mu(
	Locstate	state)
{
	float	gam = Gamma(state);
	return sqrt((gam-1.0)/(gam+1.0));
}	/*end Mu*/

LOCAL	float	Coef1(
	Locstate	state)
{
	float	gam = Gamma(state);
	return 0.5*(gam + 1.0);
}	/*end Coef1*/

LOCAL	float	Coef2(
	Locstate	state)
{
	float	gam = Gamma(state);
	return 0.5*(gam - 1.0);
}	/*end Coef2*/

LOCAL	float	Coef3(
	Locstate	state)
{
	float	gam = Gamma(state);
	return 0.5*(gam - 1.0)/gam;
}	/*end Coef3*/

LOCAL	float	Coef4(
	Locstate	state)
{
	float	gam = Gamma(state);
	return (gam - 1.0)/(gam + 1.0);
}	/*end Coef4*/

LOCAL	float	Coef5(
	Locstate	state)
{
	float	gam = Gamma(state);
	return 2.0/(gam + 1.0);
}	/*end Coef5*/

LOCAL	float	Coef6(
	Locstate	state)
{
	float	gam = Gamma(state);
	return 1.0/(gam - 1.0);
}	/*end Coef6*/

LOCAL	float	Coef7(
	Locstate	state)
{
	float	gam = Gamma(state);
	return gam/(gam - 1.0);
}	/*end Coef6*/

LOCAL	float	Coef8(
	Locstate	state)
{
	float	gam = Gamma(state);
	return 2.0 * gam/(gam + 1.0);
}	/*end Coef8*/

LOCAL	float	 R(
	Locstate	state)
{
	return MPOLY_Eos(state)->R/molecular_weight(state);
}	/*end R*/

LOCAL	float	molecular_weight(
	Locstate	state)
{
	float	n;
	float	*rho;
	float	*M;
	int	i, n_comps = Num_gas_components(state);

	rho = pdens(state);
	M = M(state);
	n_comps = Num_gas_components(state);

	for (n = 0.0, i = 0; i < n_comps; ++i)
		n += rho[i]/M[i];

	return Dens(state)/n;
}	/*end molecular_weight*/

LOCAL	float	*Vec_Gas_Gamma(
	Vec_Gas *vst,
	int offset,
	int vsize)
{
	Gas_param *params = Params(vst->state[offset]);
	float	*rho = vst->rho + offset;
	float	**rho0 = vst->rho0;
	float	*gam = ((MPOLY_EOS *) params->eos)->_gamma;
	float	*M = ((MPOLY_EOS *) params->eos)->_M;
	int	nc = params->n_comps;
	float	w, cv;
	float	mole_dens;
	int	i, k;
	static	float	*gm = NULL;
	static	int	size_gm = 0;

	if (gm == NULL || (vsize > size_gm))
	{
	    if (gm != NULL)
	    	free(gm);
	    size_gm = vsize;
	    uni_array(&gm,size_gm,FLOAT);
	}

	for (k = 0; k < vsize; ++k)
	{
	    mole_dens = 0.0;
	    cv = 0.0;
	    for (i = 0; i < nc; ++i)
	    {
	    	w = rho0[i][k+offset]/(M[i]*rho[k]);
	    	mole_dens += w;
	    	cv += w/(gam[i] - 1.0);
	    }
	    gm[k] = 1.0 + mole_dens/cv;
	}
	return gm;
}	/*end Vec_Gas_Gamma*/

LOCAL   float   MPOLY_specific_enthalpy_species(
        Locstate state,
        int      i)
{
        float   *gam = gamma(state);
        
        return (gam[i]/(gam[i] - 1.0)) * pressure(state) / Dens(state);
}       /* end MPOLY_specific_enthalpy_species */

LOCAL   float MPOLY_dynamic_viscosity(
        Locstate state,
        float    T)
{
        int i,j,n_comps = Num_gas_components(state);
        float Mass_frac[MAX_NCOMPS];

        /* The lennard-Jones size parameter delta_lennard[n] and the
        Lennard-Jones energy paramenter Tem_lennard*/
        float delta_lennard[MAX_NCOMPS];
        float Tem_lennard[MAX_NCOMPS];
        Tem_lennard[0] = 107.400;
        Tem_lennard[1] = 145.000;
        Tem_lennard[2] = 80.000;
        Tem_lennard[3] = 80.000;
        Tem_lennard[4] = 38.000;
        Tem_lennard[5] = 572.400;
        Tem_lennard[6] = 107.400;
        Tem_lennard[7] = 107.400;
        Tem_lennard[8] = 97.530;
        Tem_lennard[9] = 80.000;
        delta_lennard[0] = 3.458;
        delta_lennard[1] = 2.050;
        delta_lennard[2] = 2.750;
        delta_lennard[3] = 2.750;
        delta_lennard[4] = 2.920;
        delta_lennard[5] = 2.605;
        delta_lennard[6] = 3.458;
        delta_lennard[7] = 3.458;
        delta_lennard[8] = 3.621;
        delta_lennard[9] = 2.750;

        /* using flame master's formula to compute the dynamic/shear viscosity(see the techreport of this part). */
        float pvisc[MAX_NCOMPS];
        float Y[MAX_NCOMPS];
        float G[MAX_NCOMPS][MAX_NCOMPS];
        float delta[MAX_NCOMPS];
        float cl,T_x;
        float PM[MAX_NCOMPS][MAX_NCOMPS],Pmu[MAX_NCOMPS][MAX_NCOMPS];
        float mu = 0.0;
        float Plamda[MAX_NCOMPS];
        float C_P[MAX_NCOMPS];
        float *M = ((MPOLY_EOS *)(Params(state)->eos))->_M;

        for(i =0; i<n_comps;i++)
        {
            Mass_frac[i] = pdens(state)[i]/Dens(state);
            T_x = T/Tem_lennard[i];
            cl = Collision_intergral_2(T_x);
            pvisc[i] = 0.001*2.6693*0.00001*sqrt(M[i]*T)/(delta_lennard[i]*delta_lennard[i]*cl);
        }

        for( i = 0; i< n_comps; i++)
        {
            for(j=0;j<n_comps;j++)
            {
                PM[i][j] = M[i]/M[j];
                Pmu[i][j] =  pvisc[i]/pvisc[j];
            }
        }

        for(i=0;i<n_comps;i++)
            for(j=0;j<n_comps;j++)
                G[i][j] = pow(1 + pow(PM[j][i],0.25)*pow(Pmu[i][j],-0.5),2)*pow(1+PM[i][j],-0.5)*pow(2,-1.5);

        for(i=0;i<n_comps;i++)
        {
            delta[i] = 0;
            for(j=0;j<n_comps;j++)
            {
                delta[i] += G[i][j]*PM[i][j]*Mass_frac[j];
            }
        }

        for(i =0; i<n_comps;i++)
            mu += Mass_frac[i]*pvisc[i]/delta[i];

        return mu;
}

LOCAL   float MPOLY_dynamic_thermal_conductivity(
        Locstate state,
        float    T)
{
        int i,j,n_comps = Num_gas_components(state);
        float Mass_frac[MAX_NCOMPS];

        /* The lennard-Jones size parameter delta_lennard[n] and the
        Lennard-Jones energy paramenter Tem_lennard*/
        float delta_lennard[MAX_NCOMPS];
        float Tem_lennard[MAX_NCOMPS];
        Tem_lennard[0] = 107.400;
        Tem_lennard[1] = 145.000;
        Tem_lennard[2] = 80.000;
        Tem_lennard[3] = 80.000;
        Tem_lennard[4] = 38.000;
        Tem_lennard[5] = 572.400;
        Tem_lennard[6] = 107.400;
        Tem_lennard[7] = 107.400;
        Tem_lennard[8] = 97.530;
        Tem_lennard[9] = 80.000;
        delta_lennard[0] = 3.458;
        delta_lennard[1] = 2.050;
        delta_lennard[2] = 2.750;
        delta_lennard[3] = 2.750;
        delta_lennard[4] = 2.920;
        delta_lennard[5] = 2.605;
        delta_lennard[6] = 3.458;
        delta_lennard[7] = 3.458;
        delta_lennard[8] = 3.621;
        delta_lennard[9] = 2.750;

        /* using flame master's formula to compute the dynamic/shear viscosity(see the techreport of this part). */
        float pvisc[MAX_NCOMPS];
        float Y[MAX_NCOMPS];
        float G[MAX_NCOMPS][MAX_NCOMPS];
        float delta[MAX_NCOMPS];
        float cl,T_x;
        float PM[MAX_NCOMPS][MAX_NCOMPS],Pmu[MAX_NCOMPS][MAX_NCOMPS];
        float lamda = 0.0;
        float Plamda[MAX_NCOMPS];
        float C_P[MAX_NCOMPS];
        float *M = ((MPOLY_EOS *)(Params(state)->eos))->_M;

        for(i =0; i<n_comps;i++)
        {
            Mass_frac[i] = pdens(state)[i]/Dens(state);
            T_x = T/Tem_lennard[i];
            cl = Collision_intergral_2(T_x);
            pvisc[i] = 0.001*2.6693*0.00001*sqrt(M[i]*T)/(delta_lennard[i]*delta_lennard[i]*cl);
        }

        for( i = 0; i< n_comps; i++)
        {
            for(j=0;j<n_comps;j++)
            {
                PM[i][j] = M[i]/M[j];
                Pmu[i][j] =  pvisc[i]/pvisc[j];
            }
        }

        for(i=0;i<n_comps;i++)
            for(j=0;j<n_comps;j++)
                G[i][j] = pow(1 + pow(PM[j][i],0.25)*pow(Pmu[i][j],-0.5),2)*pow(1+PM[i][j],-0.5)*pow(2,-1.5);

        for(i=0;i<n_comps;i++)
        {
            delta[i] = 0;
            for(j=0;j<n_comps;j++)
            {
                delta[i] += G[i][j]*PM[i][j]*Mass_frac[j];
            }
        }

        for(i = 0;i<n_comps;i++)
        {
            C_P[i] = C_P_species(state,T,i);
            Plamda[i] = pvisc[i]*(C_P[i] + 5.0*MPOLY_Eos(state)->R/(4.0*M[i]));
            lamda += Mass_frac[i]*Plamda[i]/delta[i];
        }
 
        return lamda;
}

LOCAL   void MPOLY_dynamic_viscosity_thermalconduct(
        Locstate state,
        float    T,
	float 	*MUrans,
	float 	*Krans)
{
	printf("1\n");
        int i,j,n_comps = Num_gas_components(state);
        float Mass_frac[MAX_NCOMPS];

        /* using flame master's formula to compute the dynamic/shear viscosity(see the techreport of this part). */
        float pvisc[MAX_NCOMPS];
        float Y[MAX_NCOMPS];
        float G[MAX_NCOMPS][MAX_NCOMPS];
        float delta[MAX_NCOMPS];
        float cl,T_x;
        float PM[MAX_NCOMPS][MAX_NCOMPS],Pmu[MAX_NCOMPS][MAX_NCOMPS];
	float mu = 0.0;	
        float lamda = 0.0;
        float Plamda[MAX_NCOMPS];
        float C_P[MAX_NCOMPS];
        float *M = ((MPOLY_EOS *)(Params(state)->eos))->_M;

        /* The lennard-Jones size parameter delta_lennard[n] and the
        Lennard-Jones energy paramenter Tem_lennard*/
        float delta_lennard[MAX_NCOMPS];
        float Tem_lennard[MAX_NCOMPS];
        Tem_lennard[0] = 107.400;
        Tem_lennard[1] = 145.000;
        Tem_lennard[2] = 80.000;
        Tem_lennard[3] = 80.000;
        Tem_lennard[4] = 38.000;
        Tem_lennard[5] = 572.400;
        Tem_lennard[6] = 107.400;
        Tem_lennard[7] = 107.400;
        Tem_lennard[8] = 97.530;
        Tem_lennard[9] = 80.000;
        delta_lennard[0] = 3.458;
        delta_lennard[1] = 2.050;
        delta_lennard[2] = 2.750;
        delta_lennard[3] = 2.750;
        delta_lennard[4] = 2.920;
        delta_lennard[5] = 2.605;
        delta_lennard[6] = 3.458;
        delta_lennard[7] = 3.458;
        delta_lennard[8] = 3.621;
        delta_lennard[9] = 2.750;

	printf("2\n");
        for(i =0; i<n_comps;i++)
        {
            Mass_frac[i] = pdens(state)[i]/Dens(state);
            T_x = T/Tem_lennard[i];
            cl = Collision_intergral_2(T_x);
            pvisc[i] = 0.001*2.6693*0.00001*sqrt(M[i]*T)/(delta_lennard[i]*delta_lennard[i]*cl);
        }
	printf("2.05\n");
        for( i = 0; i< n_comps; i++)
        {
            for(j=0;j<n_comps;j++)
            {
                PM[i][j] = M[i]/M[j];
                Pmu[i][j] =  pvisc[i]/pvisc[j];
            }
        }
	printf("2.1\n");
        for(i=0;i<n_comps;i++)
            for(j=0;j<n_comps;j++)
                G[i][j] = pow(1 + pow(PM[j][i],0.25)*pow(Pmu[i][j],-0.5),2)*pow(1+PM[i][j],-0.5)*pow(2,-1.5);
	printf("2,2\n");
        for(i=0;i<n_comps;i++)
        {
            delta[i] = 0;
            for(j=0;j<n_comps;j++)
            {
                delta[i] += G[i][j]*PM[i][j]*Mass_frac[j];
            }
        }
	printf("2.3\n");
	for(i =0; i<n_comps;i++)
            mu += Mass_frac[i]*pvisc[i]/delta[i];
	*MUrans = mu;
	printf("3\n");
        for(i = 0;i<n_comps;i++)
        {
            C_P[i] = C_P_species(state,T,i);
            Plamda[i] = pvisc[i]*(C_P[i] + 5.0*MPOLY_Eos(state)->R/(4.0*M[i]));
            lamda += Mass_frac[i]*Plamda[i]/delta[i];
        }	
	*Krans = lamda;
	printf("4\n");
}


LOCAL float Collision_intergral_2(float T)
{
  float a = log(T);
  return exp(0.45667 + a*(-0.53955 + a*(0.18265 + a*(-0.03629 + a*(0.00241)))));
}

LOCAL float Collision_intergral_1(float T)
{
  float a = log(T);
  return exp(0.347 + a*(-0.444 + a*(0.093 + a*(-0.010))));
}

LOCAL   float C_P_species(
        Locstate state,
        float    T,
        int i)
{
       float   *gam = gamma(state);

        return R(state)*gam[i]/(gam[i] - 1.0);
}
#endif /* defined(MULTI_COMPONENT) */
