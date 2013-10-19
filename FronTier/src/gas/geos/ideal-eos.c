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
*				ideal-eos.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	File for the implementation of an ideal gas equation of state model
*
*/

#define	DEBUG_STRING	"ideal"
#include <geos/geosdecs.h>
#include <geos/ideal.h>
#define T_crit 1000 //K

/* Generic EOS prototypes */

	/* PRIMARY THERMODYNAMIC FUNCTIONS */
LOCAL	float	IDEAL_internal_energy(Locstate);
LOCAL	float	IDEAL_pressure(Locstate);
LOCAL   float   IDEAL_density(Locstate);
LOCAL	float	IDEAL_sound_speed_squared(Locstate);
LOCAL	float	IDEAL_acoustic_impedance_squared(Locstate);
LOCAL	float	IDEAL_specific_internal_energy(Locstate);

	/* SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS */
LOCAL	float	IDEAL_specific_enthalpy(Locstate);
LOCAL	float	IDEAL_temperature(Locstate);
LOCAL	float	IDEAL_entropy(Locstate);
LOCAL	float	IDEAL_adiabatic_gamma(Locstate);
LOCAL	float	IDEAL_gruneisen_gamma(Locstate);
LOCAL	float	IDEAL_fundamental_derivative(Locstate);
LOCAL	float	IDEAL_C_V(Locstate);
LOCAL	float	IDEAL_C_P(Locstate);
LOCAL	float	IDEAL_K_S(Locstate);
LOCAL	float	IDEAL_K_T(Locstate);
LOCAL   float   IDEAL_specific_enthalpy_species(Locstate,int);
LOCAL   float   IDEAL_dynamic_viscosity(Locstate,float);
LOCAL   float   IDEAL_dynamic_thermal_conductivity(Locstate,float);
LOCAL   void    IDEAL_dynamic_viscosity_thermalconduct(Locstate,float,float*,float*);
LOCAL   float   Collision_intergral_2(float);
LOCAL   float   Collision_intergral_1(float);
LOCAL   float   C_P_species(Locstate,float,int);

	/* OTHER SUPPORTING FUNCTIONS */
LOCAL	float Gamma(Locstate);
LOCAL 	float low_C_P(Locstate);
LOCAL	float high_C_P(Locstate);
LOCAL	float C_P_in_range(Locstate, float);
LOCAL double IDEAL_EOS_newton_step(Locstate, float*, float*, double, double);
LOCAL double IDEAL_EOS_newt_rootFind(Locstate, float*, float*, double, double, double);
LOCAL double IDEAL_EOS_bis_rootFind(Locstate, float*, float*, double, double, double, double);
	/* OTHER SUPPORTING FUNCTIONS */

	/* MATERIAL PROPERTY FUNCTIONS */
LOCAL	float	IDEAL_bulk_viscosity(Locstate);
LOCAL	float	IDEAL_shear_viscosity(Locstate);
LOCAL 	float   IDEAL_heat_coeff(Locstate);

	/* VECTORIZED THERMODYNAMIC FUNCTIONS */
LOCAL	void    IDEAL_single_eos_load_pressure_and_sound_speed2(Vec_Gas*,
								  int,int);
LOCAL	void	IDEAL_single_eos_load_pressure_and_gammas(Vec_Gas*,int,int);
LOCAL	void	IDEAL_single_eos_load_pressure(Vec_Gas*,int,int);
LOCAL	void	IDEAL_single_eos_load_sound_speed2(Vec_Gas*,int,int);

	/* RIEMANN SOLUTIONS UTILITY FUNCTIONS */
	/* Purely Thermodynamic Hugoniot Functions */
LOCAL	float	IDEAL_dens_Hugoniot(float,Locstate);
LOCAL	void	IDEAL_state_w_pr_on_Hugoniot(Locstate,float,Locstate,int);
LOCAL	bool	IDEAL_state_w_mf_sqr_on_Hugoniot(Locstate,float,Locstate,int);

	/* Velocity Related Hugoniot Functions */
LOCAL	float	IDEAL_pr_normal_vel_wave_curve(float,Locstate);

	/* Purely Thermodynamic Adiabatic Wave Curve Functions */
LOCAL	float	IDEAL_dens_rarefaction(float,Locstate);
LOCAL	float	IDEAL_pressure_rarefaction(float,Locstate);
LOCAL	void	IDEAL_state_on_adiabat_with_pr(Locstate,float,Locstate,int);
LOCAL	void	IDEAL_state_on_adiabat_with_dens(Locstate,float,Locstate,int);

	/* General Wave Curve Functions */
LOCAL	float	IDEAL_mass_flux(float,Locstate);
LOCAL	float	IDEAL_mass_flux_squared(float,Locstate);

	/* Functions for the Evaluation of Riemann Solutions */
LOCAL	float	IDEAL_oned_fan_state(float,Locstate,Locstate,Locstate,
				       int,bool*);

	/* Functions to Compute Riemann Solutions */
LOCAL	float	IDEAL_riemann_wave_curve(Locstate,float);
LOCAL	void	IDEAL_set_state_for_find_mid_state(Locstate,Locstate);
LOCAL	float	IDEAL_eps_for_Godunov(Locstate,float,float);
LOCAL	void	IDEAL_initialize_riemann_solver(Locstate,Locstate,float*,
						  float*,float,float*,float*,
						  bool(*)(Locstate,Locstate,
							  float,float*,float*,
							  float*,float*,float*,
							  float*,
							  RIEMANN_SOLVER_WAVE_TYPE*,
							  RIEMANN_SOLVER_WAVE_TYPE*));


	/* TWO DIMENSIONAL RIEMANN SOLUTION UTILTITY FUNCTIONS */
LOCAL	bool	IDEAL_steady_state_wave_curve(float,float,float*,Locstate);
LOCAL	float	IDEAL_pressure_at_sonic_point(float,Locstate);
LOCAL	bool	IDEAL_pr_at_max_turn_angle(float*,float,Locstate);
LOCAL	float	IDEAL_state_in_prandtl_meyer_wave(float,float,Locstate,float,
						    Locstate,Locstate,int);

#if defined(COMBUSTION_CODE)
	/* DETONATION SPECIFIC UTILITY FUNCTIONS */
LOCAL	float	IDEAL_CJ_state(Locstate,int,Locstate,int,int);
LOCAL	void	IDEAL_progress_state(float,Locstate,Locstate,float);
LOCAL	void	IDEAL_fprint_combustion_params(FILE*,Gas_param*);
#endif /* defined(COMBUSTION_CODE) */

	/* METHOD OF CHARACTERISTIC FUNCTIONS FOR W_SPEED */
LOCAL	void	IDEAL_neumann_riem_inv_moc(float*,Locstate,float,float,
	                                     Locstate,SIDE,Locstate,float,
	                                     float*,Front*);
LOCAL	void	IDEAL_shock_ahead_state_riem_inv_moc(float*,Locstate,
						       Locstate,Locstate,
						       Locstate,Locstate,
						       float,float,float,float,
						       float*,float*,int,
						       float,Front*);
LOCAL	bool	IDEAL_shock_moc_plus_rh(float*,Locstate,Locstate,Locstate,
					  Locstate,float,float*,float*,int,
					  Front*);

	/* INITIALIZATION UTILITY FUNCTIONS */
LOCAL	void	IDEAL_prompt_for_state(Locstate,int,Gas_param*,const char*);
LOCAL	void	IDEAL_prompt_for_thermodynamics(Locstate,Gas_param*,
						  const char*);
LOCAL	void	IDEAL_fprint_EOS_params(FILE*,Gas_param*);
LOCAL	void	IDEAL_read_print_EOS_params(INIT_DATA*,const IO_TYPE*,
                                              Gas_param*);
LOCAL	EOS*	IDEAL_free_EOS_params(EOS*);
LOCAL	void	IDEAL_prompt_for_EOS_params(INIT_DATA*,Gas_param*,
					      const char*,const char*);

	/* Problem Type Specific Initialization Functions */
LOCAL	float	IDEAL_RT_RS_f(float,Locstate,float,float,float);
LOCAL	void	IDEAL_RT_single_mode_perturbation_state(Locstate,
							  float*,float,
							  Locstate,float,
							  float,MODE*,float);
LOCAL	void	IDEAL_KH_single_mode_state(Locstate,float*,float,Locstate,
					     float,float,float,MODE*);
LOCAL	void	IDEAL_compute_isothermal_stratified_state(Locstate,float,
							    float,Locstate);
LOCAL	void	IDEAL_compute_isentropic_stratified_state(Locstate,float,
							    float,Locstate);
LOCAL	void	IDEAL_compute_constant_density_stratified_state(Locstate,
							          float,float,
								  Locstate);
LOCAL 	float 	IDEAL_MGamma(Locstate);

	/* Equation of state domain functions */
LOCAL	float	IDEAL_Min_energy(Locstate);
LOCAL	float	IDEAL_Min_pressure(Locstate);
LOCAL	float	IDEAL_Vacuum_dens(Locstate);
LOCAL	float	IDEAL_Raref_press(Locstate);
#if defined(COMBUSTION_CODE)
LOCAL	float	IDEAL_Tol_alpha(Locstate);
LOCAL	float	IDEAL_Tol_pressure(Locstate);
#endif /* defined(COMBUSTION_CODE) */

	/* LOCAL function prototypes */
LOCAL	bool	Mach_number_squared_behind_oblique_shock(float,float*,Locstate);
LOCAL	bool	mta_aux(float,float*,Locstate);
LOCAL	bool	oned_fan_aux(float,float*,POINTER);
LOCAL	bool	pmsr(double,double*,double*,int,Locstate);
LOCAL	bool	pr_normal_vel_wave_curve_aux(float,float*,Locstate);
LOCAL	bool	prandtl_meyer_speed_rinv(float*,float,Locstate,float*,float,
					 Locstate,float*);
LOCAL	bool	pmw_aux(float,float*,Locstate);
LOCAL	bool	rndr(double,double*,double*,int,Locstate);
LOCAL	bool	rnpr(double,double*,double*,int,Locstate);
LOCAL	float	int_dp_over_rho_c(float,Locstate,Locstate);
LOCAL	void	set_eos_function_hooks(EOS*);
LOCAL   float   gaussian_int(int,float,float,Locstate);
LOCAL   void    legendre_init(int,float*,float*);
LOCAL   void    linear_transform(int,float*,float*,float,float);
LOCAL   void    exp_transform(int,float*,float*);
LOCAL   float   one_over_rho_c(float,Locstate);
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
float	molecular_weight(Locstate);
LOCAL	float	R(Locstate);

LOCAL void computeCoeff(Locstate, float*);
LOCAL void computeCoeffs(Locstate, float*, float*);

LOCAL	bool	ideal_Teq(float x, float *fans, POINTER prm);

EXPORT	EOS	*set_IDEAL_eos(
	EOS	*eos)
{
	if (eos == NULL)
	    scalar(&eos,sizeof(IDEAL_EOS));
	(void) set_GENERIC_eos(eos);
	set_eos_function_hooks(eos);
	return eos;
}	/*end set_IDEAL_eos*/

LOCAL	void	set_eos_function_hooks(
	EOS *eos)
{
	/* PRIMARY THERMODYNAMIC FUNCTIONS */
	eos->_internal_energy = IDEAL_internal_energy;
	eos->_pressure = IDEAL_pressure;
	eos->_density = IDEAL_density;
	eos->_sound_speed_squared = IDEAL_sound_speed_squared;
	eos->_acoustic_impedance_squared = IDEAL_acoustic_impedance_squared;
	eos->_specific_internal_energy = IDEAL_specific_internal_energy;
	
	/* MATERIAL PROPERTY FUNCTIONS */
	eos->_bulk_viscosity = IDEAL_bulk_viscosity;
	eos->_shear_viscosity = IDEAL_shear_viscosity;
	eos->_heat_coeff = IDEAL_heat_coeff;

	/* SECONDARY AND SUPPORTING THERMODYNAMIC FUNCTIONS */
	eos->_specific_enthalpy = IDEAL_specific_enthalpy;
	eos->_temperature = IDEAL_temperature;
	eos->_entropy = IDEAL_entropy;
	eos->_adiabatic_gamma = IDEAL_adiabatic_gamma;
	eos->_gruneisen_gamma = IDEAL_gruneisen_gamma;
	eos->_fundamental_derivative = IDEAL_fundamental_derivative;
	eos->_C_V = IDEAL_C_V;
	eos->_C_P = IDEAL_C_P;
	eos->_K_S = IDEAL_K_S;
	eos->_K_T = IDEAL_K_T;
        eos->_specific_enthalpy_species = IDEAL_specific_enthalpy_species;
        eos->_dynamic_viscosity = IDEAL_dynamic_viscosity;
        eos->_dynamic_thermal_conductivity = IDEAL_dynamic_thermal_conductivity;
	eos->_dynamic_viscosity_thermalconduct = IDEAL_dynamic_viscosity_thermalconduct;

	/* VECTORIZED THERMODYNAMIC FUNCTIONS */
	eos->_single_eos_load_pressure_and_sound_speed2 =
	    IDEAL_single_eos_load_pressure_and_sound_speed2;
	eos->_single_eos_load_pressure_and_gammas =
	    IDEAL_single_eos_load_pressure_and_gammas;
	eos->_single_eos_load_pressure = IDEAL_single_eos_load_pressure;
	eos->_single_eos_load_sound_speed2 =
	    IDEAL_single_eos_load_sound_speed2;

	/* RIEMANN SOLUTIONS UTILITY FUNCTIONS */
	/* Purely Thermodynamic Hugoniot Functions */
	eos->_dens_Hugoniot = IDEAL_dens_Hugoniot;
	eos->_state_w_pr_on_Hugoniot = IDEAL_state_w_pr_on_Hugoniot;
	eos->_state_w_mf_sqr_on_Hugoniot = IDEAL_state_w_mf_sqr_on_Hugoniot;

	/* Velocity Related Hugoniot Functions */
	eos->_pr_normal_vel_wave_curve = IDEAL_pr_normal_vel_wave_curve;

	/* Purely Thermodynamic Adiabatic Wave Curve Functions */
	eos->_dens_rarefaction = IDEAL_dens_rarefaction;
	eos->_pressure_rarefaction = IDEAL_pressure_rarefaction;
	eos->_state_on_adiabat_with_pr = IDEAL_state_on_adiabat_with_pr;
	eos->_state_on_adiabat_with_dens = IDEAL_state_on_adiabat_with_dens;

	/* General Wave Curve Functions */
	eos->_mass_flux = IDEAL_mass_flux;
	eos->_mass_flux_squared = IDEAL_mass_flux_squared;

	/* Functions for the Evaluation of Riemann Solutions */
	eos->_oned_fan_state = IDEAL_oned_fan_state;

	/* Functions to Compute Riemann Solutions */
	eos->_riemann_wave_curve = IDEAL_riemann_wave_curve;
	eos->_set_state_for_find_mid_state =
	    IDEAL_set_state_for_find_mid_state;
	eos->_eps_for_Godunov = IDEAL_eps_for_Godunov;
	eos->_initialize_riemann_solver = IDEAL_initialize_riemann_solver;

	/* TWO DIMENSIONAL RIEMANN SOLUTION UTILTITY FUNCTIONS */
	eos->_steady_state_wave_curve = IDEAL_steady_state_wave_curve;
	eos->_pressure_at_sonic_point = IDEAL_pressure_at_sonic_point;
	eos->_pr_at_max_turn_angle = IDEAL_pr_at_max_turn_angle;
	eos->_state_in_prandtl_meyer_wave = IDEAL_state_in_prandtl_meyer_wave;

#if defined(COMBUSTION_CODE)
	/* DETONATION SPECIFIC UTILITY FUNCTIONS */
	eos->_CJ_state = IDEAL_CJ_state;
	eos->_progress_state = IDEAL_progress_state;
	eos->_fprint_combustion_params = IDEAL_fprint_combustion_params;
#endif /* defined(COMBUSTION_CODE) */

	/* METHOD OF CHARACTERISTIC FUNCTIONS FOR W_SPEED */
	eos->_neumann_riem_inv_moc = IDEAL_neumann_riem_inv_moc;
	eos->_shock_ahead_state_riem_inv_moc =
	    IDEAL_shock_ahead_state_riem_inv_moc;
	eos->_shock_moc_plus_rh = IDEAL_shock_moc_plus_rh;

	/* INITIALIZATION UTILITY FUNCTIONS */
	eos->_prompt_for_state = IDEAL_prompt_for_state;
	eos->_prompt_for_thermodynamics = IDEAL_prompt_for_thermodynamics;
	eos->_fprint_EOS_params = IDEAL_fprint_EOS_params;
	eos->_read_print_EOS_params = IDEAL_read_print_EOS_params;
	eos->_free_EOS_params = IDEAL_free_EOS_params;
	eos->_prompt_for_EOS_params = IDEAL_prompt_for_EOS_params;

	/* Problem Type Specific Initialization Functions */
	eos->_RT_RS_f = IDEAL_RT_RS_f;
	eos->_RT_single_mode_perturbation_state = IDEAL_RT_single_mode_perturbation_state;
	eos->_KH_single_mode_state = IDEAL_KH_single_mode_state;
	eos->_compute_isothermal_stratified_state = IDEAL_compute_isothermal_stratified_state;
	eos->_compute_isentropic_stratified_state = IDEAL_compute_isentropic_stratified_state;
	eos->_compute_constant_density_stratified_state = IDEAL_compute_constant_density_stratified_state;
	eos->_multiphase_eos = NO;
	eos->_compute_ns_terms = NO;
	
	eos->_MGamma = IDEAL_MGamma;

	/* Equation of state domain functions */
	eos->_Min_energy = IDEAL_Min_energy;
	eos->_Min_pressure = IDEAL_Min_pressure;
	eos->_Vacuum_dens = IDEAL_Vacuum_dens;
	eos->_Raref_press = IDEAL_Raref_press;
#if defined(COMBUSTION_CODE)
	eos->_Tol_alpha = IDEAL_Tol_alpha;
	eos->_Tol_pressure = IDEAL_Tol_pressure;
#endif /* defined(COMBUSTION_CODE) */
}

float   molecular_weight(
        Locstate        state)
{

        float   n;
        float   *rho;
	IDEAL_EOS *eos = (IDEAL_EOS*)(Params(state)->eos);
        int     i, n_comps = Num_gas_components(state);

        rho = pdens(state);
        n_comps = Num_gas_components(state);

        for (n = 0.0, i = 0; i < n_comps; ++i)
                n += rho[i]/(eos->_M[i]);

	//TMP_XY_weno
	if (Dens(state)/n < 0) {
	    printf("In molecular_weight, Dens(state) = %e\t n = %E\n", Dens(state), n);
//	    for (int ii = 0; ii < n_comps; ++ii)
//		printf("rho[%d] = %E\t M[%d] = %E\n", ii, rho[ii], ii, eos->_M[ii]);
	}

        return Dens(state)/n;

}       /*end molecular_weight*/

LOCAL	float	*Vec_Gas_Gamma(
	Vec_Gas *vst,
	int offset,
	int vsize)
{
	Gas_param *params = Params(vst->state[offset]);
	int	nc = params->n_comps;
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
	//TMP_TK
	gm[0] = Gamma(vst->state[offset]);
	for (k = 1; k < vsize; ++k)
	    gm[k] = gm[0];	

	return gm;
}	/*end Vec_Gas_Gamma*/

LOCAL   float   Coef1(
        Locstate        state)
{
        float   gam = Gamma(state);
        return 0.5*(gam + 1.0);
}       /*end Coef1*/

LOCAL   float   Coef2(
        Locstate        state)
{
        float   gam = Gamma(state);
        return 0.5*(gam - 1.0);
}       /*end Coef2*/

LOCAL   float   Coef3(
        Locstate        state)
{
        float   gam = Gamma(state);
        return 0.5*(gam - 1.0)/gam;
}       /*end Coef3*/

LOCAL   float   Coef4(
        Locstate        state)
{
        float   gam = Gamma(state);
        return (gam - 1.0)/(gam + 1.0);
}       /*end Coef4*/

LOCAL   float   Coef5(
        Locstate        state)
{
        float   gam = Gamma(state);
        return 2.0/(gam + 1.0);
}       /*end Coef5*/

LOCAL   float   Coef6(
        Locstate        state)
{
        float   gam = Gamma(state);
        return 1.0/(gam - 1.0);
}       /*end Coef6*/

LOCAL   float   Coef7(
        Locstate        state)
{
        float   gam = Gamma(state);
        return gam/(gam - 1.0);
}       /*end Coef6*/

LOCAL   float   Coef8(
        Locstate        state)
{
        float   gam = Gamma(state);
        return 2.0 * gam/(gam + 1.0);
}       /*end Coef8*/

float    R(
        Locstate        state)
{
        //return IDEAL_Eos(state)->R/molecular_weight(state);
	IDEAL_EOS *eos = (IDEAL_EOS*)(Params(state)->eos);
	float R = IDEAL_Eos(state)->R;
	float m = molecular_weight(state);
	return R/m;
}       /*end R*/


/*
*			IDEAL_prompt_for_EOS_params():
*
*	Prompts for equation of state parameters.
*/

/*ARGSUSED*/
LOCAL	void	IDEAL_prompt_for_EOS_params(
	INIT_DATA  *init,
	Gas_param  *params,
	const char *message1,
	const char *message2)
{
	IDEAL_EOS *eos = (IDEAL_EOS*) params->eos;
	float *M = eos->_M;
	char s[120];
	FILE *file;
	int i;
	/* Read which components to be used in this simulation */
	screen("Enter the number of components: ");
        (void) Gets(s);
        (void) sscanf(s,"%d", &params->n_comps);

	/* Allocate for components in EOS */
	eos->comp = (IDEAL_EOS_Component*)malloc(params->n_comps*sizeof(IDEAL_EOS_Component));

	/* Read the components' names and mol. mass, the name will be searched in NASA file and used to fill Comp structure */
	for(i=0; i<params->n_comps; i++)
	{
		screen("Enter the component %d name (e.g. O, OH, etc.): ", i+1);
		(void) Gets(s);
		(void) sscanf(s,"%s", name(&eos->comp[i]));
		screen("Enter the molecular mass for component %d: ", i+1);
		(void) Gets(s);
		(void) sscanf(s,"%lf", M+i);
	}

	screen("Enter the ideal gas constant (in proper units): ");
        (void) Gets(s);
        (void) sscanf(s,"%lf", &eos->R);

	screen("Enter the critical temperature for coeff_small and coeff_large: ");
        (void) Gets(s);
        (void) sscanf(s,"%lf", &eos->T);

	screen("Enter the file name for NASA Thermodynamics Data: ");
	(void) Gets(s);
	/*if((s[0]=='\0') || ((file = fopen(s,"r"))==NULL))
	{
		screen("ERROR in NASA Thermodynamic Data, "
		       "can't open %s\n", s);
		clean_up(ERROR);
	}
	fclose(file);*/
	//else
	//{
			//O2
			eos->comp[0]._coeff_large[0] = 3.28253784E+00;
			eos->comp[0]._coeff_large[1] = 1.48308754E-03;
			eos->comp[0]._coeff_large[2] = -7.57966669E-07;
			eos->comp[0]._coeff_large[3] = 2.09470555E-10;
			eos->comp[0]._coeff_large[4] =-2.16717794E-14;
			eos->comp[0]._coeff_large[5] =-1.08845772E+03;
			eos->comp[0]._coeff_large[6] = 5.45323129E+00;
			eos->comp[0]._coeff_small[0] = 3.78245636E+00;
			eos->comp[0]._coeff_small[1] = -2.99673416E-03;
			eos->comp[0]._coeff_small[2] =  9.84730201E-06;
			eos->comp[0]._coeff_small[3] =-9.68129509E-09;
			eos->comp[0]._coeff_small[4] = 3.24372837E-12;
			eos->comp[0]._coeff_small[5] =-1.06394356E+03;
			eos->comp[0]._coeff_small[6] = 3.65767573E+00; 
			
			//H
			eos->comp[1]._coeff_large[0] = 2.50000001E+00;
			eos->comp[1]._coeff_large[1] =-2.30842973E-11;
			eos->comp[1]._coeff_large[2] =1.61561948E-14;
			eos->comp[1]._coeff_large[3] =-4.73515235E-18;
			eos->comp[1]._coeff_large[4] =4.98197357E-22;
			eos->comp[1]._coeff_large[5] =2.54736599E+04;
			eos->comp[1]._coeff_large[6] =-4.46682914E-01;
			eos->comp[1]._coeff_small[0] = 2.50000000E+00;
			eos->comp[1]._coeff_small[1] =7.05332819E-13;
			eos->comp[1]._coeff_small[2] =-1.99591964E-15;
			eos->comp[1]._coeff_small[3] =2.30081632E-18;
			eos->comp[1]._coeff_small[4] =-9.27732332E-22;
			eos->comp[1]._coeff_small[5] =2.54736599E+04;
			eos->comp[1]._coeff_small[6] =-4.46682853E-01;
			
			//OH
			eos->comp[2]._coeff_large[0] =3.09288767E+00;
			eos->comp[2]._coeff_large[1] =5.48429716E-04;
			eos->comp[2]._coeff_large[2] =1.26505228E-07;
			eos->comp[2]._coeff_large[3] =-8.79461556E-11;
			eos->comp[2]._coeff_large[4] =1.17412376E-14;
			eos->comp[2]._coeff_large[5] =3.85865700E+03;
			eos->comp[2]._coeff_large[6] =4.47669610E+00;
			eos->comp[2]._coeff_small[0] =3.99201543E+00;
			eos->comp[2]._coeff_small[1] =-2.40131752E-03;
			eos->comp[2]._coeff_small[2] =4.61793841E-06;
			eos->comp[2]._coeff_small[3] =-3.88113333E-09;
			eos->comp[2]._coeff_small[4] =1.36411470E-12;
			eos->comp[2]._coeff_small[5] =3.61508056E+03;
			eos->comp[2]._coeff_small[6] =-1.03925458E-01;  
			
			//O
			eos->comp[3]._coeff_large[0] =2.56942078E+00;
			eos->comp[3]._coeff_large[1] =-8.59741137E-05;
			eos->comp[3]._coeff_large[2] =4.19484589E-08;
			eos->comp[3]._coeff_large[3] =-1.00177799E-11;
			eos->comp[3]._coeff_large[4] =1.22833691E-15;
			eos->comp[3]._coeff_large[5] =2.92175791E+04;
			eos->comp[3]._coeff_large[6] =4.78433864E+00;
			eos->comp[3]._coeff_small[0] =3.16826710E+00;
			eos->comp[3]._coeff_small[1] =-3.27931884E-03;
			eos->comp[3]._coeff_small[2] =6.64306396E-06;
			eos->comp[3]._coeff_small[3] =-6.12806624E-09;
			eos->comp[3]._coeff_small[4] =2.11265971E-12;
			eos->comp[3]._coeff_small[5] =2.91222592E+04;
			eos->comp[3]._coeff_small[6] =2.05193346E+00;
			
			//H2
			eos->comp[4]._coeff_large[0] =3.33727920E+00;
			eos->comp[4]._coeff_large[1] =-4.94024731E-05;
			eos->comp[4]._coeff_large[2] =4.99456778E-07;
			eos->comp[4]._coeff_large[3] =-1.79566394E-10;
			eos->comp[4]._coeff_large[4] =2.00255376E-14;
			eos->comp[4]._coeff_large[5] =-9.50158922E+02;
			eos->comp[4]._coeff_large[6] =-3.20502331E+00;
			eos->comp[4]._coeff_small[0] =2.34433112E+00;
			eos->comp[4]._coeff_small[1] =7.98052075E-03;
			eos->comp[4]._coeff_small[2] =-1.94781510E-05;
			eos->comp[4]._coeff_small[3] =2.01572094E-08;
			eos->comp[4]._coeff_small[4] =-7.37611761E-12;
			eos->comp[4]._coeff_small[5] =-9.17935173E+02;
			eos->comp[4]._coeff_small[6] =6.83010238E-01; 
			
			//H2O
			eos->comp[5]._coeff_large[0]=3.03399249E+00;
			eos->comp[5]._coeff_large[1]=2.17691804E-03;
			eos->comp[5]._coeff_large[2]=-1.64072518E-07;
			eos->comp[5]._coeff_large[3]=-9.70419870E-11;
			eos->comp[5]._coeff_large[4]=1.68200992E-14;
			eos->comp[5]._coeff_large[5]=-3.00042971E+04;
			eos->comp[5]._coeff_large[6]=4.96677010E+00;
			eos->comp[5]._coeff_small[0]=4.19864056E+00;
			eos->comp[5]._coeff_small[1]=-2.03643410E-03;
			eos->comp[5]._coeff_small[2]=6.52040211E-06;
			eos->comp[5]._coeff_small[3]=-5.48797062E-09;
			eos->comp[5]._coeff_small[4]=1.77197817E-12;
			eos->comp[5]._coeff_small[5]=-3.02937267E+04;
			eos->comp[5]._coeff_small[6]=-8.49032208E-01;
			
			//HO2
			eos->comp[6]._coeff_large[0]=4.01721090E+00;
			eos->comp[6]._coeff_large[1]=2.23982013E-03;
			eos->comp[6]._coeff_large[2]=-6.33658150E-07;
			eos->comp[6]._coeff_large[3]=1.14246370E-10;
			eos->comp[6]._coeff_large[4]=-1.07908535E-14;
			eos->comp[6]._coeff_large[5]=1.11856713E+02;
			eos->comp[6]._coeff_large[6]=3.78510215E+00;
			eos->comp[6]._coeff_small[0]=4.30179801E+00;
			eos->comp[6]._coeff_small[1]=-4.74912051E-03;
			eos->comp[6]._coeff_small[2]=2.11582891E-05;
			eos->comp[6]._coeff_small[3]=-2.42763894E-08;
			eos->comp[6]._coeff_small[4]=9.29225124E-12;
			eos->comp[6]._coeff_small[5]=2.94808040E+02;
			eos->comp[6]._coeff_small[6]=3.71666245E+00; 
			
			//H2O2
			eos->comp[7]._coeff_large[0]=4.16500285E+00;
			eos->comp[7]._coeff_large[1]=4.90831694E-03;
			eos->comp[7]._coeff_large[2]=-1.90139225E-06;
			eos->comp[7]._coeff_large[3]=3.71185986E-10;
			eos->comp[7]._coeff_large[4]=-2.87908305E-14;
			eos->comp[7]._coeff_large[5]=-1.78617877E+04;
			eos->comp[7]._coeff_large[6]=2.91615662E+00;
			eos->comp[7]._coeff_small[0]=4.27611269E+00;
			eos->comp[7]._coeff_small[1]=-5.42822417E-04;
			eos->comp[7]._coeff_small[2]=1.67335701E-05;
			eos->comp[7]._coeff_small[3]=-2.15770813E-08;
			eos->comp[7]._coeff_small[4]=8.62454363E-12;
			eos->comp[7]._coeff_small[5]=-1.77025821E+04;
			eos->comp[7]._coeff_small[6]=3.43505074E+00;

			//N2
			eos->comp[8]._coeff_large[0]=0.0;
			eos->comp[8]._coeff_large[1]=0.0;
			eos->comp[8]._coeff_large[2]=0.0;
			eos->comp[8]._coeff_large[3]=0.0;
			eos->comp[8]._coeff_large[4]=0.0;
			eos->comp[8]._coeff_large[5]=0.0;
			eos->comp[8]._coeff_large[6]=0.0;
			eos->comp[8]._coeff_small[0]=0.0;
			eos->comp[8]._coeff_small[1]=0.0;
			eos->comp[8]._coeff_small[2]=0.0;
			eos->comp[8]._coeff_small[3]=0.0;
			eos->comp[8]._coeff_small[4]=0.0;
			eos->comp[8]._coeff_small[5]=0.0;
			eos->comp[8]._coeff_small[6]=0.0; 
			
			//OHS
			eos->comp[9]._coeff_large[0]=2.75582920E+00;
			eos->comp[9]._coeff_large[1]=1.39848756E-03;
			eos->comp[9]._coeff_large[2]=-4.19428493E-07;
			eos->comp[9]._coeff_large[3]=6.33453282E-11;
			eos->comp[9]._coeff_large[4]=-3.56042218E-15;
			eos->comp[9]._coeff_large[5]=5.09751756E+04;
			eos->comp[9]._coeff_large[6]=5.62581429E+00;
			eos->comp[9]._coeff_small[0]=3.46084428E+00;
			eos->comp[9]._coeff_small[1]=5.01872172E-04;
			eos->comp[9]._coeff_small[2]=-2.00254474E-06;
			eos->comp[9]._coeff_small[3]=3.18901984E-09;
			eos->comp[9]._coeff_small[4]=-1.35451838E-12;
			eos->comp[9]._coeff_small[5]=5.07349466E+04;
			eos->comp[9]._coeff_small[6]=1.73976415E+00; 

	//}
}		/*end IDEAL_prompt_for_EOS_params*/

/*
*			IDEAL_read_print_EOS_params():
*
*	Reads the equation of state data from a NASA database.
*/

/*ARGSUSED*/
LOCAL	void	IDEAL_read_print_EOS_params(
	INIT_DATA     *init,
	const IO_TYPE *io_type,
	Gas_param     *params)
{
	FILE      *file = io_type->file;
	IDEAL_EOS *eos = (IDEAL_EOS *)params->eos;
	float	  *M = eos->_M;
	int	  c, i, nc;

	(void) fgetstring(file,"number of components = ");
	(void) fscanf(file,"%d",&nc);
	params->n_comps = nc;

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

	eos->R = fread_float("R = ",io_type);
	eos->T = fread_float("Tcrit = ",io_type);

	eos->comp = (IDEAL_EOS_Component*)malloc(params->n_comps*sizeof(IDEAL_EOS_Component));
	//O2
	eos->comp[0]._coeff_large[0] = 3.28253784E+00;
	eos->comp[0]._coeff_large[1] = 1.48308754E-03;
	eos->comp[0]._coeff_large[2] = -7.57966669E-07;
	eos->comp[0]._coeff_large[3] = 2.09470555E-10;
	eos->comp[0]._coeff_large[4] =-2.16717794E-14;
	eos->comp[0]._coeff_large[5] =-1.08845772E+03;
	eos->comp[0]._coeff_large[6] = 5.45323129E+00;
	eos->comp[0]._coeff_small[0] = 3.78245636E+00;
	eos->comp[0]._coeff_small[1] = -2.99673416E-03;
	eos->comp[0]._coeff_small[2] =  9.84730201E-06;
	eos->comp[0]._coeff_small[3] =-9.68129509E-09;
	eos->comp[0]._coeff_small[4] = 3.24372837E-12;
	eos->comp[0]._coeff_small[5] =-1.06394356E+03;
	eos->comp[0]._coeff_small[6] = 3.65767573E+00; 
			
	//H
	eos->comp[1]._coeff_large[0] = 2.50000001E+00;
	eos->comp[1]._coeff_large[1] =-2.30842973E-11;
	eos->comp[1]._coeff_large[2] =1.61561948E-14;
	eos->comp[1]._coeff_large[3] =-4.73515235E-18;
	eos->comp[1]._coeff_large[4] =4.98197357E-22;
	eos->comp[1]._coeff_large[5] =2.54736599E+04;
	eos->comp[1]._coeff_large[6] =-4.46682914E-01;
	eos->comp[1]._coeff_small[0] = 2.50000000E+00;
	eos->comp[1]._coeff_small[1] =7.05332819E-13;
	eos->comp[1]._coeff_small[2] =-1.99591964E-15;
	eos->comp[1]._coeff_small[3] =2.30081632E-18;
	eos->comp[1]._coeff_small[4] =-9.27732332E-22;
	eos->comp[1]._coeff_small[5] =2.54736599E+04;
	eos->comp[1]._coeff_small[6] =-4.46682853E-01;
	
	//OH
	eos->comp[2]._coeff_large[0] =3.09288767E+00;
	eos->comp[2]._coeff_large[1] =5.48429716E-04;
	eos->comp[2]._coeff_large[2] =1.26505228E-07;
	eos->comp[2]._coeff_large[3] =-8.79461556E-11;
	eos->comp[2]._coeff_large[4] =1.17412376E-14;
	eos->comp[2]._coeff_large[5] =3.85865700E+03;
	eos->comp[2]._coeff_large[6] =4.47669610E+00;
	eos->comp[2]._coeff_small[0] =3.99201543E+00;
	eos->comp[2]._coeff_small[1] =-2.40131752E-03;
	eos->comp[2]._coeff_small[2] =4.61793841E-06;
	eos->comp[2]._coeff_small[3] =-3.88113333E-09;
	eos->comp[2]._coeff_small[4] =1.36411470E-12;
	eos->comp[2]._coeff_small[5] =3.61508056E+03;
	eos->comp[2]._coeff_small[6] =-1.03925458E-01;  
	
	//O
	eos->comp[3]._coeff_large[0] =2.56942078E+00;
	eos->comp[3]._coeff_large[1] =-8.59741137E-05;
	eos->comp[3]._coeff_large[2] =4.19484589E-08;
	eos->comp[3]._coeff_large[3] =-1.00177799E-11;
	eos->comp[3]._coeff_large[4] =1.22833691E-15;
	eos->comp[3]._coeff_large[5] =2.92175791E+04;
	eos->comp[3]._coeff_large[6] =4.78433864E+00;
	eos->comp[3]._coeff_small[0] =3.16826710E+00;
	eos->comp[3]._coeff_small[1] =-3.27931884E-03;
	eos->comp[3]._coeff_small[2] =6.64306396E-06;
	eos->comp[3]._coeff_small[3] =-6.12806624E-09;
	eos->comp[3]._coeff_small[4] =2.11265971E-12;
	eos->comp[3]._coeff_small[5] =2.91222592E+04;
	eos->comp[3]._coeff_small[6] =2.05193346E+00;
	
	//H2
	eos->comp[4]._coeff_large[0] =3.33727920E+00;
	eos->comp[4]._coeff_large[1] =-4.94024731E-05;
	eos->comp[4]._coeff_large[2] =4.99456778E-07;
	eos->comp[4]._coeff_large[3] =-1.79566394E-10;
	eos->comp[4]._coeff_large[4] =2.00255376E-14;
	eos->comp[4]._coeff_large[5] =-9.50158922E+02;
	eos->comp[4]._coeff_large[6] =-3.20502331E+00;
	eos->comp[4]._coeff_small[0] =2.34433112E+00;
	eos->comp[4]._coeff_small[1] =7.98052075E-03;
	eos->comp[4]._coeff_small[2] =-1.94781510E-05;
	eos->comp[4]._coeff_small[3] =2.01572094E-08;
	eos->comp[4]._coeff_small[4] =-7.37611761E-12;
	eos->comp[4]._coeff_small[5] =-9.17935173E+02;
	eos->comp[4]._coeff_small[6] =6.83010238E-01; 
	
	//H2O
	eos->comp[5]._coeff_large[0]=3.03399249E+00;
	eos->comp[5]._coeff_large[1]=2.17691804E-03;
	eos->comp[5]._coeff_large[2]=-1.64072518E-07;
	eos->comp[5]._coeff_large[3]=-9.70419870E-11;
	eos->comp[5]._coeff_large[4]=1.68200992E-14;
	eos->comp[5]._coeff_large[5]=-3.00042971E+04;
	eos->comp[5]._coeff_large[6]=4.96677010E+00;
	eos->comp[5]._coeff_small[0]=4.19864056E+00;
	eos->comp[5]._coeff_small[1]=-2.03643410E-03;
	eos->comp[5]._coeff_small[2]=6.52040211E-06;
	eos->comp[5]._coeff_small[3]=-5.48797062E-09;
	eos->comp[5]._coeff_small[4]=1.77197817E-12;
	eos->comp[5]._coeff_small[5]=-3.02937267E+04;
	eos->comp[5]._coeff_small[6]=-8.49032208E-01;
	
	//HO2
	eos->comp[6]._coeff_large[0]=4.01721090E+00;
	eos->comp[6]._coeff_large[1]=2.23982013E-03;
	eos->comp[6]._coeff_large[2]=-6.33658150E-07;
	eos->comp[6]._coeff_large[3]=1.14246370E-10;
	eos->comp[6]._coeff_large[4]=-1.07908535E-14;
	eos->comp[6]._coeff_large[5]=1.11856713E+02;
	eos->comp[6]._coeff_large[6]=3.78510215E+00;
	eos->comp[6]._coeff_small[0]=4.30179801E+00;
	eos->comp[6]._coeff_small[1]=-4.74912051E-03;
	eos->comp[6]._coeff_small[2]=2.11582891E-05;
	eos->comp[6]._coeff_small[3]=-2.42763894E-08;
	eos->comp[6]._coeff_small[4]=9.29225124E-12;
	eos->comp[6]._coeff_small[5]=2.94808040E+02;
	eos->comp[6]._coeff_small[6]=3.71666245E+00; 
	
	//H2O2
	eos->comp[7]._coeff_large[0]=4.16500285E+00;
	eos->comp[7]._coeff_large[1]=4.90831694E-03;
	eos->comp[7]._coeff_large[2]=-1.90139225E-06;
	eos->comp[7]._coeff_large[3]=3.71185986E-10;
	eos->comp[7]._coeff_large[4]=-2.87908305E-14;
	eos->comp[7]._coeff_large[5]=-1.78617877E+04;
	eos->comp[7]._coeff_large[6]=2.91615662E+00;
	eos->comp[7]._coeff_small[0]=4.27611269E+00;
	eos->comp[7]._coeff_small[1]=-5.42822417E-04;
	eos->comp[7]._coeff_small[2]=1.67335701E-05;
	eos->comp[7]._coeff_small[3]=-2.15770813E-08;
	eos->comp[7]._coeff_small[4]=8.62454363E-12;
	eos->comp[7]._coeff_small[5]=-1.77025821E+04;
	eos->comp[7]._coeff_small[6]=3.43505074E+00;

	//N2
	eos->comp[8]._coeff_large[0]=0.0;
	eos->comp[8]._coeff_large[1]=0.0;
	eos->comp[8]._coeff_large[2]=0.0;
	eos->comp[8]._coeff_large[3]=0.0;
	eos->comp[8]._coeff_large[4]=0.0;
	eos->comp[8]._coeff_large[5]=0.0;
	eos->comp[8]._coeff_large[6]=0.0;
	eos->comp[8]._coeff_small[0]=0.0;
	eos->comp[8]._coeff_small[1]=0.0;
	eos->comp[8]._coeff_small[2]=0.0;
	eos->comp[8]._coeff_small[3]=0.0;
	eos->comp[8]._coeff_small[4]=0.0;
	eos->comp[8]._coeff_small[5]=0.0;
	eos->comp[8]._coeff_small[6]=0.0; 
	
	//OHS
	eos->comp[9]._coeff_large[0]=2.75582920E+00;
	eos->comp[9]._coeff_large[1]=1.39848756E-03;
	eos->comp[9]._coeff_large[2]=-4.19428493E-07;
	eos->comp[9]._coeff_large[3]=6.33453282E-11;
	eos->comp[9]._coeff_large[4]=-3.56042218E-15;
	eos->comp[9]._coeff_large[5]=5.09751756E+04;
	eos->comp[9]._coeff_large[6]=5.62581429E+00;
	eos->comp[9]._coeff_small[0]=3.46084428E+00;
	eos->comp[9]._coeff_small[1]=5.01872172E-04;
	eos->comp[9]._coeff_small[2]=-2.00254474E-06;
	eos->comp[9]._coeff_small[3]=3.18901984E-09;
	eos->comp[9]._coeff_small[4]=-1.35451838E-12;
	eos->comp[9]._coeff_small[5]=5.07349466E+04;
	eos->comp[9]._coeff_small[6]=1.73976415E+00; 


}		/*end IDEAL_read_print_EOS_params*/

/* Solve for T: (a+bT+cT^2+dT^3+eT^4-1)(RT/M) - E=0 */
LOCAL double IDEAL_EOS_eqn(Locstate state, float* c_small, float* c_large, double E, double T)
{
        double* coeff;
	IDEAL_EOS *eos = (IDEAL_EOS*)(Params(state)->eos);
	double R_GAS = R(state);
	double T_MID = eos->T;
	float M = molecular_weight(state);

        if(T<=T_MID)
		coeff = c_small;
        else
		coeff = c_large;
	double term1 = (coeff[0]+T*(coeff[1]+T*(coeff[2]+T*(coeff[3]+T*coeff[4])))-1)*R_GAS*T;
        return (coeff[0]+T*(coeff[1]+T*(coeff[2]+T*(coeff[3]+T*coeff[4])))-1)*R_GAS*T-E;
}


LOCAL double IDEAL_EOS_newton_step(Locstate state, float* c_small, float* c_large, double E, double T)
{
	double *c;
	IDEAL_EOS *eos = (IDEAL_EOS*)(Params(state)->eos);
	double T_MID = eos->T;
	
	if (T <= T_MID)
		c = c_small;
	else
		c = c_large;
	
	double f, df;
	f = R(state)*T*( c[0] + (c[1]*T) + (c[2]*T*T) + (c[3]*T*T*T) + (c[4]*T*T*T*T) - 1 ) - E;
	df = R(state)*(c[0] + (2*c[1]*T) + (3*c[2]*T*T) + (4*c[3]*T*T*T) + (5*c[4]*T*T*T*T) - 1);

	return T - (f/df) ;
}

LOCAL double IDEAL_EOS_newt_rootFind(Locstate state, float* c_small, float* c_large, double E, double Tinit, double err)
{
	float Tn, T;
	T = Tinit;
	Tn = IDEAL_EOS_newton_step(state, c_small, c_large, E, T);
	int counter = 0;

	do
	{
		T = Tn;
		Tn = IDEAL_EOS_newton_step(state, c_small, c_large, E, T);
		counter++;
		if( counter > 10000 )
		{
			return IDEAL_EOS_bis_rootFind(state, c_small, c_large, E, 0, 3500, err);
			break;
		}
	}while(fabs(Tn-T) > err);
	
	return .5*(Tn+T);
	
}

LOCAL double IDEAL_EOS_bis_rootFind(Locstate state, float* c_small, float* c_large, double E, double Tmin, double Tmax, double err)
{
        double l = IDEAL_EOS_eqn(state, c_small, c_large, E, Tmin);
        double r = IDEAL_EOS_eqn(state, c_small, c_large, E, Tmax);
	double m;
	
	while( fabs(Tmax-Tmin) > err)
	{
		m = IDEAL_EOS_eqn(state, c_small, c_large, E, (Tmin+Tmax)*.5 );
		if (m < 0)
		{
			Tmin = (Tmin+Tmax)*.5;
			l = m;
		}else{
			Tmax = (Tmin+Tmax)*.5;
			r = m;
		}
	}
	
	return .5*(Tmax+Tmin);
}

/*
*			IDEAL_fprint_EOS_params():
*
*	Prints the parameters that define the given equation of state.
*	NOTE:  This is not really an initialization function,  but it is
*	convenient to locate it next the the corresponding read function.
*/

/*ARGSUSED*/
LOCAL	void	IDEAL_fprint_EOS_params(
	FILE *file,
	Gas_param *params)
{
	float	R = ((IDEAL_EOS *) params->eos)->R;
	float 	Tcrit = ((IDEAL_EOS *) params->eos)->T;
	IDEAL_EOS *eos = (IDEAL_EOS *)params->eos;
	int	i, nc = params->n_comps;
	float	M[nc];

	(void) fprintf(file,"\tEquation of state = %d MULTI_COMP_IDEAL_GAS\n",
		MULTI_COMP_IDEAL);
	(void) fprintf(file,"\tnumber of components = %d\n",params->n_comps);
	(void) fprintf(file,"\tM = ");
	for (i = 0; i < nc; ++i)
		M[i] = eos->_M[i];
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
	(void) fprintf(file," R = ");
	if (is_binary_output() == YES)
	{
	    (void) fprintf(file,"\f%c",1);
	    (void) fwrite((const void *) &R,FLOAT,1,file);
	}
	else
	    (void) fprintf(file,"%"FFMT,R);

	(void) fprintf(file," Tcrit = ");
	if (is_binary_output() == YES)
	{
	     (void) fprintf(file, "\f%c",1);
	     (void) fwrite((const void *) &Tcrit, FLOAT,1,file);
	}
	else
		(void)	fprintf(file, "%"FFMT,Tcrit);


	(void) fprintf(file,"\n");
}		/*end IDEAL_fprint_EOS_params */


/*
*                       IDEAL_free_EOS_params():
*
*       Frees the storage allocated for an equation of state parameter
*       function.
*/

LOCAL   EOS*    IDEAL_free_EOS_params(
        EOS *eos)
{
        free(eos);
        return NULL;
}               /*end IDEAL_free_EOS_params*/


	/*end IDEAL_prompt_for_EOS_params*/

/*
*			IDEAL_fundamental_derivative():
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

LOCAL	float	IDEAL_fundamental_derivative(
	Locstate state)
{
	screen("ERROR in IDEAL_fundamental_derivative(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end IDEAL_fundamental_derivative*/

/*
*			IDEAL_C_P():
*
*	Specific heat at constant pressure.
*
*
*                        dS  |
*		C_P = T ---- |
*                        dT  | P
*		For ideal gas, C_P is calculated using the NASA formula. Refer to
*		Jones and Hawkins, Engineering Thermodynamics. 
*/

/* Computes specific heat at constant pressure for states with temperatures 
   below the lower limit of the NASA tables.  It computes specific heat
   at T = 200. */
LOCAL	float low_C_P(
	Locstate state)
{
	float cp, T;
	float coeff[7];
	
	//lower bound for NASA data.
	T = 200;
	computeCoeff(state, coeff);

	cp = R(state)*(coeff[0]+T*(coeff[1]+T*(coeff[2]+T*(coeff[3]+T*(coeff[4])))));
	
	return cp;
}

LOCAL	float high_C_P(
	Locstate state)
{
    	float cp, T;
	float coeff[7];

	//higher bound for NASA data.
	T = 3500;
	computeCoeff(state, coeff);

	cp = R(state)*(coeff[0]+T*(coeff[1]+T*(coeff[2]+T*(coeff[3]+T*(coeff[4])))));
	return cp;
}

LOCAL	float	C_P_in_range(
	Locstate state,
	float	T)
{
    	float coeff[7];
	computeCoeff(state, coeff);
	float cp = R(state)*(coeff[0]+coeff[1]*T+coeff[2]*T*T+coeff[3]*T*T*T+coeff[4]*T*T*T*T);
	return cp;
}

LOCAL	float	IDEAL_C_P(
	Locstate state)
{
	float coeff[7];
	float T = temperature(state);
	
	computeCoeff(state, coeff);

	float cp = R(state)*(coeff[0]+coeff[1]*T+coeff[2]*T*T+coeff[3]*T*T*T+coeff[4]*T*T*T*T);

	return cp;
}
	/* end IDEAL_C_P */


LOCAL	float C_P_species(
	Locstate state,
        float    T,
	int i)
{
	float coeff[7];
	float r;
        int j;
	
	IDEAL_EOS *eos = (IDEAL_EOS *)(Params(state)->eos);
	r = (IDEAL_Eos(state)->R);
	r = r/ (eos->_M[i]);

        if( T < eos->T)
        {
                for ( j=0; j<7; j++ )
                        coeff[j] = coeff_small(IDEAL_Eos_Component(state,i),j);

        }else{
                for ( j=0; j<7; j++ )
                        coeff[j] = coeff_large(IDEAL_Eos_Component(state,i),j);
        }

	return r*(coeff[0]+coeff[1]*T+coeff[2]*T*T+coeff[3]*T*T*T+coeff[4]*T*T*T*T);
}


/*
*			IDEAL_C_V():
*
*		BASIC TEMPERATURE DEPENDENT FUNCTION
*
*	Specific heat at constant volume.
*
*                        dS  |
*		C_V = T ---- |
*                        dT  | V
*		For ideal gas, C_V is calculated using the formula C_V = C_P  - R. Refer to
*		Jones and Hawkins, Engineering Thermodynamics. 
*/

/*ARGSUSED*/
LOCAL	float	IDEAL_C_V(
	Locstate state)
{
	return IDEAL_C_P(state) - R(state);
}	/* end IDEAL_C_V */

/*
*			IDEAL_K_S():
*
*	Isentropic compressibility.
*
*                        1   dV  |
*		K_S = - --- ---- |
*                        V   dP  | S
*/

LOCAL	float	IDEAL_K_S(
	Locstate state)
{
	screen("ERROR in IDEAL_K_S(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}	/* end IDEAL_K_S */

/*
*			IDEAL_K_T():
*
*	Isothermal compressibility.
*
*                        1   dV  |
*		K_T = - --- ---- |
*                        V   dP  | T
*/

LOCAL	float	IDEAL_K_T(
	Locstate state)
{
	screen("ERROR in IDEAL_K_T(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}	/*end IDEAL_K_T */

LOCAL float Gamma(
      Locstate state)
{
	float cp, cv, T;
	T = temperature(state);
	if (T < 200)
		cp = low_C_P(state);
	else if (T > 3500)
	    	cp = high_C_P(state);
	else	
		cp = C_P_in_range(state, T);

	cv = cp - R(state); 
	
	return cp/cv;
}

LOCAL float IDEAL_MGamma(
      Locstate state)
{
	return Gamma(state);
}

/*
*			IDEAL_pressure():
*
*		BASIC FUNCTION NO DEFAULT IMPLEMENTATION
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

LOCAL	float	IDEAL_pressure(
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
	    screen("ERROR in IDEAL_pressure(), no such state type\n");
	    clean_up(ERROR);
	}
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	pr = max(pr,Min_pressure(state));
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	return pr;
}		/*end IDEAL_pressure*/

LOCAL   float   IDEAL_density(
        Locstate state)
{	
        switch (state_type(state))
	{
	case	GAS_STATE:
	case	EGAS_STATE:
	case	FGAS_STATE:
	case	TGAS_STATE:
	case	VGAS_STATE:
	        return Dens(state);
	default:
	        screen("ERROR in IDEAL_density(), function unimplemented\n");
		clean_up(ERROR);
		break;
	}
	return ERROR_FLOAT;
}               /*end IDEAL_density*/

/*
*			IDEAL_sound_speed_squared():
*
*		BASIC FUNCTION NO DEFAULT IMPLEMENTATION
*
*	Returns the square of the local sound speed of the state.
*
*                        2   dP  |
*			c = ---- |
*                           drho |S
*/

/*ARGSUSED*/
LOCAL	float	IDEAL_sound_speed_squared(
	Locstate state)
{
        if (state_type(state) == VGAS_STATE)
        {
            float c = Sound_speed(state);
            return c*c;
        }
	else
//		return Gamma(state)*R(state)*temperature(state); //TMP_TK
		return Gamma(state)*pressure(state)/Dens(state); //TMP_TK
	
}		/*end IDEAL_sound_speed_squared*/

/*
*			IDEAL_specific_enthalpy():
*
*	This function computes the specific enthalpy of the given state.
*
*			H = E + P*V = E + RT
*
*	E = specific internal energy, P = pressure, V = specific volume,
*	R = Universal Gas Constant, T = Temperature.
*	Refer to Engineering Thermodynamics by Jones and Hawkins.  
*
*/

LOCAL	float	IDEAL_specific_enthalpy(
	Locstate state)
{
	float coeff[7];
	float T = temperature(state);

	computeCoeff(state, coeff);

#if defined(VERBOSE_PLUS_GAS)
	if (state_type(state) == VGAS_STATE)
	    return Enthalpy(state);
#endif /* defined(VERBOSE_PLUS_GAS) */

	return R(state)*T*(coeff[0]+coeff[1]*T/2+coeff[2]*T*T/3+coeff[3]*T*T*T/4+
				 coeff[4]*T*T*T*T/5+coeff[5]/T);
}		/*end IDEAL_specific_enthalpy*/
	

/*
*			IDEAL_specific_internal_energy():
*
*		BASIC FUNCTION NO DEFAULT IMPLEMENTATION
*
*	Returns the specific internal energy = internal energy per unit
*	mass of the state.
*/
LOCAL   float   IDEAL_specific_internal_energy(
        Locstate state)
{
        float rho = Dens(state);
	float tmp;

        switch (state_type(state))
        {
        case GAS_STATE:
            return      (Energy(state) - kinetic_energy(state))/rho;

        case EGAS_STATE:
            return      Energy(state);

        case TGAS_STATE:
            return      Coef6(state)* Press(state) / rho;

        case FGAS_STATE:
            return      R(state)*Temperature(state)*Coef6(state);

        case VGAS_STATE:
            return      Int_en(state);

        default:
            screen("ERROR in IDEAL_specific_internal_energy(), "
                   "no such state type\n");
            clean_up(ERROR);
            break;
        }
        return ERROR_FLOAT;

}               /*end IDEAL_specific_internal_energy*/


/*
*			IDEAL_gruneisen_gamma():
*
*		BASIC FUNCTION NO DEFAULT IMPLEMENTATION
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

/*ARGSUSED*/
LOCAL	float	IDEAL_gruneisen_gamma(
	Locstate state)
{
	return Gamma(state) - 1.0;
}		/*end IDEAL_gruneisen_gamma*/

/*
*		BASIC TEMPERATURE DEPENDENT FUNCTIONS
*		NO DEFAULT IMPLEMENTATION AVAILABLE,
*		BUT MY NOT BE NEEDED FOR ALL APPLICATIONS
*/

LOCAL void computeCoeffs(Locstate state, float* c_small, float* c_large)
{
	IDEAL_EOS *eos = (IDEAL_EOS*)(Params(state)->eos);
        float* rho = pdens(state);
	int n, x;

	for(x=0; x<7; x++)
		c_small[x] = c_large[x] = 0;
	
	for(n=0; n<Params(state)->n_comps; n++)
	{
		float mf = rho[n]/Dens(state);
		if(mf!=0)
		{
			for(x=0; x<7; x++)
			{
				c_small[x] += coeff_small(IDEAL_Eos_Component(state, n), x)*mf;
				c_large[x] += coeff_large(IDEAL_Eos_Component(state, n), x)*mf;
			}
		}
	}
}

LOCAL void computeCoeff(Locstate state, float* coeff)
{
	float c_small[7];
	float c_large[7];
	IDEAL_EOS *eos = (IDEAL_EOS*)(Params(state)->eos);

	computeCoeffs(state, c_small, c_large);
	if(temperature(state)<eos->T)
		memcpy(coeff, c_small, 7*sizeof(float));
	else
		memcpy(coeff, c_large, 7*sizeof(float));
}

/*
*			IDEAL_temperature():
*
*		BASIC TEMPERATURE DEPENDENT FUNCTION
*
*	Returns the thermodynamic temperature of a state.
*
*                            dE |
*			T = --- |
*                            dS |V
*/

LOCAL	float	IDEAL_temperature(
	Locstate state)
{
	float T;
	if (state_type(state) == FGAS_STATE)
	{
		T = Temperature(state);
	}else if (state_type(state) == TGAS_STATE)
	{
		T = Press(state) / (Dens(state)*R(state)) ;
	}
	else
	{
		float c_small[7];
		float c_large[7];
		float y, px;
		int i;
		float E = internal_energy(state);

		computeCoeffs(state, c_small, c_large);
		E = E/Dens(state);

		//T = IDEAL_EOS_bis_rootFind(state, c_small, c_large, E, 200, 3500, 100);
		//T = IDEAL_EOS_newt_rootFind(state, c_small, c_large, E, T, 1e-6);
			
		IDEAL_RT_PRMS RTprms; // root parameters

		uni_array(&RTprms.c_small,7,sizeof(float));
		uni_array(&RTprms.c_large,7,sizeof(float));

		for(i=0; i<7; i++)
		{
			RTprms.c_small[i] = c_small[i];
			RTprms.c_large[i] = c_large[i];
		}
	
		RTprms.E = E; 
		RTprms.R = R(state);
	
		y = 0.0;
		find_root(ideal_Teq, (POINTER)&RTprms, y, &px, 3500.0, 200, 1.0e-2, 1.0e-2); 	
		T = px;
	
		free(RTprms.c_small);
		free(RTprms.c_large);
	}
	return T;
}		/*end IDEAL_temperature*/

/*
*			IDEAL_entropy():
*
*		BASIC TEMPERATURE DEPENDENT FUNCTION
*
*	Returns the specific entropy of a state.
*/

LOCAL	float	IDEAL_entropy(
	Locstate state)
{
	if (state_type(state) == VGAS_STATE)
            return Entropy(state);
#if !defined(UNRESTRICTED_THERMODYNAMICS)
        if (Dens(state) < Vacuum_dens(state))
            return 0.0;
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */

        //return (log(pressure(state)) - Gamma(state)*log(Dens(state))) *
        return (log(pressure(state)) - (Gamma(state)*log(Dens(state))))*Coef6(state)*R(state);
}		/*end IDEAL_entropy*/


/***************PRIMARY THERMODYNAMIC FUNCTIONS ****************************/

/*
*			IDEAL_internal_energy():
*
*	Returns the internal energy per unit volume of a state.
*/

LOCAL	float	IDEAL_internal_energy(
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
	    return 	Dens(state) * Int_en(state);
		
	case TGAS_STATE:
	    return	Press(state) / (Gamma(state) - 1.0);

	case FGAS_STATE:
	    return 	( R(state)*Dens(state)*Temperature(state) )/(Gamma(state) - 1.0);

	default:
	    screen("ERROR: in IDEAL_internal_energy(), no such state type\n");
	    clean_up(ERROR);
	}
	return ERROR_FLOAT;
}		/*end IDEAL_internal_energy*/




/*
*		IDEAL_acoustic_impedance_squared():
*
*	Returns the square of the local acoustic impedance of the state.
*
*                        2     dP  |
*			i = - ---- |
*                              dV  |S
*/

LOCAL	float	IDEAL_acoustic_impedance_squared(
	Locstate state)
{
        if (state_type(state) == VGAS_STATE)
        {
            float i = Dens(state)*Sound_speed(state);
            return i*i;
        }
	else
		return Gamma(state)*pressure(state)*Dens(state);
}		/*end IDEAL_acoustic_impedance_squared*/






/*
*			IDEAL_adiabatic_gamma():
*
*	Returns the dimensionless sound speed
*
*		gamma = - d(log P)/d(log V) | .
*					     S
*	As usual P = thermodynamic pressure,  V = specific volume
*	and S = specific entropy.
*/

LOCAL	float	IDEAL_adiabatic_gamma(
	Locstate state)
{
	screen("ERROR in IDEAL_adiabatic_gamma(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end IDEAL_adiabatic_gamma*/





/************** MATERIAL PROPERTY FUNCTIONS ********************************/
/*ARGSUSED*/
LOCAL	float	IDEAL_shear_viscosity(
	Locstate state)
{
	return 0.0;
}	/*end IDEAL_shear_viscosity */

/*ARGSUSED*/
LOCAL	float	IDEAL_bulk_viscosity(
	Locstate state)
{
	return 0.0;
}	/*end IDEAL_bulk_viscosity */

/*ARGSUSED*/
LOCAL   float   IDEAL_heat_coeff(
        Locstate state)
{
	return 0.0;
}               /*end IDEAL_heat_coeff*/
/************** END MATERIAL PROPERTY FUNCTIONS ****************************/

/************** VECTORIZED THERMODYNAMIC FUNCTIONS *************************/

/*
*		IDEAL_single_eos_load_pressure_and_sound_speed2():
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

LOCAL	void	IDEAL_single_eos_load_pressure_and_sound_speed2(
	Vec_Gas *vst,
	int offset,
	int vsize)
{
	static Locstate tmpst = NULL;
	Locstate        *state = vst->state + offset;
	float           *rho = vst->rho + offset;
	float           *p = vst->p + offset, *c2 = vst->c2 + offset;
	float           *e = vst->e + offset;
	float           *FD = NULL;
	int             k;

	if (vst->FD != NULL)
	    FD = vst->FD + offset;

	if (tmpst == NULL)
	{
	    (*Params(state[0])->_alloc_state)(&tmpst,Params(state[0])->sizest);
	    zero_scalar(tmpst,Params(state[0])->sizest);
	    set_type_of_state(tmpst,EGAS_STATE);
	}
	for (k = 0; k < vsize; ++k)
	{
	    Dens(tmpst) = rho[k];
	    Energy(tmpst) = e[k];
	    Set_params(tmpst,state[k]);
	    p[k] = pressure(tmpst);
	    c2[k] = sound_speed_squared(tmpst);
	    if (FD != NULL)
	        FD[k] = fundamental_derivative(tmpst);
	}
}		/*end IDEAL_single_eos_load_pressure_and_sound_speed2*/

/*
*		IDEAL_single_eos_load_pressure_and_gammas():
*
*	Loads the pressure, adiabatic exponent, and Gruneisen
*	coefficient uni_arrays of the Vec_Gas state vst.
*	This function assumes that the specific internal energy
*	uni_array vst->e is already loaded.
*
*	NOTE :
*	Only callable via the function wrapper load_pressure_and_gammas.
*	Assumes that the specific internal energy field is set.
*	This function could be written in terms of the locstate
*	thermodynamic functions,  but is provided in primitive
*	form for increased efficiency of execution time code.
*/
LOCAL	void	IDEAL_single_eos_load_pressure_and_gammas(
	Vec_Gas *vst,
	int offset,
	int vsize)
{
	float *rho = vst->rho + offset;
	float *p = vst->p + offset;
	float *c2 = vst->c2 + offset, *GAM = vst->GAM + offset;
	Gas_param *params = Params(vst->state[offset]);
	float *gm = Vec_Gas_Gamma(vst,offset,vsize);
	int   k;

	if (Vec_Gas_field_set(vst,re))
	{
	    float *re = vst->re + offset;
	    for (k = 0; k < vsize; ++k)
	    {
	        GAM[k] = gm[k] - 1.0;
		//p[k] = pressure(vst->state[offset+k]);
	        p[k] = (gm[k]-1.0)*re[k];
	    }
	}
	else
	{
	    float *e = vst->e + offset;
	    for (k = 0; k < vsize; ++k)
	    {
	        GAM[k] = gm[k] - 1.0;
		//p[k] = pressure(vst->state[offset+k]);
	        p[k] = (gm[k]-1.0)*rho[k]*e[k];
	    }
	}
	if (vst->FD != NULL)
	{
	    float *FD = vst->FD + offset;
	    for (k = 0; k < vsize; ++k)
	        FD[k] = 0.5*(gm[k]+1.0);
	}
#if !defined(UNRESTRICTED_THERMODYNAMICS)
	limit_pressure(p,vst->min_pressure + offset,vsize);
#endif /* !defined(UNRESTRICTED_THERMODYNAMICS) */
	for (k = 0; k < vsize; ++k)
	    c2[k] = gm[k]*p[k]/rho[k];
}		/*end IDEAL_single_eos_load_pressure_and_gammas*/


/*
*			IDEAL_single_eos_load_pressure():
*
*	Loads a vector of pressures into the appropriate field of the 
*	Vec_Gas structure.
*
*	NOTE :
*	Only callable via the function wrapper load_pressure.
*	Assumes that the specific internal energy field is set.
*	This function could be written in terms of the locstate
*	thermodynamic functions,  but is provided in primitive
*	form for increased efficiency of execution time code.
*/

/*ARGSUSED*/
LOCAL	void	IDEAL_single_eos_load_pressure(
	Vec_Gas *vst,
	int offset,
	int vsize)
{
	static Locstate tmpst = NULL;
	Locstate *state = vst->state + offset;
	float *rho = vst->rho + offset;
	float *p = vst->p + offset;
	float *e = vst->e + offset;
	int i, k;

	float *gm = Vec_Gas_Gamma(vst, offset, vsize);

	Gas_param *params = Params(vst->state[offset]);
	int nc = params->n_comps;
	float **rho0 = vst->rho0;
	if (tmpst == NULL)
	{
	    (*Params(state[0])->_alloc_state)(&tmpst,Params(state[0])->sizest);
	    zero_scalar(tmpst,Params(state[0])->sizest);
	    set_type_of_state(tmpst,EGAS_STATE);
	}
	for (k = 0; k < vsize; ++k)
	{
	    /*Dens(tmpst) = rho[k];
	    Energy(tmpst) = e[k];
	    for (i = 0; i < nc; ++i)
		pdens(tmpst)[i] = pdens(state[k])[i];
	    Set_params(tmpst,state[k]);
	    p[k] = pressure(tmpst);*/
	    //p[k] = pressure(state[k]);
	    p[k] = (gm[k] - 1.0) * rho[k] * e[k];

	}
}		/*end IDEAL_single_eos_load_pressure*/

/*
*		IDEAL_single_eos_load_sound_speed2():
*
*	Loads a vector of sound speeds into the	appropriate fields of the
*	Vec_Gas structure.
*
*	NOTE :
*	Only callable via the function wrapper load_sound_speed.
*	Assumes that the specific internal energy field is set.
*	This function could be written in terms of the locstate
*	thermodynamic functions,  but is provided in primitive
*	form for increased efficiency of execution time code.
*/

LOCAL	void	IDEAL_single_eos_load_sound_speed2(
	Vec_Gas *vst,
	int offset,
	int vsize)
{
	screen("ERROR in IDEAL_single_eos_load_sound_speed2(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return;
}		/*end IDEAL_single_eos_load_sound_speed2*/

/***************END VECTORIZED THERMODYNAMIC FUNCTIONS *********************/

/***************RIEMANN SOLUTIONS UTILITY FUNCTIONS ************************/

/***************Purely Thermodynamic Hugoniot Functions*********************/

/*
*			IDEAL_dens_Hugoniot():
*
*	Given the state state0 on one side of an oblique shock and the pressure
*	p1 on the other side, this function returns the density rho1 of the
*	state with pressure p1.  Rho1 is found by solving the Hugoniot relation
*
*		(p1 + p0)*(1/rho0 - 1/rho1) = 2*(e1 - e0)
*
*	where e0 and e1 are the specific internal energies of the two
*	respective states.  For a given equation of state the specific
*	internal energy can be expressed as a function of the
*	pressure and density.  Thus the above equation can be solved to
*	give rho1 as a function of state0 and p1.
*
*
*	Reference: Courant and Friedrichs page 302 ff.
*/


LOCAL	float	IDEAL_dens_Hugoniot(
	float    p1,
	Locstate st0)
{
	float   p0, c4;

        p0 = pressure(st0);
        c4 = Coef4(st0);
        return Dens(st0)*(p1 + p0*c4)/ (p0 + p1*c4);
}	    /*end IDEAL_dens_Hugoniot*/


/*
*			IDEAL_state_w_pr_on_Hugoniot():
*
*	Given the state state0 on one side of an oblique shock and the pressure
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
*	give rho1 and e1 as a function of state0 and p1.  The internal
*	energy is then given by E1 = r1 * e1.
*
*	IMPORTANT NOTE:
*		If stype1 == GAS_STATE the energy in state1 is
*		the internal energy.  The kinetic energy must
*		be added separately.  The reason for this is
*		that this function is a purely theromdynamic
*		function and is independent of the state
*		velocities.
*
*	Reference: Courant and Friedrichs page 302 ff.
*/

LOCAL	void	IDEAL_state_w_pr_on_Hugoniot(
	Locstate state0,
	float p1,
	Locstate state1,
	int stype1)
{
	screen("ERROR in IDEAL_state_w_pr_on_Hugoniot(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return;
}		/*end IDEAL_state_w_pr_on_Hugoniot*/

/*
*			IDEAL_state_w_mf_sqr_on_Hugoniot():
*
*	Given the state state0 on one side of an oblique shock and the square
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
*		If stype1 == GAS_STATE the energy in state1 is
*		the internal energy.  The kinetic energy must
*		be added separately.  The reason for this is
*		that this function is a purely theromdynamic
*		function and is independent of the state
*		velocities.
*
*/

/*ARGSUSED*/
LOCAL	bool	IDEAL_state_w_mf_sqr_on_Hugoniot(
	Locstate state0,
	float m2,
	Locstate ans,
	int stype1)
{
	screen("ERROR in IDEAL_state_w_mf_sqr_on_Hugoniot(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end IDEAL_state_w_mf_sqr_on_Hugoniot*/


/***************End Purely Thermodynamic Hugoniot Functions*****************/
/***************Velocity Related Hugoniot Functions*************************/

/*
*			IDEAL_pr_normal_vel_wave_curve():
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

LOCAL	float	IDEAL_pr_normal_vel_wave_curve(
	float du,	/* normal velocity change across shock = (u1 - u0)*/
	Locstate state0)
{
	screen("ERROR in IDEAL_pr_normal_vel_wave_curve(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end IDEAL_pr_normal_vel_wave_curve*/

LOCAL	bool	pr_normal_vel_wave_curve_aux(
	float		p,
	float		*du,
	Locstate	state0)
{
	screen("ERROR in one_over_rho_c(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end pr_normal_vel_wave_curve_aux*/


/***************End Velocity Related Hugoniot Functions*********************/
/***************Purely Thermodynamic Adiabatic Wave Curve Functions*********/


/*	
*			IDEAL_dens_rarefaction():
*
*	Given the state state0 and the pressure on the other side  
*	of a simple wave in steady irrotational flow, this function 
* 	returns the density on the other side.
*
*	The answer is give by the solution of the ordinary differential
*	equations 
*
*                 dh    |
*	   	------- |    = V,  h(p0) = h0;
*		  dp    | s
*       and     
*
*		d (rho)     |           1 
*            ------------   |    =   -------,  rho(p0) = rho0; 
*               d p         | s        c^2
*
*       where h is specific enthalpy,  p is pressure, V is specific 
*	volume (V = 1/rho), rho is density, s is entropy, and c is 
*	sound speed. The derivatives are taken at constant entropy. 
*	
*	Two equations are necessary, because c^2 depends on p which is
*	computed by the EOS-specific pressure fcn using e (energy) which is
*	obtained from h, which is the solution of the first equation.
*/

/*ARGSUSED*/
LOCAL	bool	rndr(
	double p,
	double *y,
	double *f,
	int n,
	Locstate st)
{
	screen("ERROR in rndr(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/* end rndr*/

LOCAL	float	IDEAL_dens_rarefaction(
	float p1,
	Locstate st0)
{
	return Dens(st0) * pow(p1/pressure(st0),1.0/Gamma(st0));
}		/*end IDEAL_dens_rarefaction*/

/*	
*			IDEAL_pressure_rarefaction():
*
*	Given the state state0 and the density on the other side of
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

/*ARGSUSED*/
LOCAL	bool	rnpr(
	double rho,
	double *y,
	double *f,
	int n,
	Locstate st)
{
	screen("ERROR in rnpr(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/* end rnpr*/

LOCAL	float	IDEAL_pressure_rarefaction(
	float rho1,
	Locstate state0)
{
	screen("ERROR in IDEAL_pressure_rarefaction(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end IDEAL_pressure_rarefaction*/


/*	
*			IDEAL_state_on_adiabat_with_pr():
*
*	Given the state state0 and the pressure on the other side of
*	a simple wave in steady irrotational flow, this function returns
*	the thermodynamic variable on the other side.
*
*	IMPORTANT NOTE:
*		If stype1 == GAS_STATE the energy in state1 is
*		the internal energy.  The kinetic energy must
*		be added separately.  The reason for this is
*		that this function is a purely theromdynamic
*		function and is independent of the state
*		velocities.
*
*/

LOCAL	void	IDEAL_state_on_adiabat_with_pr(
	Locstate st0,
	float p1,
	Locstate st1,
	int stype1)
{
	int     k;

        zero_state_velocity(st1,Params(st0)->dim);
        Set_params(st1,st0);
        set_type_of_state(st1,stype1);
        Dens(st1) = Dens(st0)*pow(p1/pressure(st0),1.0/Gamma(st0));
	
	//----NEW
	Locstate tmpst;
	size_t sizest = Params(st0)->sizest;
	g_alloc_state(&tmpst,sizest);
	zero_state_velocity(tmpst,Params(st0)->dim);
	Set_params(tmpst,st1);
	set_type_of_state(tmpst,TGAS_STATE);
	Dens(tmpst) = Dens(st0)*pow(p1/pressure(st0),1.0/Gamma(st0));
	Press(tmpst) = p1;
        //----
	
	if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            for(k = 0; k < Params(st0)->n_comps; k++)
            {
                pdens(st1)[k] = (pdens(st0)[k]/Dens(st0))*Dens(st1);
                /* New 051005 */
                if(fabs(pdens(st1)[k]) < 10.0*MACH_EPS && pdens(st1)[k] < 0.0)
                    pdens(st1)[k] = 0.0;
                /* End of New 051005 */
		
		pdens(tmpst)[k] = pdens(st1)[k];
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
                //Energy(st1) = p1*Coef6(st1)/Dens(st1);
		//----NEW 3/4/2012
		Energy(st1) = ( p1 / (Gamma(tmpst) -1.0) ) / Dens(st1);

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

	free(tmpst);
}		/*end IDEAL_state_on_adiabat_with_pr*/

/*	
*			IDEAL_state_on_adiabat_with_dens():
*
*	Given the state state0 and the density on the other side of
*	a simple wave in steady irrotational flow, this	function returns
*	the pressure and internal energy on the other side.
*
*	IMPORTANT NOTES:
*		1.  If stype1 == GAS_STATE the energy in state1 is
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

LOCAL	void	IDEAL_state_on_adiabat_with_dens(
	Locstate state0,
	float rho1,
	Locstate state1,
	int stype1)
{
	screen("ERROR in IDEAL_state_on_adiabat_with_dens(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return;
}		/*end IDEAL_state_on_adiabat_with_dens*/




/***************End Purely Thermodynamic Adiabatic Wave Curve Functions*****/
/***************General Wave Curve Functions********************************/

/*
*			IDEAL_mass_flux():
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

LOCAL	float	IDEAL_mass_flux(
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

}		/*end IDEAL_mass_flux*/

/*
*			IDEAL_mass_flux_squared():
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

LOCAL	float	IDEAL_mass_flux_squared(
	float p,
	Locstate state0)
{
	screen("ERROR in IDEAL_mass_flux_squared(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end IDEAL_mass_flux_squared*/


/***************End General Wave Curve Functions****************************/
/***************Functions for the Evaluation of Riemann Solutions***********/

/*
*			IDEAL_oned_fan_state():
*
*	This is a utility function provided for the evaluation of states
*	in a simple wave.   Given sta, it solves for stm using the
*	equation:
*
*	                     / p_m        |            / p_m        |
*	                    /             |           /             |
*	                    \       dP    |           \       G dP  |
*	    w = c_m - c_a +  \    -----   |         =  \     ------ |
*	                      \   rho c   |             \    rho c  |
*	                      /           |             /           |
*	                     /p_a         | S = S_a    / p_a        | S = S_a
*
*	                                               / c_m        |
*	                                              /             |
*	                                              \        dc   |
*	                                            =  \     ------ |
*	                                                \     mu^2  |
*	                                                /           |
*	                                               / c_a        | S = S_a
*
*
*	here c is the sound speed,  rho the density,  S the specific entropy,
*	p the pressure,  and mu^2 = (G - 1)/G,  where G is the fundamental
*	derivative of gas dynamics.  The returned state1 contains only
*	the thermodyanics of the state in the rarefaction fan.  In particular
*	state1 can be used to evaluate the pressure, density, and sound speed
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

struct FAN_AUX	{
	Locstate sta;
	Locstate stm;
	float	ca, cm;
};

LOCAL	float	IDEAL_oned_fan_state(
	float    w,
	Locstate sta,
	Locstate stb,
	Locstate stm,
	int      stype_m,
	bool     *vacuum)
{
	float   c_a, c_m, p_a;
        float   c2, c3, c4;
        int     k;

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
}		/* end IDEAL_oned_fan_state*/

LOCAL	bool	oned_fan_aux(
	float	p,
	float	*w,
	POINTER	prms)
{
	screen("ERROR in oned_fan_aux(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end oned_fan_aux*/
/***************End Functions for the Evaluation of Riemann Solutions********/



/***************Functions to Compute Riemann Solutions**********************/


/*
*			IDEAL_riemann_wave_curve():
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
LOCAL   float   IDEAL_riemann_wave_curve(
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
}               /*end IDEAL_riemann_wave_curve*/

/*
*		IDEAL_set_state_for_find_mid_state():
*
*	Copies the Gas state st into the thermodynamic
*	version Tst, for some EOSs a VGas state is set.
*
*	Technical function added for enhanced performance.
*/

LOCAL	void	IDEAL_set_state_for_find_mid_state(
	Locstate Tst,
	Locstate st)
{
	set_state(Tst,TGAS_STATE,st);
	//set_state(Tst,VGAS_STATE,st);
}		/*end IDEAL_set_state_for_find_mid_state*/

/*
*			IDEAL_eps_for_Godunov():
*
*	Returns a tolerance to be used to determine convergence of the
*	of the Riemann solver.
*
*	Technical function added for enhanced performance.
*/

/*ARGSUSED*/
LOCAL	float	IDEAL_eps_for_Godunov(
	Locstate state,
	float pstar,
	float r_eps)
{
	return r_eps;
}		/*end IDEAL_eps_for_Godunov*/


/*
*			IDEAL_initialize_riemann_solver()
*
*	Computes the epsilons and the initial guess for the pressure
*	in the secant iteration of find_mid_state.
*
*	Technical function added for enhanced performance.
*/

/*ARGSUSED*/
LOCAL   void    IDEAL_initialize_riemann_solver(
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
}               /*end IDEAL_initialize_riemann_solver*/


/***************End Functions to Compute Riemann Solutions******************/
/***************END RIEMANN SOLUTIONS UTILITY FUNCTIONS ********************/

/***************TWO DIMENSIONAL RIEMANN SOLUTION UTILTITY FUNCTIONS*********/

/*
*			IDEAL_steady_state_wave_curve():
*
*	Calculates the steady state wave curve through the state
*	state0 with steady state flow speed q0,  parameterized
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
*	where gamma0 = adiabatic_gamma(state0), Mn = m/(rho0*c0) =
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

LOCAL	bool	IDEAL_steady_state_wave_curve(
	float p1,
	float M0sq,
	float *theta,
	Locstate state0)
{
	screen("ERROR in IDEAL_steady_state_wave_curve(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end IDEAL_steady_state_wave_curve*/


/*
*			IDEAL_pressure_at_sonic_point():
*
*	Returns the pressure at the sonic point of the shock polar
*	through the state state0 with steady state Mach number M0.
*/

LOCAL	float	IDEAL_pressure_at_sonic_point(
	float M0sq,
	Locstate state0)
{
	screen("ERROR in IDEAL_pressure_at_sonic_point(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end IDEAL_pressure_at_sonic_point*/

LOCAL	bool	Mach_number_squared_behind_oblique_shock(
	float		p,
	float		*M2,
	Locstate	st0)
{
	screen("ERROR in Mach_number_squared_behind_oblique_shock(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}

/*
*			IDEAL_pr_at_max_turn_angle():
*
*	Given state0 and the Mach number (squared) of state0 in the frame
*	of a shock, this function calculates the pressure at the point of
*	maximum turning angle on the turning angle pressure shock polar
*	through state0.
*
*	Returns FUNCTION_SUCCEEDED if sucessful,  FUNCTION_FAILED otherwise.
*/

LOCAL	bool	IDEAL_pr_at_max_turn_angle(
	float *prm,
	float M0sq,	/* Mach number of state0 in the frame of the shock */
	Locstate state0)
{
	screen("ERROR in IDEAL_pr_at_max_turn_angle(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end IDEAL_pr_at_max_turn_angle*/

LOCAL	bool	mta_aux(
	float		p,
	float		*y,
	Locstate	st0)
{
	screen("ERROR in mta_aux(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end mta_aux*/


/*
*		IDEAL_state_in_prandtl_meyer_wave():
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
*	The returned state1 contains only
*	the thermodyanics of the state in the rarefaction fan.  In particular
*	state1 can be used to evaluate the pressure, density, and sound speed
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
LOCAL	float	IDEAL_state_in_prandtl_meyer_wave(
	float w,
	float A_a,
	Locstate sta,
	float A_b,
	Locstate stb,
	Locstate stm,
	int stype_m)
{
	screen("ERROR in IDEAL_state_in_prandtl_meyer_wave(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/* end IDEAL_state_in_prandtl_meyer_wave*/

LOCAL	bool	pmw_aux(
	float		p,
	float		*y,
	Locstate	st0)
{
	screen("ERROR in pmw_aux(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end pmw_aux*/

/***************END TWO DIMENSIONAL RIEMANN SOLUTION UTILTITY FUNCTIONS*****/

#if defined(COMBUSTION_CODE)
/***************DETONATION SPECIFIC UTILITY FUNCTIONS*********************/

/*
*			IDEAL_CJ_state():
*
* 	This routine finds the state behind a CJ_detonation.
*	The inputs are the initial state "start"
*	and the side (l_or_r, -1 or 1) we are on.
*/

/*ARGSUSED*/
LOCAL	float	IDEAL_CJ_state(
	Locstate CJ,
	int st_type_CJ,
	Locstate start,
	int l_or_r,
	int avail)
{
	screen("ERROR in IDEAL_CJ_state(), function unimplemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}	/* end IDEAL_CJ_state*/


/*
*	 		IDEAL_progress_state(): 
*
*	Finds the gas state as a function of reaction progress
*	in the steady frame.  
*/
	
/*ARGSUSED*/
LOCAL	void	IDEAL_progress_state(
	float prog,		/* reaction progress */
	Locstate init,		/* init = state behind front */
	Locstate ans,		/* TGas states, init = state behind front */
	float max_vol)		/* maximum allowed volume of reacted state */
{
	screen("ERROR in IDEAL_progress_state(), function unimplemented\n");
	clean_up(ERROR);
}	/* end IDEAL_progress_state*/

/*
*			IDEAL_fprint_combustion_params():
*
*	Prints combustion related parameters.
*/

/*ARGSUSED*/
LOCAL	void	IDEAL_fprint_combustion_params(
	FILE *file,
	Gas_param *params)
{
	screen("ERROR in IDEAL_fprint_combustion_params(), "
	       "function unimplemented\n");
	clean_up(ERROR);
}
/***************END DETONATION SPECIFIC UTILITY FUNCTIONS*****************/
#endif /* defined(COMBUSTION_CODE) */


/***************METHOD OF CHARACTERISTIC FUNCTIONS FOR W_SPEED**************/

/*
*			IDEAL_neumann_riem_inv_moc():
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
*		                            = nor[0] for cylindrical symmetry
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

/*ARGSUSED*/
LOCAL	void	IDEAL_neumann_riem_inv_moc(
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
	screen("ERROR in IDEAL_neumann_riem_inv_moc(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	//return ERROR_FLOAT;
	return;
}		/*end IDEAL_neumann_riem_inv_moc*/

/*
*		IDEAL_shock_ahead_state_riem_inv_moc():
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
*		                   /ans(position = pt + W*dt)
*		                  /  + 0 -   
*		                 /    +  0  -   
*		                /      +   0   -   
*		          shock/        +    0    -   
*		         front/          +     0     -   
*		             /            +      0      -   
*			____/st0__________st3____st2____st1____
*		           pt
*
*		side = NEGATIVE_SIDE
*	                                  \
*			                ans\(position = pt + W*dt)
*			             + 0 -  \
*			          +  0  -    \shock
*			       +   0   -      \front
*			    +    0    -        \
*		         +     0     -          \
*                     +      0      -            \
*		____st3____st2____st1__________st0\______
*	                                         pt
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
LOCAL	void	IDEAL_shock_ahead_state_riem_inv_moc(
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
	screen("ERROR in IDEAL_shock_ahead_state_riem_inv_moc(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return;
}		/*end IDEAL_shock_ahead_state_riem_inv_moc*/


/*
*			IDEAL_shock_moc_plus_rh():
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
*		nor	shock normal
*		W	first prediction of shock front velocity
*		front   Front structure
*	Output variables:
*		ans	times updated state behind shock
*		W	updated shock speed
*
*		The answer state ans is returned in Gas format.
*/

/*ARGSUSED*/
LOCAL	bool	IDEAL_shock_moc_plus_rh(
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
	screen("ERROR in IDEAL_shock_moc_plus_rh(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end IDEAL_shock_moc_plus_rh*/

/***************END METHOD OF CHARACTERISTIC FUNCTIONS FOR W_SPEED**********/

/***************INITIALIZATION UTILITY FUNCTIONS****************************/

/*
*			IDEAL_prompt_for_state():
*
*	Prompts for a hydrodynamical state.  The form of
*	the states depends of the Eos. 	The type of the state
*	is returned.
*/

LOCAL	void	IDEAL_prompt_for_state(
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
	set_type_of_state(state,FGAS_STATE);
	Params(state) = params;
	screen("Enter the density, temperature");
	for (i = 0; i < dim; ++i)
	{
	    screen(", ");
	    if (i == (dim - 1))
		screen("and ");
	    screen("%s",velmesg[i]);
	}
	screen("%s: ",mesg);
	(void) Scanf("%f %f",&Dens(state),&Temperature(state));
	for (i = 0; i < dim; ++i)
	    (void) Scanf("%f",&Vel(state)[i]);
	(void) getc(stdin); /*read trailing newline*/

	set_state(state,stype,state);
}		/* end IDEAL_prompt_for_state */

/*
*			IDEAL_prompt_for_thermodynamics():
*
*	Prompts for the thermodynamic variables in a state.  Returns
*	a state with the appropriate thermodynamic state and zero velocity.
*	The return status gives the state type representation of the state.
*/

LOCAL	void	IDEAL_prompt_for_thermodynamics(
	Locstate   state,
	Gas_param  *params,
	const char *mesg)
{
	if (params == NULL)
	{
	    g_obstacle_state(state,g_sizest());
	    return;
	}
	set_type_of_state(state,FGAS_STATE);
	zero_state_velocity(state,MAXD);
	Params(state) = params;
	screen("Enter the density and temperature%s: ",mesg);
	(void) Scanf("%f %f\n",&Dens(state),&Temperature(state));
}		/* end IDEAL_prompt_for_thermodynamics */


/***************Problem Type Specific Initialization Functions**************/

/*
*			IDEAL_RT_RS_f():
*
*	Support function for the computation of a solution to the linearized
*	Rayleigh-Taylor flow.
*
*	NEEDED:  More complete documentation
*/

/*ARGSUSED*/
LOCAL	float	IDEAL_RT_RS_f(
	float		s,
	Locstate	amb_st,
	float		dz,
	float		k_sqr,
	float		g_z)
{
	screen("ERROR in IDEAL_RT_RS_f(), function unimplemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end IDEAL_RT_RS_f*/


/*
*			IDEAL_RT_single_mode_perturbation_state():
*
*	Computes the perturbation term for the solution to the linearized
*	Euler equations in the Rayleigh-Taylor problem.  See the appendix to
*
*	``The Dynamics of Bubble Growth for
*				Rayleigh-Taylor Unstable Interfaces''
*	C. L. Gardner, J. Glimm, O. McBryan, R. Menikoff, D. H. Sharp,
*	and Q. Zhang, Phys. Fluids 31 (3), 447-465 (1988).
*
*	for an explanation of the formulas for the case where both fluids
*	are stiffened polytropic gases.
*
*	Note that ans is only the perturbation to the ambient state.
*/

/*ARGSUSED*/
LOCAL	void	IDEAL_RT_single_mode_perturbation_state(
	Locstate	ans,
	float		*coords,
	float		t,
	Locstate	amb_st,
	float		z_intfc,
	float		z_bdry,
	MODE		*mode,
	float		g_z)
{
	screen("ERROR in IDEAL_RT_single_mode_perturbation_state(), "
	       "function unimplemented\n");
	clean_up(ERROR);
}		/*end IDEAL_RT_single_mode_perturbation_state*/

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

/*ARGSUSED*/
LOCAL	void	IDEAL_KH_single_mode_state(
	Locstate	ans,
	float		*coords,
	float		t,
	Locstate	amb_st,
	float		stream_velocity,
	float		z_intfc,
	float		z_bdry,
	MODE		*mode)
{
	screen("ERROR in IDEAL_KH_single_mode_state(), "
	       "function unimplemented\n");
	clean_up(ERROR);
}		/*end IDEAL_KH_single_mode_state*/


/*
*		IDEAL_compute_isothermal_stratified_state():
*
*	Solves for the state at height dz above the reference state
*	ref_st in an isothermal one dimensional steady flow.
*
*	The solution is computed by solving the differential
*	equation:
*
*	    P_z = rho gz,	P(0) = P_R, rho(0) = rho_r, T = T_r.
*
*	The solution is given implicitly by
*
*	    f = f_R + g*dz
*	    T = T_r
*
*	where f = e + P/rho - T*S.
*/

/*ARGSUSED*/
LOCAL	void	IDEAL_compute_isothermal_stratified_state(
	Locstate	ans,
	float		dz,	/* distance from reference position */
	float		gz,	/* gravity */
	Locstate	ref_st)
{
	screen("ERROR in IDEAL_compute_isothermal_stratified_state(), "
	       "function unimplemented\n");
	clean_up(ERROR);
}		/*end IDEAL_compute_isothermal_stratified_state */

/*
*		IDEAL_compute_isentropic_stratified_state():
*
*	Solves for the state at height dz above the reference state
*	ref_st in an isentropic one dimensional steady flow.
*
*	The solution is computed by solving the differential equation:
*
*	    P_z = rho gz,	P(0) = P_R, rho(0) = rho_r, S = S_r.
*
*	The solution is given implicitly by
*
*	    h =  h_r + g*gz
*	    S = S_r
*
*	where h = e + P/rho is the specific enthalpy.
*
*/

/*ARGSUSED*/
LOCAL	void	IDEAL_compute_isentropic_stratified_state(
	Locstate	ans,
	float		dz,	/* distance from reference position */
	float		gz,	/* gravity */
	Locstate	ref_st)
{
	screen("ERROR in IDEAL_compute_isentropic_stratified_state(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return;
}	/*end IDEAL_compute_isentropic_stratified_state*/

/*
*		IDEAL_compute_constant_density_stratified_state():
*
*	Solves for the state at height dz above the reference state
*	ref_st in a constant density one dimensional steady flow.
*
*	The solution is computed by solving the differential
*	equation:
*
*		P_z = rho gz,	P(0) = P_R
*/

LOCAL	void	IDEAL_compute_constant_density_stratified_state(
	Locstate	ans,
	float		dz,	/* distance from reference position */
	float		gz,	/* gravity */
	Locstate	ref_st)
{
	screen("ERROR in IDEAL_compute_constant_density_stratified_state(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return;
}	/*end IDEAL_compute_constant_density_stratified_state*/

/***************End Problem Type Specific Initialization Functions**********/
/***************END INITIALIZATION UTILITY FUNCTIONS************************/

/***************EQUATION OF STATE DOMAIN FUNCTIONS**************************/

LOCAL	float	IDEAL_Min_energy(
	Locstate	state)
{
	return Params(state)->min_energy;
}	/*end IDEAL_Min_energy*/

LOCAL	float	IDEAL_Min_pressure(
	Locstate	state)
{
	return Params(state)->min_pressure;
}	/*end IDEAL_Min_pressure*/

LOCAL	float	IDEAL_Vacuum_dens(
	Locstate	state)
{
	return Params(state)->vacuum_dens;
}	/*end IDEAL_Vacuum_dens*/

LOCAL	float	IDEAL_Raref_press(
	Locstate	state)
{
	return Params(state)->raref_press;
}	/*end IDEAL_Raref_press*/

#if defined(COMBUSTION_CODE)
LOCAL	float	IDEAL_Tol_alpha(
	Locstate	state)
{
	return Params(state)->tol_alpha;
}	/*end IDEAL_Tol_alpha*/

LOCAL	float	IDEAL_Tol_pressure(
	Locstate	state)
{
	return Params(state)->tol_press;
}	/*end IDEAL_tol_pressure*/
#endif /* defined(COMBUSTION_CODE) */

/***************END EQUATION OF STATE DOMAIN FUNCTIONS*************************/

/*
*			prandtl_meyer_speed_rinv():
*
*	This function computes the purely speed dependent part of
*	the Riemann invariant for a steady irrotational two dimensional
*	supersonic flow.  The full riemann invariants are given by
*
*		theta - rinv for the GAMMA+ family
*
*	and
*
*		theta + rinv for the GAMMA- family.
*
*	This function returns rinv as the value of the integral
*
*	/ p0					/ q1
*      /                   2		       /                   2
*      \     sqrt(1 - (c/q)  ) * dp      =     \     sqrt(1 - (c/q)  ) * dq
*	\		        ----            \                       ----
*	/		        q*c*rho         /			 c
*      / p1				       / q0
*
*	Where the sound speed c and the flow speed q are related by
*	the Bernoulli relation
*
*			0.5*q*q + i = 0.5*sqr(q0) + i0
*
*	where i is the specific internal enthalpy and qhat is the
*	constant critical speed.
*
*	It is assumed that state1 is has correctly ft_assigned thermodynamic
*	data so that state1 and state0 have the same entropy and that
*	q0, state0, q1 and state1 satisfy the Bernoulli relation above.
*
*	In addition this function computes and returns the Mach angles 
*	A0, and A1 for the flow given by
*
*			sin(Ai) = qi/ci, i = 0,1.
*
*	The flow speeds of the two states may be passed either directly
*	or as Mach numbers.  If the flag mach_numbers_given is YES
*	the quantities x0sq and x1sq are assumed to be the steady
*	state Mach numbers of the two states respectively.  Otherwise
*	these numbers are assumed to be the steady state flow speeds.
*/

LOCAL	bool	prandtl_meyer_speed_rinv(
	float		*A0,
	float		x0sq,
	Locstate	state0,
	float		*A1,
	float		x1sq,
	Locstate	state1,
	float		*rinv)
{
	screen("ERROR in prandtl_meyer_speed_rinv(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end prandtl_meyer_speed_rinv*/


/*ARGSUSED*/
LOCAL	bool	pmsr(
	double p,
	double *y,
	double *f,
	int n,
	Locstate state0)
{
	screen("ERROR in pmsr(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end pmsr*/


/*
*			int_dp_over_rho_c():
*
*
*	Returns the value of the integral:
*
*		 
*                        / p1        |
*                       /            |
*                       \      dP    |
*                        \   ------  |
*                         \   rho c  |
*                         /          |
*                        / p0        | S
*
*	Also returns the state on the adiabat containing state0 with
*	pressure P.
*
*/

LOCAL	float	int_dp_over_rho_c(
	float		p1,
	Locstate	state0,
	Locstate	state1)
{
	screen("ERROR in int_dp_over_rho_c(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end int_dp_over_rho_c*/

LOCAL   float           gaussian_int(
	int             np,
	float           a,
	float           b,
	Locstate	state)
{
	screen("ERROR in gaussian_int(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end gaussian_int*/

LOCAL   void    legendre_init(
	int     np,
	float*  point,
	float*  weight)
{
	screen("ERROR in legendre_init(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return;
}		/*end legendre_init*/

LOCAL   void    linear_transform(
	int     np,
	float*  point,
	float*  weight,
	float   a,
	float   b)
{
	screen("ERROR in linear_transform(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return;
}		/*end linear_transform*/       

LOCAL   void    exp_transform(
	int     np,
	float*  point,
	float*  weight)
{
	screen("ERROR in exp_transform(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return;
}		/*end exp_transform*/       

LOCAL   float   one_over_rho_c(
	float		p,
	Locstate	state)
{
	screen("ERROR in one_over_rho_c(), "
	       "function unimplemented. Look at generic-eos.c"
		"to see how the function can be implemented\n");
	clean_up(ERROR);
	return ERROR_FLOAT;
}		/*end one_over_rho_c*/


LOCAL   float   IDEAL_specific_enthalpy_species(
        Locstate state,
        int      i)
{
        float coeff[7];
        float r, T;
	int j;
        
        IDEAL_EOS *eos = (IDEAL_EOS *)(Params(state)->eos);
        r = (IDEAL_Eos(state)->R);
        r = r/ (eos->_M[i]);        
        T = temperature(state);

	if( T < eos->T)
	{
		for ( j=0; j<7; j++ )
			coeff[j] = coeff_small(IDEAL_Eos_Component(state,i),j);

	}else{
		for ( j=0; j<7; j++ )
			coeff[j] = coeff_large(IDEAL_Eos_Component(state,i),j);

	}

        return (r*T*(coeff[0]+coeff[1]*T/2+coeff[2]*T*T/3+coeff[3]*T*T*T/4+
                                 coeff[4]*T*T*T*T/5+coeff[5]/T));
}       /* end IDEAL_specific_enthalpy_species */

LOCAL	bool	ideal_Teq(float x, float *fans, POINTER prm)
{
	float   coeff[7];
	float   E, M, R;
	int	i;

	if( x < T_crit)
	{	
		for(i=0; i<7; i++)	
			coeff[i]  = ((IDEAL_RT_PRMS *) prm)->c_small[i];
	}
	else
	{
		for(i=0; i<7; i++)	
			coeff[i]  = ((IDEAL_RT_PRMS *) prm)->c_large[i];
	}

	E  = ((IDEAL_RT_PRMS *) prm)->E;
	R  = ((IDEAL_RT_PRMS *) prm)->R;

	*fans = (coeff[0]+x*(coeff[1]+x*(coeff[2]+x*(coeff[3]+x*coeff[4])))- 1)*R*x - E;

	return FUNCTION_SUCCEEDED;
}

LOCAL   float IDEAL_dynamic_viscosity(
        Locstate state,
        float    T)
{
        int i,j,n_comps = Num_gas_components(state);
        float Mass_frac[MAX_NCOMPS];
        IDEAL_EOS *eos = (IDEAL_EOS *)(Params(state)->eos);

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

        for(i =0; i<n_comps;i++)
        {
            Mass_frac[i] = pdens(state)[i]/Dens(state);
            T_x = T/Tem_lennard[i];
            cl = Collision_intergral_2(T_x);
            pvisc[i] = 0.001*2.6693*0.00001*sqrt((eos->_M[i])*T)/(delta_lennard[i]*delta_lennard[i]*cl);
        }

        for( i = 0; i< n_comps; i++)
        {
            for(j=0;j<n_comps;j++)
            {
                PM[i][j] = eos->_M[i]/eos->_M[j];
                Pmu[i][j] =  pvisc[i]/pvisc[j];
            }
        }

        for(i=0;i<n_comps;i++)
            for(j=0;j<n_comps;j++)
	    {
		double tmp = sqrt(sqrt(PM[j][i]));
		double tmp1 = 1.0 / sqrt(Pmu[i][j]);
		double tmp2 = tmp * tmp1;
		double tmp3 = 1.0 / sqrt(1+PM[i][j]);
		G[i][j] = 0.353553 * sqr(1+tmp2) * tmp3;
	    }
//                G[i][j] = pow(1 + pow(PM[j][i],0.25)*pow(Pmu[i][j],-0.5),2)*pow(1+PM[i][j],-0.5)*pow(2,-1.5);

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

LOCAL   float IDEAL_dynamic_thermal_conductivity(
        Locstate state,
        float    T)
{
        int i,j,n_comps = Num_gas_components(state);
        float Mass_frac[MAX_NCOMPS];
        IDEAL_EOS *eos = (IDEAL_EOS *)(Params(state)->eos);

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

        for(i =0; i<n_comps;i++)
        {
            Mass_frac[i] = pdens(state)[i]/Dens(state);
            T_x = T/Tem_lennard[i];
            cl = Collision_intergral_2(T_x);
            pvisc[i] = 0.001*2.6693*0.00001*sqrt((eos->_M[i])*T)/(delta_lennard[i]*delta_lennard[i]*cl);
        }

        for( i = 0; i< n_comps; i++)
        {
            for(j=0;j<n_comps;j++)
            {
                PM[i][j] = eos->_M[i]/eos->_M[j];
                Pmu[i][j] =  pvisc[i]/pvisc[j];
            }
        }

        for(i=0;i<n_comps;i++)
            for(j=0;j<n_comps;j++)
	    {
		double tmp = sqrt(sqrt(PM[j][i]));
		double tmp1 = 1.0 / sqrt(Pmu[i][j]);
		double tmp2 = tmp * tmp1;
		double tmp3 = 1.0 / sqrt(1+PM[i][j]);
		G[i][j] = 0.353553 * sqr(1+tmp2) * tmp3;
	    }
//                G[i][j] = pow(1 + pow(PM[j][i],0.25)*pow(Pmu[i][j],-0.5),2)*pow(1+PM[i][j],-0.5)*pow(2,-1.5);

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
            Plamda[i] = pvisc[i]*(C_P[i] + (5.0*IDEAL_Eos(state)->R/(4.0*eos->_M[i])));
            lamda += Mass_frac[i]*Plamda[i]/delta[i];
        }

        return lamda;
}

LOCAL   void IDEAL_dynamic_viscosity_thermalconduct(
        Locstate state,
        float    T,
	float 	*MUrans,
	float 	*Krans)
{
        int i,j,n_comps = Num_gas_components(state);
        float Mass_frac[MAX_NCOMPS];
        IDEAL_EOS *eos = (IDEAL_EOS *)(Params(state)->eos);

        /* using flame master's formula to compute the dynamic/shear viscosity(see the techreport of this part). */
        float pvisc[MAX_NCOMPS];
        float Y[MAX_NCOMPS];
        float G[MAX_NCOMPS][MAX_NCOMPS];
        float delta[MAX_NCOMPS];
        float cl,T_x;
        float PM[MAX_NCOMPS][MAX_NCOMPS],Pmu[MAX_NCOMPS][MAX_NCOMPS];
	float mu = 0;
        float lamda = 0;
        float Plamda[MAX_NCOMPS];
        float C_P[MAX_NCOMPS];
	float *M = eos->_M;
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

        i = 0;
        bool pureo2;
      if( 1 - pdens(state)[i]/Dens(state) < 0.001)
      {
        pureo2 = 1;
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
	    {
		double tmp = sqrt(sqrt(PM[j][i]));
		double tmp1 = 1.0 / sqrt(Pmu[i][j]);
		double tmp2 = tmp * tmp1;
		double tmp3 = 1.0 / sqrt(1+PM[i][j]);
		G[i][j] = 0.353553*sqr(1+tmp2) * tmp3;
	    }

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
	*MUrans = mu;

        for(i = 0;i<n_comps;i++)
        {
            C_P[i] = C_P_species(state,T,i);
            Plamda[i] = pvisc[i]*(C_P[i] + (5.0*IDEAL_Eos(state)->R/(4.0*M[i])));
            lamda += Mass_frac[i]*Plamda[i]/delta[i];
        }
	*Krans = lamda;
      }
      else
      {
        pureo2 = 0;
        i = 0;
        j = 0;
            Mass_frac[i] = pdens(state)[i]/Dens(state);
            T_x = T/Tem_lennard[i];
            cl = Collision_intergral_2(T_x);
            pvisc[i] = 0.001*2.6693*0.00001*sqrt(M[i]*T)/(delta_lennard[i]*delta_lennard[i]*cl);

                PM[i][j] = M[i]/M[j];
                Pmu[i][j] =  pvisc[i]/pvisc[j];

		double tmp = sqrt(sqrt(PM[j][i]));
		double tmp1 = 1.0 / sqrt(Pmu[i][j]);
		double tmp2 = tmp * tmp1;
		double tmp3 = 1.0 / sqrt(1+PM[i][j]);
		G[i][j] = 0.353553*sqr(1+tmp2) * tmp3;

            delta[i] = 0;
            delta[i] += G[i][j]*PM[i][j]*Mass_frac[j];

            mu += Mass_frac[i]*pvisc[i]/delta[i];
	    *MUrans = mu;

            C_P[i] = C_P_species(state,T,i);
            Plamda[i] = pvisc[i]*(C_P[i] + (5.0*IDEAL_Eos(state)->R/(4.0*M[i])));
            lamda += Mass_frac[i]*Plamda[i]/delta[i];
	    *Krans = lamda;
      }
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
	
