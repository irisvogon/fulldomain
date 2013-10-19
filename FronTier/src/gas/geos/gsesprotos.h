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
*				gsesprotos.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
 */

#if !defined(_GSESPROTOS_H)
#define _GSESPROTOS_H

#if defined(SESAME_CODE) && defined(TWOD)

#if defined(cray)
#define	cllib		CLLIB
#define	oplib_f		OPLIB_F
#define	oplib_uf	OPLIB_UF
#define	oplib_f		OPLIB_F
#define	oplib_uf	OPLIB_UF
#define	s4get		S4GET
#define	s2get		S2GET
#define	s2eos		S2EOS
#endif /* defined(cray) */

#include <geos/sesame.h>

	/* geos EXPORTED Function Declarations */

/*	geos/gseshyp.c */
IMPORT	const SESAME_TABLE_TYPE *sesame_table_numbers(void);
IMPORT	const char **sesame_table_names(void);
IMPORT	const char **sesame_headers(void);
IMPORT	const size_t *sesame_size_of_states(void);
IMPORT	int	init_ses_hyp_soln_func(Wave*,Front*);
IMPORT	void	init_sesame_hyp_tri_solns(SESAME_EOS*,POINTER,POINTER);
IMPORT	void	set_default_ses_wave_and_front(INIT_DATA*,Front*,Wave*,
					       size_t,bool);
IMPORT	void	set_ses_hooks(Front*,SESAME_TABLE_TYPE);
IMPORT	void	set_ses_inv_hyp_solns(SESAME_EOS*);
IMPORT	void	set_user_hooks_for_sesame(void);

/*	geos/gsesintrp.c */
IMPORT	bool	ses_tri_interpolator(float*,LINEAR_ELEMENT*,TRI_SOLN*,
				     Locstate);
IMPORT	bool	ses_tri_lin_comb_states(float,float,float,float*,Locstate,
					float*,Locstate,float*,
					Locstate,RECT_GRID*,Locstate);
IMPORT	void	ses_lin_comb_states(float,float,float*,Locstate,float*,
				    Locstate,RECT_GRID*,Locstate);
IMPORT	void	ses_quad_interpolator(float*,BILINEAR_ELEMENT*,TRI_SOLN*,
				      Locstate);
IMPORT	void	set_ses_intrp_flag(int,SESAME_TABLE_TYPE);
IMPORT	void	set_ses_intrp_flag_all(SESAME_TABLE_TYPE);
#if defined(PHASE_CODE)

/*	geos/gsesphase.c */
IMPORT	int	get_phase_hyp_state(float,SESAME_EOS*,float*,float*,float*);
IMPORT	int	get_phase_state(float,SESAME_EOS*,PHASE_BDRY*,
				float*,float*,float*);
IMPORT	void	cold_PE_spline(SESAME_EOS*,COLD_CURVE*);
IMPORT	void	get_phase_temp_state(float,Front*,float*,
				     float*,float*,int*,SESAME_EOS*);
IMPORT	void	init_RT_interior_states(Wave*,Front*,SESAME_EOS*,
					COLD_CURVE*,PHASE_BDRY*);
IMPORT	void	init_new_phase_bound(Wave*,Front*,SESAME_EOS*,COLD_CURVE*,
				     PHASE_BDRY*);
IMPORT	void	lookspl(SESAME_EOS*,PHASE_BDRY*,float,int,int,int,int,
			float*,float*,float*,float*);
IMPORT	void	phase_spline(SESAME_EOS*,COLD_CURVE*,PHASE_BDRY*);
IMPORT	void	ses_phase_states(Front*,SESAME_EOS*,PHASE_BDRY*);
IMPORT	void	set_boundary_states(Front*,SESAME_EOS*,PHASE_BDRY*,
				    COLD_CURVE*);

/*	geos/gsesspline.c */
IMPORT	void	get_temp_state(float,int,SESAME_EOS*,PHASE_BDRY*,float*,float*,
			       float*,int*);
IMPORT	void	init_PE_phase_grid(SESAME_EOS*,float,PHASE_BDRY*,int*);
IMPORT	void	init_PE_spline(SESAME_EOS*,PHASE_BDRY*);
IMPORT	void	set_cross_states(Front*,SESAME_EOS*,PHASE_BDRY*,float*,
				 float*,float*,float*,float*);
#endif /* defined(PHASE_CODE) */

/*	geos/gsesprint.c */
IMPORT	void	fprint_SESAME_params(FILE*,SESAME_EOS*);
IMPORT	void	print_ps_tri_soln(FILE*,SESAME_EOS*);
IMPORT	void	print_re_tri_soln(FILE*,SESAME_EOS*);
IMPORT	void	print_rt_tri_soln(FILE*,SESAME_EOS*);
IMPORT	void	print_rs_tri_soln(FILE*,SESAME_EOS*);
IMPORT	void	print_title_for_sesame(FILE*,SESAME_EOS*);
IMPORT	void	print_vp_tri_soln(FILE*,SESAME_EOS*);
IMPORT	void	read_print_SESAME_params(SESAME_EOS*,const IO_TYPE*);
IMPORT	void	ses_rt_fprint_state_data(FILE*,Locstate,INTERFACE*);
IMPORT	void	ses_re_fprint_state_data(FILE*,Locstate,INTERFACE*);
IMPORT	void	ses_rs_fprint_state_data(FILE*,Locstate,INTERFACE*);
IMPORT	void	ses_ps_fprint_state_data(FILE*,Locstate,INTERFACE*);
IMPORT	void	ses_vp_fprint_state_data(FILE*,Locstate,INTERFACE*);
IMPORT	void	verbose_ses_show_intfc_states(INTERFACE*,SESAME_TABLE_TYPE,
					      SESAME_EOS*);

/*	geos/gsesinout.c */
IMPORT	bool	restart_sesame(INIT_DATA*,SESAME_EOS*);
IMPORT	void	init_ses_prt(Printplot*,Front*,int);
IMPORT	void	ses_printout(CHART*,Printplot*,bool,int);

/*	geos/sesinv.c */
IMPORT	bool	zero_temperature_cold_curve(SESAME_EOS*);
IMPORT	float	sesame_rt_sound_speed_squared(float,float,float*,float*);
IMPORT	void	s2eos_lookup(float,float,float*,float*,SESAME_EOS*);
IMPORT	void	set_RT_entropy_from_cold_curve(Front*,Wave*,SESAME_EOS*);
IMPORT	void	set_RT_entropy_from_mid_point(Front*,Wave*,SESAME_EOS*);
IMPORT	void	setrt(float,float,Locstate,SESAME_EOS*);

/*	geos/sesspln.c */
IMPORT	int	splcomp(float*,float*,float*,int,float*,float*,
				      float*,float*,float*);
IMPORT	int	splcomp2(float*,float*,float*,int,float*,float*,float*,float*,
			 float*,float,float);
IMPORT	void	spline(float,float,float,float,float,float,float,float*,float*);

/*	geos/sesstate.c */
IMPORT	void	phbnd(float*,float,float,float,SESAME_EOS*);
IMPORT	void	sets(float,float*,int,int*,SESAME_EOS*);
IMPORT	void	setspb(float,float*,int,float,SESAME_EOS*);

#if defined(__cplusplus)
extern "C" {
#endif /* defined(__cplusplus) */
/*	geos/gsestoc.F */
    FORTRAN	void	FORTRAN_NAME(cllib)(int*);
#   if defined(cray)
    FORTRAN	void	SFORTRAN_NAME(oplib_f)(int*,_fcd);
    FORTRAN	void	SFORTRAN_NAME(oplib_uf)(int*,_fcd);
#   else /* defined(cray) */
    FORTRAN	void	SFORTRAN_NAME(oplib_f)(int*,const char*,int);
    FORTRAN	void	SFORTRAN_NAME(oplib_uf)(int*,const char*,int);
#   endif /* defined(cray) */

/*	geos/sesadd.F */
    FORTRAN	void	FORTRAN_NAME(s4get)(int*,int*,float*,int*,int*,int*);

/*	geos/sesame.F */
    FORTRAN	void	FORTRAN_NAME(s2get)(int*,int*,float*,int*,int*,int*);
    FORTRAN	void	FORTRAN_NAME(s2eos)(int*,float*,float*,float*,
                                            float*,float*);
#if defined(__cplusplus)
}
#endif /* defined(__cplusplus) */

#endif /* defined(SESAME_CODE) && defined(TWOD) */
#endif /* !defined(_GSESPROTOS_H) */
