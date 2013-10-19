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
*				gnodeprotos.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#if !defined(_GNODEPROTOS_H)
#define _GNODEPROTOS_H

	/* gnode EXPORTED Function Declarations */

#if defined(TWOD)

/*	gnode/gRP.c */
IMPORT	RP_DATA	*allocate_RP_DATA_structure(size_t,bool,int);
IMPORT	void	copy_RP_DATA_structure(RP_DATA*,RP_DATA*);
IMPORT	void	free_RP_DATA_structure(RP_DATA*);

/*	gnode/gbnode.c */
IMPORT	int	B_reflect_node_propagate(Front*,Wave*,NODE*,NODE*,RPROBLEM**,
					 float,float*,NODE_FLAG);
IMPORT	int	attached_b_node_propagate(Front*,Wave*,NODE*,NODE*,
					  RPROBLEM**,float,float*,NODE_FLAG);
IMPORT	int	is_bow_shock_attached(Locstate,float,Locstate,float*);
IMPORT	int	is_regular_reflection(float*,Front*,RP_DATA*);

/*	gnode/gccnode.c */
IMPORT	int	cc_node_propagate(Front*,Wave*,NODE*,NODE*,RPROBLEM**,
				  float,float*,NODE_FLAG);
IMPORT	void	set_cc_normal_function(NORMAL_FUNCTION*);
IMPORT	void	set_cc_tangent_function(TANGENT_FUNCTION*);

/*      gnode/gmdnode.c */
IMPORT	int	reg_to_mach_diff_reconfigure(O_CURVE**,O_CURVE**,NODE*);

/*	gnode/gnode.c */
IMPORT	bool	is_bad_state_at_node(const char*,NODE*);
IMPORT	int	g_node_propagate(Front*,POINTER,NODE*,NODE*,RPROBLEM**,
				 float,float*,NODE_FLAG,POINTER);
IMPORT	void	node_warning(const char*,const char*,const char*,NODE_FLAG);

/*	gnode/gnodesub.c */
IMPORT ADJUST_ANGLE_VALUE g_adjust_angle_len(NODE*,O_CURVE*,O_CURVE*,Locstate,
					     Locstate,float*,float,float*,
					     float*,Front*,Wave*);
IMPORT	float	sonic_radius(Locstate,float*,float,int);
IMPORT	bool	g_check_delete_redundant_node(NODE*,CURVE*,CURVE*);
IMPORT	int	adjust_angle_at_node(NODE*,O_CURVE*,O_CURVE*,Locstate,Locstate,
				     float*,float,float,Front*,Wave*);
IMPORT	ANGLE_DIRECTION	g_find_i_to_prop_dir(Front*,POINTER,NODE*,CURVE*,
				     	     ORIENTATION,float,
				     	     COMPONENT*,POINT*,float*);
IMPORT	int	modify_curves_at_node(POINT*,BOND**,NODE*,O_CURVE**,O_CURVE**,
				      int,bool*,Front*,Wave*,float,
				      RP_DATA*,NODE_FLAG);
IMPORT	int	propagate_curve_near_node(NODE*,NODE*,O_CURVE*,O_CURVE*,BOND*,
					  Locstate,Locstate,bool,float,
					  Front*,Wave*,float,NODE_FLAG);
IMPORT	void	find_adjacent_curves(O_CURVE*,ANGLE_DIRECTION*,O_CURVE*,
				     O_CURVE*,COMPONENT*);
IMPORT	void	find_curve_with_status(NODE*,CURVE**,ORIENTATION*,int);
IMPORT	void	identify_curves_with_status(NODE*,O_CURVE*,O_CURVE*,int);

/*	gnode/gpcnode.c */
IMPORT	int	precursor_shock_rr_propagate(Front*,Wave*,NODE*,NODE*,
					     RPROBLEM**,float,float*,NODE_FLAG,
					     O_CURVE**,O_CURVE**,POINT*,
					     RP_DATA*,BOND**);

/*	gnode/gsndnode.c */
IMPORT	int	g_snd_node_propagate(Front*,Front*,POINTER,INTERFACE*,
				     NODE*,NODE*,float);
/*	gnode/gssnode.c */
IMPORT	int	Mach_node_propagate(Front*,Wave*,NODE*,NODE*,RPROBLEM**,
				    float,float*,NODE_FLAG);
IMPORT	int	cross_node_propagate(Front*,Wave*,NODE*,NODE*,RPROBLEM**,
				     float,float*,NODE_FLAG);
IMPORT	int	overtake_node_propagate(Front*,Wave*,NODE*,NODE*,RPROBLEM**,
					float,float*,NODE_FLAG);
IMPORT	void	ramp_reflection_corner_posn(float*,int,int);
IMPORT	void	tangent_at_degenerate_node(CURVE*,ORIENTATION,
					   CURVE*,ORIENTATION,
					   float*,Front*);

/*	gnode/gssnsts.c */
IMPORT	int	find_Mach_node_states(NODE*,O_CURVE*,O_CURVE*,O_CURVE*,
				      O_CURVE*,O_CURVE*,BOND**,BOND**,POINT*,
				      Wave*,Front*,RPROBLEM**,float,float*,
				      RP_DATA*,NODE_FLAG);
IMPORT	int	find_cross_node_states(float*,RP_DATA*,bool*,bool);
IMPORT	int	find_overtake_node_states(float*,RP_DATA*,bool*,bool);
IMPORT	int	find_steady_ahead_speed(ANGLE_DIRECTION,Locstate,Locstate,
					Locstate,float,float*,float*,
					float*,float*);
IMPORT	void	set_temp_mnode_normal(NORMAL_FUNCTION*);
#endif /* defined(TWOD) */

#if defined(FULL_PHYSICS) && defined(TWOD)
/*	gnode/ganom.c */
IMPORT	int	anomalous_reflection_propagate(Front*,Wave*,NODE*,NODE*,
					       RPROBLEM**,float,float*,
					       NODE_FLAG,O_CURVE**,O_CURVE**,
					       bool*,RP_DATA*,POINT*,float*,
					       BOND**);
IMPORT	int	refl_curve_overtakes_incident_shock(O_CURVE**,O_CURVE**,BOND**,
						    POINT*,Front*,Wave*,
						    RPROBLEM**,float,float*,
						    int,ANGLE_DIRECTION,
						    NODE_FLAG);

/*	gnode/gscnode.c */
IMPORT	ANGLE_DIRECTION	incident_shock_orientation(int,ORIENTATION);
IMPORT	bool	curves_at_shock_diffraction(NODE*,O_CURVE**,ANGLE_DIRECTION*,
					    bool);
IMPORT	int	diffraction_node_propagate(Front*,Wave*,NODE*,NODE*,
					   RPROBLEM**,float,float*,NODE_FLAG,
					   POINTER);
IMPORT	int	modify_diffraction_node(POINT*,NODE*,NODE*,O_CURVE**,O_CURVE**,
					BOND**,Front*,Wave*,float,
					RP_DATA*,NODE_FLAG);
IMPORT	int	trans_node_parameter(void);
IMPORT	int	set_trans_node_parameter(int);
IMPORT	int	sonic_incident_shock_at_diffraction(NODE*,O_CURVE**,O_CURVE**,
						    BOND**,POINT*,float*,
						    RP_DATA*,Front*,Wave*,
						    RPROBLEM**,float,float*,
						    NODE_FLAG);
IMPORT	int	transmission_node_propagate(Front*,Wave*,NODE*,NODE*,
					    RPROBLEM**,float,float*,NODE_FLAG);
IMPORT	void	set_temp_tnode_normal(NORMAL_FUNCTION*);

/*	gnode/gscnsts.c */
IMPORT	int	find_transmission_node_states(float*,float**,RP_DATA*,int,
					      int,NODE_FLAG);
IMPORT	int	is_regular_diffraction_node(float*,float*,float*,float**,
					    RP_DATA*,RP_DATA*,bool*,Front*,
					    int,NODE_FLAG);
IMPORT	int	is_small_inc_ang_reg_diff_node(float**,RP_DATA*,bool*,NODE_FLAG);

#endif /* defined(FULL_PHYSICS) && defined(TWOD) */

/*	gnode/gcurve.c */
#if defined(THREED)
IMPORT	void	g_curve_propagate_3d(Front*,POINTER,CURVE*,CURVE*,float);
#endif /* defined(THREED) */

#endif /* !defined(_GNODEPROTOS_H) */
