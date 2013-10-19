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
*				iprotos.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#if !defined(_IPROTOS_H)
#define _IPROTOS_H

#include <intfc/int.h>

 /* libintfc.a EXPORTED Function Declarations*/

/*	comp.c*/
IMPORT	COMPONENT	component(float*,INTERFACE*);
IMPORT	COMPONENT	long_component(float*,INTERFACE*);
IMPORT	COMPONENT	nearest_interior_comp(bool,COMPONENT,
					      float*,INTERFACE*);
IMPORT	COMPONENT	new_component(COMPONENT);
IMPORT	bool	is_excluded_comp(COMPONENT,INTERFACE*);
IMPORT	int	check_comps_at_nodes(INTERFACE*,O_NODE**);
IMPORT	int	comps_consistent_at_node(NODE*);
IMPORT	void	exclude_comp(COMPONENT,INTERFACE*);
IMPORT	void	set_topological_grid(INTERFACE*,RECT_GRID*);
IMPORT	void	show_COMP(FILE*,INTERFACE*);

/*	iblkb.c*/
IMPORT	int construct_bdry_blk(BLK_CRX*,BLK_TRI*);

/*	iecomps.c*/
IMPORT	bool	equivalent_comps(COMPONENT,COMPONENT,INTERFACE*);
IMPORT	const COMPONENT	*equivalent_components_list(COMPONENT,int*,INTERFACE*);
IMPORT  void	fprint_comp_equiv_lists(FILE*,INTERFACE*);
IMPORT  void	fprint_comp_equiv_lists(FILE*,INTERFACE*);
IMPORT	void	set_equivalent_comps(COMPONENT,COMPONENT,INTERFACE*);

/*	ifourier.c*/
IMPORT  FOURIER_POLY    *get_fourier_bubble(float*,float*,int,const char*);
IMPORT  FOURIER_POLY    *get_fourier_multi_mode(float*,float*,int,const char*);
IMPORT  FOURIER_POLY    *get_fourier_mixed(float*,float*,int,const char*);
IMPORT  FOURIER_POLY    *get_fourier_random(float*,float*,int,const char*);
IMPORT	FOURIER_POLY	*allocate_fourier_poly(int,int,ALIGN*);
IMPORT	FOURIER_POLY	*get_fourier_coeffs(float*,float*,int,const char*);
IMPORT	LEGENDRE_POLY	*get_legendre_coeffs(float,const char*);
IMPORT  LEGENDRE_POLY   *get_legendre_mixed(float,const char*);
IMPORT  LEGENDRE_POLY   *get_legendre_multi_mode(float,const char*,const char*);
IMPORT  LEGENDRE_POLY   *get_legendre_random(float,const char*);
IMPORT  LEGENDRE_POLY   *allocate_legendre_poly(int,ALIGN*);
IMPORT	float		fourier_poly(float*,FOURIER_POLY*);
IMPORT	float		legendre_poly(float,LEGENDRE_POLY*);
IMPORT	int		random_bubble_num_modes(const char*,int*,int*,int);
IMPORT	int		spherical_num_modes(const char*,int*,int*,float*);
IMPORT	void		init_random_modes(int,int,int,int,FOURIER_POLY*,
					  float*,float*);

/*	intfc.c*/
IMPORT	BOND	*i_Bond(POINT*,POINT*);
IMPORT	CURVE	**split_curve(POINT*,BOND*,CURVE*,
			      COMPONENT,COMPONENT,COMPONENT,COMPONENT);
IMPORT	CURVE	*curve_of_bond(BOND*,INTERFACE*);
IMPORT	CURVE	*i_copy_curve(CURVE*,NODE*,NODE*);
IMPORT	CURVE	*i_make_curve(COMPONENT,COMPONENT,NODE*,NODE*);
IMPORT	CURVE	*join_curves(CURVE*,CURVE*,COMPONENT,COMPONENT,BOND**);
IMPORT	INTERFACE	*i_copy_interface(INTERFACE*);
IMPORT	INTERFACE	*current_interface(void);
IMPORT	INTERFACE	*make_interface(int);
IMPORT	INTERFACE	*read_print_interface(INIT_DATA*,const IO_TYPE*,
                                              bool,int*);
IMPORT	INTERFACE	*read_interface(void);
IMPORT	HYPER_SURF *i_make_hypersurface(COMPONENT,COMPONENT);
IMPORT	HYPER_SURF_BDRY	*i_make_hypersurface_boundary(void);
IMPORT	NODE	*i_copy_node(NODE*);
IMPORT	NODE	*i_make_node(POINT*);
IMPORT	NODE	*node_of_point(POINT*,INTERFACE*);
IMPORT	POINT	*i_Point(float*);
IMPORT	POINT	*i_Static_point(INTERFACE*);
IMPORT	POINT	*i_copy_point(POINT*);
IMPORT	POINT	*i_make_point(float*,COMPONENT,COMPONENT);
IMPORT	POINTER	init_table_Store(size_t,INIT_DATA*);
IMPORT	POINTER	Store(size_t);
IMPORT	POINTER	store(size_t);
IMPORT	bool	exists_interface(INTERFACE*);
IMPORT	bool	i_delete_end_of_bond(BOND*,CURVE*);
IMPORT	bool	i_delete_node(NODE*);
IMPORT	bool	i_delete_start_of_bond(BOND*,CURVE*);
IMPORT	bool	i_insert_point_in_bond(POINT*,BOND*,CURVE*);
IMPORT	bool	next_bond(INTERFACE*,BOND**,CURVE**);
IMPORT	bool	next_curve(INTERFACE*,CURVE**);
IMPORT	bool	next_hypersurface1d(INTERFACE*,HYPER_SURF**);
IMPORT	const char *i_boundary_type_as_string(int);
IMPORT	float	i_random01(INTERFACE*);
IMPORT	int	i_delete_curve(CURVE*);
IMPORT	int	i_delete_interface(INTERFACE*);
IMPORT	int	i_delete_point(POINT*);
IMPORT	struct	Table	*interface_table_list(void);
IMPORT	struct	Table	*table_of_interface(INTERFACE*);
IMPORT	uint64_t	bond_number(BOND*,INTERFACE*);
IMPORT	uint64_t	curve_number(CURVE*);
IMPORT	uint64_t	hypersurface_boundary_number(HYPER_SURF_BDRY*);
IMPORT	uint64_t	hypersurface_element_number(HYPER_SURF_ELEMENT*,
						    INTERFACE*);
IMPORT	uint64_t	hypersurface_number(HYPER_SURF*);
IMPORT	uint64_t	interface_number(INTERFACE*);
IMPORT	uint64_t	node_number(NODE*);
IMPORT	uint64_t	point_number(POINT*);
IMPORT	uint64_t	table_number(struct Table*);
IMPORT	int	i_read_boundary_type_from_string(const char*);
IMPORT	void	delete_from_cross_list(CROSS*);
IMPORT	void	fprint_hypersurface(FILE*,HYPER_SURF*);
IMPORT	void	i_fprint_interface(FILE*,INTERFACE*);
IMPORT	void	print_hypersurface(HYPER_SURF*);
IMPORT	void	print_hypersurface_boundaries(HYPER_SURF_BDRY**);
IMPORT	void	set_current_interface(INTERFACE*);
IMPORT	POINT *i_average_points(bool,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,
	                             POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*);
IMPORT	ORIENTATION orientation_of_bond_at_tri(BOND*,TRI*);
IMPORT	bool	delete_side_of_tri(TRI*,SURFACE*,int);
IMPORT	bool	delete_vertex_of_tri(POINT*pt,TRI*,SURFACE*);
EXPORT	bool	flip_diagonal(TRI*,int);
EXPORT	bool	retriangulate_polygon(POINT**,int,POINT**,int,const float*,
				      TRI**,int,SURFACE*,BOND_TRI**,TRI**,
				      TRI***,int*);
IMPORT	int side_of_tri_with_bond(BOND*,TRI*);
IMPORT	uint64_t	bond_tri_number(BOND_TRI*,INTERFACE*);
IMPORT	uint64_t	surface_number(SURFACE*);
IMPORT	uint64_t	tri_number(TRI*,INTERFACE*);
IMPORT  void 	order_interface(INTERFACE *); 
IMPORT  void 	delete_scn(SURFACE *); 
IMPORT  bool    change_node_of_closed_curve(POINT*, CURVE*);
IMPORT  void    reset_nodes_posn(INTERFACE *);


/*	irefl.c*/
IMPORT	void	i_reflect_node(NODE*,float*,float*);
IMPORT	void	i_reflect_curve(CURVE*,float*,float*);
IMPORT	void	i_reflect_point(POINT*,float*,float*,INTERFACE*);
IMPORT	void	i_reflect_surface(SURFACE*,float*,float*);

/*	iscatter.c*/
IMPORT	int	domain_id(int*,int*,int);
IMPORT	int	neighbor_id(int*,int*,int,int,PP_GRID*);
IMPORT	void	find_Cartesian_coordinates(int,PP_GRID*,int*);
IMPORT	void	print_PP_GRID_structure(PP_GRID*);
IMPORT  PP_GRID *set_pp_grid(INIT_DATA*,RECT_GRID*);

/*	isub.c*/
IMPORT  const char *grid_direction_name(GRID_DIRECTION);
IMPORT  void    init_seg_crx_lists(INTERFACE*,int,int);
IMPORT	BDRY_SIDE rect_bdry_side_for_curve(int*,int*,CURVE*,RECT_GRID*);
IMPORT	P_LINK	*add_to_hash_table(POINTER,POINTER,P_LINK*,int);
IMPORT	POINTER	find_from_hash_table(POINTER,P_LINK*,int);
IMPORT	bool	i_is_comp_in_list(COMPONENT,COMP_LIST*);
IMPORT	float	scalar_product_on_bonds(BOND*,BOND*,int);
IMPORT	float	scaled_bond_length(BOND*,float*,int);
IMPORT	void	scaled_tri_params(TRI*,float*,float*,float*);
#if !defined(separation)
IMPORT float separation(POINT*,POINT*,int);
#endif /* !defined(separation) */
IMPORT	int	curve_in_curve_list(CURVE*,CURVE**);
IMPORT	bool	is_c_on_intfc(CURVE*);
IMPORT	bool	is_b_on_curve(CURVE*,BOND*);
IMPORT	bool    pointer_in_list(POINTER,int,POINTER*);
IMPORT	bool    integer_in_list(int,int,int*);
IMPORT	void	i_add_comp_to_list(COMPONENT,COMP_LIST*,INTERFACE*);
IMPORT	void	rect_bdry_side_for_hyper_surf(int*,int*,HYPER_SURF*,
					      RECT_GRID*);
IMPORT	void	reset_hash_table(P_LINK*,int);
IMPORT	void	reset_intfc_num_points(INTERFACE*);
IMPORT	void	vector_product_on_bonds(BOND*,BOND*,int,float*);
IMPORT	float	area_of_closed_curve(CURVE*);
IMPORT	ANGLE_DIRECTION	c1_to_c2_direction(O_CURVE*,O_CURVE*);
IMPORT	bool	robust_cross_bonds(BOND*,int,BOND*,int,float*,float*,
				   RECT_GRID*,POINT*);
IMPORT	int	cross_sign(BOND*,BOND*);
IMPORT	int	intersect_bond_with_curve_segment(BOND*,BOND*,BOND*,O_CURVE*,
					  BOND**,POINT*,RECT_GRID*);
IMPORT	int	is_short_curve(CURVE*,ORIENTATION,RECT_GRID*,float);
IMPORT	int	robust_cross_bond_circle(BOND*,POINT*,float,float*,POINT*);
IMPORT	int	robust_extend_bond_to_cross_bond(BOND*,ORIENTATION,BOND*,
					float*,float*,POINT*,float*,int);
IMPORT	int	robust_extend_bonds_to_cross(BOND*,ORIENTATION,
					     int,BOND*,ORIENTATION,int,
					     POINT*,float*,float*,POINT*,
					     RECT_GRID*);
IMPORT	int	robust_quad_roots_in_interval(double*,double,double,double,
					      double,double,double);
IMPORT	void	big_angle(BOND*,CURVE*,BOND*,CURVE*,float*,float*,RECT_GRID*);
IMPORT 	int  	seg_index2d(int,int,GRID_DIRECTION,int*);
IMPORT	int	set_tri_list_around_point(POINT*,TRI*,TRI***,INTERFACE*);
IMPORT	int	tri_list_computed_by_normal(POINT*,TRI*,TRI***,INTERFACE*);
IMPORT	const float* const *side_vector(const TRI*);
IMPORT	const float        *vector_on_tri_side(const TRI*,int,float*);
IMPORT	const float        *length_side(const TRI*);
IMPORT	const float	   *Tri_normal(const TRI*);
IMPORT	float              length_of_tri_side(const TRI*,int);
IMPORT	float              sqr_norm(const TRI*);
IMPORT	void	area_weighted_normal3d(POINT*,HYPER_SURF_ELEMENT*,
                                       HYPER_SURF*,float*);
IMPORT	void	sine_weighted_normal3d(POINT*,HYPER_SURF_ELEMENT*,
                                       HYPER_SURF*,float*);
IMPORT	void	average_position_of_surface(float*,SURFACE*);
IMPORT	void	omit_vertex_in_plane_fit(void);
IMPORT	void	plane_fit_normal3d(POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,
				   float*);
IMPORT	void	reset_normal_on_intfc(INTERFACE*);
IMPORT	void	set_normal_of_tri(TRI*);
IMPORT	void	PointArrayRing1(POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,int*,
					POINT**);
IMPORT	void	PointArrayRing2(POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,int*,
					int*, POINT**,POINT**);
IMPORT  int     two_points_share_side(POINT*,TRI*,POINT*,INTERFACE*);
IMPORT 	int  	seg_index3d(int,int,int,GRID_DIRECTION,int*);
IMPORT 	int  	face_index3d(int,int,int,GRID_DIRECTION,int*);
IMPORT  void    init_face_crx_lists(INTERFACE*,int,int);
IMPORT  void    set_normal_from_tris(POINT*,TRI**,int,float*);
IMPORT  void    plane_fit_normal3d_along_wall(float*, POINT*, TRI**,int,TRI**,int);

IMPORT  void print_crxings(CRXING*,bool);
IMPORT  const char *crossing_direction_name(CROSSING_DIRECTION);
IMPORT bool WLSP_compute_normal2d(POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*);
IMPORT bool WLSP_compute_normal3d(POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*);
IMPORT bool WLSP_compute_normal3d0(POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*);

/*	iuserintfc.c*/
IMPORT	I_USER_INTERFACE	*i_user_hook(int);
IMPORT	void	SetChunkSize(size_t);
IMPORT	void	i_preserve_user_hooks(int,PRESERVE_USER_HOOKS);

/*	ppcopy.c*/
IMPORT	INTERFACE	*i_receive_interface(int);
IMPORT	POINTER	_new_address(INTERFACE*,POINTER,POINTER*,POINTER*,int);
IMPORT	void	i_reconstruct_bond_pointers(BOND*,INTERFACE*,INTERFACE*,
					    POINTER*,POINTER*,int);
IMPORT	void	i_reconstruct_curve_pointers(CURVE*,INTERFACE*,INTERFACE*,
					     POINTER*,POINTER*,int);
IMPORT	void	i_reconstruct_interface_pointers(INTERFACE*,struct Table*,
						 POINTER*,POINTER*);
IMPORT	void	i_reconstruct_node_pointers(NODE*,INTERFACE*,INTERFACE*,
					    POINTER*,POINTER*,int);
IMPORT	void	i_reconstruct_point_pointers(POINT*,INTERFACE*,INTERFACE*,
					     POINTER*,POINTER*,int);
IMPORT	void	i_send_interface(INTERFACE*,int);

/*	shift.c*/
IMPORT	INTERFACE *remap_interface(INTERFACE*,
				   void (*)(POINT*,BOND*,CURVE*,POINT*,BOND*,
					   CURVE*,bool,RECT_GRID*,POINTER),
				   void (*)(INTERFACE*,INTERFACE*,
					   void (*remap)(POINT*,BOND*,CURVE*,
							POINT*,BOND*,CURVE*,
							bool,RECT_GRID*,
							POINTER),
					   POINTER),
				   POINTER);

/*	top.c*/
IMPORT	BDRY_SIDE boundary_side(const float*,const RECT_GRID*,float);
IMPORT	BDRY_SIDE nearest_boundary(float*,RECT_GRID*);
IMPORT	BOND	*random_bond_on_curve(CURVE*);
IMPORT	CURVE	*adjacent_curve(CURVE*,ORIENTATION,ANGLE_DIRECTION,
				ORIENTATION*);
IMPORT	CURVE	*i_attach_curve_to_node(CURVE*,POINT*,BOND*,NODE*);
IMPORT	CURVE	*i_make_fourier_curve(int,float,float,FOURIER_POLY*,
				      COMPONENT,COMPONENT);
IMPORT	CURVE	*read_curve_from_file(COMPONENT,COMPONENT,NODE*,NODE*,char*);
IMPORT	O_NODE	*make_onode(NODE*);
IMPORT	bool	bonds_in_strict_order(BOND*,BOND*);
IMPORT	bool	points_in_strict_order(POINT*,BOND*,POINT*,BOND*,int);
IMPORT	const char *bdry_side_as_string(BDRY_SIDE);
IMPORT	float	curve_length(CURVE*);
IMPORT	int	intersect_ray_with_boundary(float*,float*,float*,float*,
					    float*,int);
IMPORT	int	intersect_ray_with_curve(POINT*,float*,BOND*,BOND*,
					 CURVE*,ORIENTATION,BOND**,POINT*);
IMPORT	int	intersect_ray_with_sector(POINT*,POINT*,float*,float**,
					  float*,int);
IMPORT	bool	i_move_closed_loop_node(CURVE*,BOND*);
IMPORT	int	num_curves_at_node(NODE*,int*,int*);
IMPORT	void	change_node_of_curve(CURVE*,ORIENTATION,NODE*);
IMPORT	void	i_cut_curve(POINT*,BOND*,CURVE*,ORIENTATION);
IMPORT	void	copy_o_curve(O_CURVE*,O_CURVE*);
IMPORT	void	delete_list(POINTER**);
IMPORT	void	i_invert_curve(CURVE*);
IMPORT	void	i_reverse_curve(CURVE*);
IMPORT	void	invert_bond(BOND*);
IMPORT	void	merge_and_delete_nodes(NODE*,NODE*);
IMPORT	void	nearest_boundary_point(float*,float*,RECT_GRID*);
IMPORT	void	plot_interface(INTERFACE*,const char*,int*,int*,const char*);
IMPORT	void	print_bdry_side(const char*,BDRY_SIDE,const char*);
IMPORT	void	print_bond(BOND*);
IMPORT	void	print_bond_list(CURVE*);
IMPORT	void	print_curve_with_orient(CURVE*,ORIENTATION);
IMPORT	void	print_int_vector(const char*,const int*,int,const char*);
IMPORT	void	print_o_curve(O_CURVE*);
IMPORT	void	print_o_curve_family(O_CURVE_FAMILY*);
IMPORT	void	print_onode(O_NODE*);
IMPORT	void	print_onode_list(O_NODE**);
IMPORT	void	rbl_after_move_node(NODE*);
IMPORT	void	set_bond_length(BOND*,int);
IMPORT	void	set_point_of_bond(POINT*,BOND*,ORIENTATION,int);
IMPORT	void	update_num_points(INTERFACE*);
IMPORT	int	num_points_on_curve(CURVE*);

/*	zoom.c*/
IMPORT	INTERFACE	*i_zoom_interface(INTERFACE*,RECT_GRID*,
					  float*,float*,float**);
IMPORT	void    clip_interface2d(INTERFACE*);
IMPORT	void	rotate_and_zoom_rect_grid(RECT_GRID*,float*,float*,float**);

/*	comp1d.c*/
IMPORT	bool consistent_components1d(INTERFACE*);
IMPORT	bool	make_point_comp_lists(INTERFACE*);
IMPORT	void	show_point_comp_lists(INTERFACE*);
IMPORT  void    reset_intfc_components(INTERFACE*);

/*	setb1d.c*/


/*	comp2d.c*/
IMPORT	bool	make_bond_comp_lists(INTERFACE*);
IMPORT	float	shortest_distance2d(float*,BOND*,POINT**,float*,SIDE*);
IMPORT  void    show_box_comp_crx(int*,int*,int*,COMPONENT*,int*);

/*	cross2d.c*/
IMPORT	int	find_companion_cross(CROSS*,CROSS**,ORIENTATION*,ORIENTATION*,
				     ORIENTATION*,ORIENTATION*);
IMPORT	void	add_to_cross_list(CROSS**,CURVE*,BOND*,CURVE*,BOND*,POINT*);
IMPORT	void	insert_in_cross_list(CROSS*,CROSS*);
IMPORT	void	print_cross(CROSS*);
IMPORT	void	print_cross_list(CROSS*);
IMPORT	void	rcl_after_delete_bond_fragment_at_node(CROSS*,POINT*,
						       CURVE*,ORIENTATION);
IMPORT	void	rcl_after_insert_point(CROSS*,POINT*,BOND*);
IMPORT	void	rcl_after_join(CROSS*,CURVE*,CURVE*,CURVE*);
IMPORT	void	rcl_after_split(CROSS*,POINT*,BOND*,CURVE*,CURVE**);

/*	isect2d.c*/
IMPORT	bool	bond_crosses_curve(BOND*,CURVE*,POINT*,BOND**,ORIENTATION);
IMPORT	bool	cross_bonds(BOND*,BOND*,POINT*);

/*	isurgery.c */
IMPORT	void	identify_detached_surface_curve_pair(INTERFACE*);
IMPORT	int 	next_null_sided_tri(TRI*,POINT*,TRI**);
IMPORT	int 	prev_null_sided_tri(TRI*,POINT*,TRI**);

/*	setb1d.c*/
IMPORT  bool    i_set_boundary1d(INTERFACE*,RECT_GRID*,COMPONENT,float);
/*	setb2d.c*/
IMPORT	bool	i_set_boundary2d(INTERFACE*,RECT_GRID*,COMPONENT,float);
/*      igview.c */
IMPORT  void    geomview_intfc_plot2d(const char*,INTERFACE*,RECT_GRID*);
IMPORT  void    tecplot_interface(const char*,INTERFACE*);
IMPORT  void    tecplot_blk_intfc_plot(const char *,BLK_TRI *);
IMPORT  void    tecplot_curve(const char*,FILE*,CURVE*);
IMPORT  void    tecplot_surface(const char*,FILE*,SURFACE*);
IMPORT  void    tecplot_surface_in_ball(const char*,SURFACE*);
IMPORT  void    tecplot_box(const char*,FILE*,float*,float*);
IMPORT	void	tecplot_tris(const char *,TRI **, int);
IMPORT  void    tecplot_show_tris(const char*, TRI **, int, FILE *);
IMPORT	void	tecplot_show_box_tri(const char*,RECT_BOX*,TRI**,int,FILE*);
IMPORT	void	tecplot_show_box_tris(const char*,TRI**,int,RECT_GRID*,int*);
IMPORT	void    tecplot_interface_in_box(const char*,FILE*,int*,int*,
				INTERFACE*);
IMPORT	void	set_shift_for_tecplot(float, float, float);
IMPORT	void    tecplot_surface_in_ball(const char *, SURFACE*);
IMPORT	void	set_tst_posn(float*);
IMPORT	void	vtk_interface_plot(const char*,INTERFACE*,bool);
IMPORT	void	sdl_interface_plot(const char*,INTERFACE*);
IMPORT  void    gd_initplot(char*,char*,float,float,float,float,int);
IMPORT  void    gd_plotdata(int,float*,float*);
IMPORT  void    gd_plotframe(char*);
IMPORT  void    gd_closeplot();

/*	comp3d.c*/
IMPORT	bool	make_tri_comp_lists(INTERFACE*);
IMPORT	bool	within_tri(const float*,const float*,const float*,
	                   const float*,const float*,float);
IMPORT	bool	i_insert_point_in_tri(POINT*,TRI*,SURFACE*);
IMPORT	bool	i_insert_point_in_tri_side(POINT*,int,TRI*,SURFACE*);
IMPORT	bool	i_undo_insert_point_in_tri(POINT*,TRI*,SURFACE*);
IMPORT	bool	i_undo_insert_point_in_tri_side(POINT*,int,TRI*,SURFACE*);
IMPORT  bool    tri_edge_crossing(TRI*,float*,float*,int,int*,int*,float*);
IMPORT	void	assign_tri_icoords(RECT_GRID*,TRI*);
IMPORT  bool    line_tri_crossing(float*,TRI*,float*,float*,float);
IMPORT  bool    line_point_projection(float*,int*,float*,float*,float*,float);
IMPORT  bool    is_tri_outside(INTERFACE*,TRI*,RECT_GRID*);
IMPORT  bool    is_tri_outside_box(TRI*, float **);
IMPORT  bool    is_outside_surface(INTERFACE*,SURFACE*,RECT_GRID*);
IMPORT  bool    is_outside_surfaces(INTERFACE*,RECT_GRID*);
IMPORT  void    delete_outside_surface(INTERFACE *);

/* 	idiagnostic.c*/
IMPORT	int  	index_of_pointer(POINTER*,POINTER);
IMPORT  int 	points_on_surface(SURFACE*);
IMPORT 	void 	print_blk_tri(BLK_TRI*);
IMPORT	void 	print_bond(BOND*);
IMPORT  void 	points_of_interface(INTERFACE*);
IMPORT	bool	the_tri(TRI*);
IMPORT	bool	the_side(TRI*);
IMPORT	bool	the_BOND(BOND*);
IMPORT	bool	the_point(POINT*);
IMPORT	bool	the_pt(float *);
IMPORT 	bool	check_tri_and_neighbor(TRI*);
IMPORT	void	print_tri_coords(TRI*);
IMPORT  void 	find_blk_tri(BLK_TRI *);
	/*igview.c*/
IMPORT	void	gview_bounding_box(FILE*,const float*,const float*,
	                           int,const char*);
IMPORT  void    gview_cube(FILE*,const float*,const float*);
IMPORT  void 	gview_local_surface(SURFACE*,const char*,const char*,
				    SURFACE_COLOR,const float*,float);
IMPORT	void	gview_plot_axes(const char*,const char*,const float*,
				const float*,const float*,const float*);
IMPORT	void	gview_plot_coord_sys(const char*,const char*,const float*,
	                             const float*,const float*,const float*,
				     const float*,const float*);
IMPORT	void	gview_plot_polyline(const char*,const char*,POINT**,
	                            int,bool,float,float,float,float,
				    const float*,const float*);
IMPORT	void	gview_plot_triangle_list(const char*,const char*,TRI**,int,
	                                 float,float,float,float,float,float,
					 float,const float*,const float*);
IMPORT	void	gview_plot_vertices(const char*,const char*,POINT**,
	                            int,const float*,const float*);
IMPORT	void	gview_plot_c_curve(const C_CURVE*,int,const char*);
IMPORT  void    gview_plot_curve(const CURVE*,const char*,const char*,
				 SURFACE_COLOR,int);
IMPORT	void	gview_plot_tri_and_point_list(const char*,TRI**,
	                                      const float*,int,float* const*,
					      SURFACE_COLOR,float,int,int);
IMPORT	void	gview_plot_tri_list(const char*,TRI**,int);
IMPORT	void	gview_intfc_within_range(const char*,INTERFACE*,
					float*,float);
IMPORT  void    gview_polyline(const char*,const char*,float* const*,
			       int,SURFACE_COLOR,int);
IMPORT  void 	gview_surface(SURFACE*,const char*,SURFACE_COLOR);
IMPORT	void	geomview_interface_plot(const char*,INTERFACE*,RECT_GRID*);
IMPORT	void	set_point_list_bounding_box(POINT**,int,float*,
	                                    float*,bool,bool);
IMPORT	void	set_tri_list_bounding_box(TRI**,int,float*,float*,
	                                  bool,bool);
IMPORT	void	set_vector_bounding_box(const float*,const float*,float,
	                                float*,float*,bool,bool);
IMPORT	void	gview_point_tri_rings(const char*,POINT*);

/*	int3d.c*/
IMPORT	BOND_TRI *i_link_tri_to_bond(BOND_TRI*,TRI*,SURFACE*,BOND*,CURVE*);
IMPORT	SURFACE	*i_copy_surface(SURFACE*,CURVE**,CURVE**,bool);
IMPORT	SURFACE *i_join_surfaces(CURVE*);
IMPORT	SURFACE	*i_make_surface(COMPONENT,COMPONENT,CURVE**,CURVE**);
IMPORT	SURFACE	*i_read_surface(INTERFACE*,int);
IMPORT	TRI	*i_make_tri(POINT*,POINT*,POINT*,POINTER,POINTER,POINTER,int);
IMPORT	TRI	*Next_tri_at_vertex(TRI*,POINT*);
IMPORT	TRI	*Prev_tri_at_vertex(TRI*,POINT*);
IMPORT	bool    curve_is_in_surface_bdry(SURFACE*,CURVE*,ORIENTATION*);
IMPORT	bool    find_surfaces_to_join_at_curve(CURVE*,SURFACE**,SURFACE**);
IMPORT	bool    i_sort_bond_tris(INTERFACE*);
IMPORT	bool	next_tri(INTERFACE*,TRI**,SURFACE**);
IMPORT	bool	remove_curve_from_surface_bdry(SURFACE*,CURVE*,ORIENTATION);
IMPORT	bool	tris_on_side_of_bond_for_join(BOND*,TRI**,int*,TRI**,int*);
IMPORT	bool	link_neighbor_tris(TRI*,TRI*);
IMPORT	int	i_delete_surface(SURFACE*);
IMPORT	void 	i_reverse_bond(BOND*);
IMPORT	void 	i_reorder_curve_link_list(CURVE*);
IMPORT	void	insert_tri_at_head_of_list(TRI*,SURFACE*);
IMPORT	void	insert_tri_at_tail_of_list(TRI*,SURFACE*);
IMPORT	void	install_curve_in_surface_bdry(SURFACE*,CURVE*,ORIENTATION);
IMPORT	void	link_tri_list_to_surface(TRI*,TRI*,SURFACE*);
IMPORT	void	null_tri_array_numbers(INTERFACE*);
IMPORT	void	print_tri(TRI*,INTERFACE*);
IMPORT	void	remove_tri_from_surface(TRI*,SURFACE*,bool);
IMPORT	void	rotate_triangle(TRI*,int);
IMPORT  bool    assign_btri_states(BOND_TRI*, BOND_TRI*);
IMPORT  SURFACE *detach_one_surface(SURFACE *);
IMPORT  void    print_wall_crx(const char*,int*,int,int,CRXING*);
IMPORT  void    print_wall_curve_crx(const char*,int*,int,int,CRXING*);
IMPORT  void    print_wall_curve_crx0(const char*,POINT *, int,CRXING*);

/*	iprt3d.c*/
IMPORT	void	print_c_bond(C_BOND*,INTERFACE*);
IMPORT	void	print_c_curve(C_CURVE*,INTERFACE*);
IMPORT	void	print_c_surf(C_SURF*,INTERFACE*);
IMPORT	void	print_c_surf_flag(C_SURF_FLAG*);


/*	isect3d.c*/
IMPORT	C_BOND	*i_CBond(C_BOND*,POINT*,POINT*,TRI*,TRI*);
IMPORT	void	i_print_intersections3d(CROSS*,INTERFACE*);

/*	map.c*/
IMPORT  int     NumOfCurvePoints(CURVE*);
IMPORT  int     NumOfSurfPoints(SURFACE*);
IMPORT  int     NumOfIntfcPoints(INTERFACE*);
IMPORT  int     NumOfNodes(INTERFACE*);
IMPORT  int     NumOfCurveBonds(CURVE*);
IMPORT  int     NumOfIntfcBonds(INTERFACE*);
IMPORT  int     NumOfCurves(INTERFACE*);
IMPORT  int     NumOfSurfaces(INTERFACE*);
IMPORT  int     NumOfSurfTris(SURFACE*);
IMPORT  int     NumOfIntfcTris(INTERFACE*);
IMPORT  int     NumOfInteriorPoints(INTERFACE*);

IMPORT  void    ArrayOfCurves(INTERFACE*intfc,CURVE**);
IMPORT  void    ArrayOfSurfaces(INTERFACE*,SURFACE**);
IMPORT  void    ArrayOfCurvePoints(CURVE*,float*);
IMPORT  void    ArrayOfIntfcPoints(INTERFACE*intfc,float*);
IMPORT  void    ArrayOfSurfTris_FT(SURFACE*,TRI**);
IMPORT  void    ArrayOfSurfTris(SURFACE*,float*,int*);
IMPORT  void    ArrayOfIntfcTris_FT(INTERFACE*,TRI**);
IMPORT  void    ArrayOfIntfcTris(INTERFACE*,float*,int*);
IMPORT	int	GridSegCrossing(CRXING**,int*,GRID_DIRECTION,INTERFACE*);
IMPORT	COMPONENT *GridIntfcComp(INTERFACE*);
IMPORT	bool	IntfcGetPointChain(POINT*,POINT**,int);

IMPORT	int	GridSegCrossing(CRXING**,int*,GRID_DIRECTION,INTERFACE*);
IMPORT	COMPONENT *GridIntfcComp(INTERFACE*);

/*	setb3d.c*/
IMPORT	bool	i_set_boundary3d(INTERFACE*,RECT_GRID*,COMPONENT,float);

/*	trisurf.c*/
IMPORT	void	oblique_planar_surface_triangulation(SURFACE*,RECT_GRID*); 
IMPORT	void	planar_hole_surface_triangulation(SURFACE*,RECT_GRID*,
						  POINT*,POINT*,POINT*,POINT*);
IMPORT	void	planar_surface_triangulation(SURFACE*,RECT_GRID*,const bool);

/*	icheck3d.c*/
IMPORT	void	null_sides_are_consistent(void);
IMPORT	bool	i_consistent_interface(INTERFACE*);
IMPORT	void	check_double_cone_point(INTERFACE*);

/*	userhooks.c */
IMPORT	BOND	*Bond(POINT*,POINT*);
IMPORT	CURVE	*attach_curve_to_node(CURVE*,POINT*,BOND*,NODE*);
IMPORT	CURVE	*copy_curve(CURVE*,NODE*,NODE*);
IMPORT	CURVE	*make_curve(COMPONENT,COMPONENT,NODE*,NODE*);
IMPORT	CURVE	*read_curve(INTERFACE*,int);
IMPORT	HYPER_SURF *make_hypersurface(COMPONENT,COMPONENT);
IMPORT	HYPER_SURF_BDRY	*make_hypersurface_boundary(void);
IMPORT	INTERFACE	*copy_interface(INTERFACE*);
IMPORT	INTERFACE	*zoom_interface(INTERFACE*,RECT_GRID*,
					float*,float*,float**);
IMPORT	NODE	*copy_node(NODE*);
IMPORT	NODE	*make_node(POINT*);
IMPORT	NODE	*read_node(INTERFACE*,int);
IMPORT	POINT	*Point(float*);
IMPORT	POINT	*Static_point(INTERFACE*);
IMPORT	POINT	*average_points(bool,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,
	                             POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*);
IMPORT	POINT	*copy_point(POINT*);
IMPORT	bool	delete_end_of_bond(BOND*,CURVE*);
IMPORT	bool	delete_node(NODE*);
IMPORT	bool	delete_start_of_bond(BOND*,CURVE*);
IMPORT	bool	insert_point_in_bond(POINT*,BOND*,CURVE*);
IMPORT	bool	intersections(INTERFACE*,CROSS**,const bool);
IMPORT	bool	is_subdomain_boundary(HYPER_SURF*);
IMPORT	bool	is_subdomain_node(NODE*);
IMPORT	bool	is_virtual_fixed_node(NODE*);
IMPORT	bool	long_nearest_interface_point(float*,COMPONENT,INTERFACE*,
					     USE_BOUNDARIES,HYPER_SURF*,float*,
					     float*,HYPER_SURF_ELEMENT**,
					     HYPER_SURF**);
IMPORT	bool	long_nearest_similar_interface_point(float*,COMPONENT,COMPONENT,
						     INTERFACE*,USE_BOUNDARIES,
						     HYPER_SURF*,float*,float*,
						     HYPER_SURF_ELEMENT**,
						     HYPER_SURF**);
IMPORT	bool	move_closed_loop_node(CURVE*,BOND*);
IMPORT	bool	nearest_interface_point(float*,COMPONENT,INTERFACE*,
					USE_BOUNDARIES,HYPER_SURF*,float*,
					float*,HYPER_SURF_ELEMENT**,
					HYPER_SURF**);
IMPORT	bool	nearest_similar_interface_point(float*,COMPONENT,COMPONENT,
						INTERFACE*,USE_BOUNDARIES,
						HYPER_SURF*,float*,float*,
						HYPER_SURF_ELEMENT**,
						HYPER_SURF**);
IMPORT	bool	make_interface_topology_lists(INTERFACE*);
IMPORT	bool	next_point(INTERFACE*,POINT**,HYPER_SURF_ELEMENT**,
			   HYPER_SURF**);
IMPORT	bool	next_hypersurface(INTERFACE*,HYPER_SURF**);
IMPORT	bool	set_boundary(INTERFACE*,RECT_GRID*,COMPONENT,float);
IMPORT	bool	user_join_curves(CURVE*,CURVE*,CURVE*);
IMPORT	bool	user_split_curve(int,POINT*,BOND*,CURVE*,CURVE**);
IMPORT	float	random01(INTERFACE*);
IMPORT	float	cross_tolerance(INTERFACE*);
IMPORT	int	delete_curve(CURVE*);
IMPORT	int	delete_interface(INTERFACE*);
IMPORT	int	print_number_of_tangles(const char*,INTERFACE*,CROSS*);
IMPORT	int	read_boundary_type_from_string(const char*,INTERFACE*);
IMPORT	int	read_print_intfc_rect_grids(const IO_TYPE*,INTERFACE*,REMAP*);
IMPORT	int	user_read_node(NODE*);
IMPORT	int	user_read_print_interface(INIT_DATA*,const IO_TYPE*,
                                          INTERFACE*,bool);
IMPORT	void	fprint_boundary_type(FILE*,const char*,int,const char*,
				     INTERFACE*);
IMPORT	void	fprint_curve(FILE*,CURVE*);
IMPORT	void	fprint_interface(FILE*,INTERFACE*);
IMPORT	void	fprint_intfc_rect_grids(FILE*,INTERFACE*);
IMPORT	void	fprint_node(FILE*,NODE*);
IMPORT	void	invert_curve(CURVE*);
IMPORT	void	reverse_curve(CURVE*);
IMPORT	void	print_curve(CURVE*);
IMPORT	void	print_interface(INTERFACE*);
IMPORT	void	print_node(NODE*);
IMPORT	void	print_intersections(CROSS*,INTERFACE*);
IMPORT	void	print_crossing_elements(CROSS*,INTERFACE*);
IMPORT	void	reflect_interface(INTERFACE*,float*,float*);
IMPORT	void	reflect_node(NODE*,float*,float*);
IMPORT	void	reflect_curve(CURVE*,float*,float*);
IMPORT	void	reflect_point(POINT*,float*,float*,INTERFACE*);
IMPORT	void	fset_hyper_surf_color(FILE*,HYPER_SURF*);
IMPORT	void	user_copy_hyper_surf(HYPER_SURF*,HYPER_SURF*);
IMPORT	void	user_fprint_curve(FILE*,CURVE*);
IMPORT	void	user_fprint_interface(FILE*,INTERFACE*);
IMPORT	void	user_fprint_intfc_rect_grids(FILE*,INTERFACE*);
IMPORT	void	user_fprint_node(FILE*,NODE*);
IMPORT	void	user_make_interface(INTERFACE*);
IMPORT	void	user_read_curve(CURVE*);
IMPORT	bool	user_read_print_curve(CURVE*,const IO_TYPE*,bool);
IMPORT	void	user_read_print_intfc_rect_grids(const IO_TYPE*,INTERFACE*,
						 bool,REMAP*);
IMPORT	void	user_read_print_node(NODE*,const IO_TYPE*,bool);
IMPORT	POINT	*make_point(float*,COMPONENT,COMPONENT);
IMPORT	POINT	*read_point(INTERFACE*,int);
IMPORT	POINT	*read_print_point(INTERFACE*,const IO_TYPE*,bool);
IMPORT	int	delete_point(POINT*);
IMPORT	void	fprint_point(FILE*,POINT*);
IMPORT	void	print_point(POINT*);
IMPORT	void	user_fprint_point(FILE*,POINT*);
IMPORT	void	user_read_print_point(POINT*,const IO_TYPE*,bool);
IMPORT	void	user_read_point(INTERFACE*,POINT*);
IMPORT	INTERFACE	*receive_interface(int);
IMPORT	void	reconstruct_bond_pointers(BOND*,INTERFACE*,INTERFACE*,
					  POINTER*,POINTER*,int);
IMPORT	void	reconstruct_curve_pointers(CURVE*,INTERFACE*,INTERFACE*,
					   POINTER*,POINTER*,int);
IMPORT	void	reconstruct_interface_pointers(INTERFACE*,struct Table*,
					       POINTER*,POINTER*);
IMPORT	void	reconstruct_node_pointers(NODE*,INTERFACE*,INTERFACE*,
					  POINTER*,POINTER*,int);
IMPORT	void	reconstruct_point_pointers(POINT*,INTERFACE*,INTERFACE*,
					   POINTER*,POINTER*,int);
IMPORT	void	send_interface(INTERFACE *,int);
IMPORT	BOND_TRI *link_tri_to_bond(BOND_TRI*,TRI*,SURFACE*,BOND*,CURVE*);
IMPORT	C_BOND	*CBond(C_BOND*,POINT*,POINT*,TRI*,TRI*);
IMPORT	SURFACE *copy_surface(SURFACE*,CURVE**,CURVE**,bool);
IMPORT	SURFACE *join_surfaces(CURVE*);
IMPORT	SURFACE *make_surface(COMPONENT,COMPONENT,CURVE**,CURVE**);
IMPORT	SURFACE *read_surface(INTERFACE*,int);
IMPORT	TRI	*make_tri(POINT*,POINT*,POINT*,POINTER,POINTER,POINTER,int);
IMPORT	bool	consistent_interface(INTERFACE*);
IMPORT	bool	insert_point_in_tri(POINT*,TRI*,SURFACE*);
IMPORT	bool	insert_point_in_tri_side(POINT*,int,TRI*,SURFACE*);
IMPORT	bool sort_bond_tris(INTERFACE*);
IMPORT	bool	undo_insert_point_in_tri(POINT*,TRI*,SURFACE*);
IMPORT	bool	undo_insert_point_in_tri_side(POINT*,int,TRI*,SURFACE*);
IMPORT	int	delete_surface(SURFACE*);
IMPORT	void 	reverse_bond(BOND*);
IMPORT	void 	reorder_curve_link_list(CURVE*);
IMPORT	void	assign_curve_boundary_flag(CURVE*);
IMPORT	void	gview_plot_interface(const char*,INTERFACE*);
IMPORT	void	fprint_surface(FILE*,SURFACE*);
IMPORT	void	print_surface(SURFACE*);
IMPORT	void	reconstruct_surface_pointers(SURFACE*,INTERFACE*,INTERFACE*,
					     POINTER*,POINTER*,int);
IMPORT	void	reconstruct_tri_pointers(TRI*,INTERFACE*,INTERFACE*,
					 POINTER*,POINTER*,int);
IMPORT	void	reflect_surface(SURFACE*,float*,float*);
IMPORT	void	user_fprint_surface(FILE*,SURFACE*);
IMPORT	void	user_install_faces(SURFACE*,int);
IMPORT	void	user_read_print_surface(SURFACE*,const IO_TYPE*,bool);
IMPORT	void	user_read_surface(SURFACE*);
IMPORT  void    check_print_intfc(const char*,const char*,char,INTERFACE*,
				int,int,bool);

/*      iredist.c */
IMPORT	void  	equi_redist_curve_seg(CURVE*,BOND*,BOND*,int,float,float,
					RECT_GRID*);
IMPORT	void	rect_bdry_curve_redist(CURVE*,ORIENTATION,RECT_GRID*,float*);
IMPORT  void    rect_bdry_redist2d(INTERFACE*,RECT_GRID*,int);
IMPORT 	bool 	i_delete_point_adjacent_to_node(CURVE*,ORIENTATION);
IMPORT  bool  	closed_curve_node_redistribute(INTERFACE*,bool);


/*      imkcurve.c */
IMPORT	CURVE	*make_elliptic_curve(ELLIP_PARAMS*,COMPONENT,COMPONENT,float);
IMPORT	void	coords_on_ellips(float*,float*,ELLIP_PARAMS*);
IMPORT  float 	multi_sine_mode_func(POINTER,float*);
IMPORT	void	prompt_make_level_curves(INTERFACE*,RECT_GRID*);
IMPORT  CURVE 	*prompt_make_linear_curve(INTERFACE*,RECT_GRID*);
IMPORT  CURVE 	*prompt_make_elliptic_curve(INTERFACE*,RECT_GRID*);
IMPORT  CURVE	*make_level_curve(RECT_GRID*,INTERFACE*,COMPONENT,COMPONENT,
			float (*func)(POINTER,float*),POINTER,bool);
IMPORT	CURVE	*make_array_curve(INTERFACE*,COMPONENT,COMPONENT,int,
			float**,bool);
IMPORT	POINTER init_ellipse_params(RECT_GRID*);
IMPORT	float 	ellipse_func(POINTER,float*);


/*      imksurf.c */
IMPORT  bool 	make_bdry_surfaces(INTERFACE*,RECT_GRID*);
IMPORT	bool    make_level_surface(RECT_GRID*,INTERFACE*,COMPONENT,COMPONENT,
                        float (*func)(POINTER,float*),POINTER,
                        SURFACE**);
IMPORT  bool    make_comp3_surfaces(RECT_GRID*,COMPONENT,COMPONENT,COMPONENT,
			float (*func1)(POINTER,float*),POINTER,
                        float (*func2)(POINTER,float*),POINTER,
                        SURFACE***,CURVE**);
IMPORT  bool    grid_line_crx_in_dir(float (*func)(POINTER,float*),
                        POINTER,int,float*,float*,float*,int);
IMPORT  bool    read_sdl_surface(INTERFACE*,COMPONENT,COMPONENT,char*,
			SURFACE**);
IMPORT  bool    read_vtk_surface(INTERFACE*,COMPONENT,COMPONENT,char*,
			SURFACE**);
IMPORT	float	dumbbell_func(POINTER,float*);
IMPORT	float	ellipsoid_func(POINTER,float*);
IMPORT	float	plane_func(POINTER,float*);
IMPORT	float	bdry_box_func(POINTER,float*);
IMPORT  SURFACE *prompt_make_level_surface(INTERFACE*,RECT_GRID*);
IMPORT  void 	prompt_make_bdry_surfaces(INTERFACE*,RECT_GRID*);
IMPORT  void    print_blk_crx(const BLK_CRX*);
IMPORT  void    reset_domain_comp(COMPONENT***,RECT_GRID);
IMPORT  void    assign_positive_comp(float (*func)(POINTER,float*),POINTER,
                        COMPONENT***,RECT_GRID,COMPONENT);
IMPORT  void    assign_negative_comp(float (*func)(POINTER,float*),POINTER,
                        COMPONENT***,RECT_GRID,COMPONENT);
IMPORT  void    assign_intersection_comp(float (*func_1)(POINTER,float*),
                        POINTER,float (*func_2)(POINTER,float*),POINTER,
                        COMPONENT***,RECT_GRID,COMPONENT,SIDE,SIDE);
IMPORT  void    make_grid_surfaces(BLK_CRX*,EG_CRX*,int*,bool);
IMPORT  void    alloc_grid_crx_mem(EG_CRX*,int*,int,bool); 
IMPORT  void    free_grid_crx_mem(EG_CRX*,bool);
IMPORT  bool    onfront_block(int,int,int,const EG_CRX*);
IMPORT  bool    is_curve_crx(COMPONENT,COMPONENT,COMPONENT,COMPONENT);
IMPORT  int     install_grid_crx(float (*func)(POINTER,float*),POINTER,
                        EG_CRX*,RECT_GRID,COMPONENT,COMPONENT);
IMPORT  int     count_crx_through_comp(int*,COMPONENT***);
IMPORT  int     make_curves_from_blk(CURVE**,int*,int*,BLK_TRI ****,CURVE **,int);
// make 3 comp surfaces from comp functions
IMPORT	void    show_comp(COMPONENT ***,RECT_GRID);
IMPORT  bool    make_surfaces_from_comp(RECT_GRID*,int (*func)(POINTER, float*),POINTER,
		SURFACE**,CURVE**,int*,int*);

//    iwallsurf.c
IMPORT  void  set_is_wall_surface(SURFACE *);
IMPORT  void  set_is_not_wall_surface(SURFACE *);
IMPORT  bool  is_wall_surface(SURFACE *);
IMPORT  void  get_default_fluid_comp(int*, int*, INTERFACE*);
IMPORT  void  set_wall_flag_for_surface(INTERFACE *);
IMPORT  void  reset_wall_flag_for_surface(INTERFACE *);
IMPORT  void  get_default_fluid_comp(int*, int*, INTERFACE*);
IMPORT  int   add_to_o_surfaces(O_SURFACE **, int*, SURFACE *sp[4]);
IMPORT  void    set_grid_for_surface_construction(RECT_GRID*,RECT_GRID*);

//	itrisset.c
IMPORT	void	set_tol_for_tri_sect(float);
IMPORT	void	swap_positions(float*,float*,int);

IMPORT	bool 	tri_recorded(TRI*,TRI**,int);
IMPORT  bool    two_tris_share_pts(TRI*,TRI*,int);
IMPORT  bool    two_tris_share_side(TRI*,TRI*,int);
IMPORT	bool    point_in_crx_tri(float*,TRI*);
IMPORT	int	merge_tris_set(TRI**,int,TRI**,int);
IMPORT	int	set_tris_set_in_box(TRI**,int,int*,int*,INTERFACE*);
IMPORT	void	move_bound_inside_grid(int*,RECT_GRID*,int);

IMPORT	void	tri_bound_block(float**,TRI*);
IMPORT	void	tris_bound_box(float**,TRI**,int);
IMPORT	bool	blocks_sect(float**,float**);
IMPORT	bool	tris_sect(TRI*,TRI*);

IMPORT	void	plane_of_tri(float*,TRI*);
IMPORT	bool 	plane_side_intersection(const float*,TRI*,int,float*,int*);
IMPORT	bool	test_tris_intersection(TRI*,TRI*);
IMPORT	int	tris_intersection(TRI**,TRI**, int);
IMPORT	bool	tangled_tris_bound_box(int*,int*,int*,int*,INTERFACE*);

IMPORT	bool	link_neighbor_null_side_tris(TRI*,TRI*);
IMPORT	int	linking_tris_with_pairs(TRI**,int,TRI**,int,TRI**,int,
				TRI**,int);
IMPORT  void    centroid_of_tri(float*,TRI*);
IMPORT  bool  	skip_bdry_tri(TRI*);
IMPORT  void    sort_tris_set(TRI**,int,POINTER);
IMPORT  int     bound_tris_set(TRI**,TRI**,int);
IMPORT  int     count_tris_in_box(int*,int*,INTERFACE*);
IMPORT  int     count_tris_in_top_box(int*,int*,INTERFACE*);
IMPORT  int     tris_set_in_top_box(TRI**, int, int*, int*, INTERFACE*);
IMPORT  bool 	tri_in_grid_block(TRI*,int*,int*,RECT_GRID*);
IMPORT  int     rect_boxes_from_tangled_tris(RECT_BOX*,INTERFACE*);
IMPORT  int     remove_single_tri_null_loop(TRI**,int,bool);
IMPORT  int     seal_all_loops_wo_constraint(TRI**,int*,TRI**,int,int,bool);
IMPORT  bool    check_valid_tris(TRI**,int,INTERFACE*);
IMPORT  bool    check_valid_intfc(const char*, INTERFACE*);
IMPORT	int     sep_common_edge_from_tris(TRI***,TRI**,int,INTERFACE*);
IMPORT	bool    sep_common_point_from_loop(TRI **,int,TRI**,int*,INTERFACE*);
IMPORT	void	set_comm_pt_fac(float);
IMPORT	float	get_comm_pt_fac();
IMPORT	void    compute_point_smooth(SMOOTH_PARA*,SMOOTH_TOL*,INTERFACE*);
IMPORT	int	linking_tris_with_pairs_fix(TRI**,int,TRI**,int,TRI**,int,
				TRI**,int,INTERFACE*);

/*      iblkc2.c */
IMPORT  BLK_CRX *alloc_blk_crx(bool);
IMPORT  int     construct_comp2_blk(BLK_CRX*,BLK_TRI*);
IMPORT  void    stitch_inside_blk(BLK_TRI*);
IMPORT  void    stitch_adj_blk(BLK_TRI*,BLK_TRI*);
IMPORT  void    remove_null_pair(BLK_TRI*,BLK_TRI*,int);
IMPORT  void    create_triangle(BLK_TRI*,POINT*,POINT*,POINT*,SURFACE*);
IMPORT  void    reorder_curve_link_list(CURVE*);
IMPORT	ORIENTATION curve_surface_orientation(SURFACE*,CURVE*);
IMPORT  void    set_debug_name(const char*);

/*      iblkc3.c */
IMPORT  int     is_surface(BLK_CRX*,SURFACE*);
IMPORT  int     construct_comp3_blk(BLK_CRX*,BLK_TRI*);

/*	ixgraph.c */
EXPORT  void 	xgraph_2d_intfc_within_range(const char*,INTERFACE*,float*,
				  float);
IMPORT  FILE 	*xgraph_file_open(const char*,const char*,
				  const COORDINATE_PLANE);
IMPORT  void 	xgraph_RECT_GRID(const char*,RECT_GRID*);
IMPORT  void 	xgraph_affine_vector(FILE*,float*,float*,
				     const COORDINATE_PLANE,const char*);
IMPORT  void 	xgraph_curve(FILE*,CURVE*,const COORDINATE_PLANE);
IMPORT  void 	xgraph_interface_curves(const char*,const char*,INTERFACE*,
					const COORDINATE_PLANE);
IMPORT  void  	xgraph_interface_nodes(const char*,const char*,INTERFACE*,
				       const COORDINATE_PLANE);
IMPORT  void 	xgraph_interface_surfaces(const char*,const char*,INTERFACE*,
					  const COORDINATE_PLANE);
IMPORT  void 	xgraph_line_segment(FILE*,float*,float*,
				    const COORDINATE_PLANE,const char*);
IMPORT  void 	xgraph_new_data_set(FILE*);
IMPORT  void 	xgraph_point(FILE*,float*,const COORDINATE_PLANE,const char*);
IMPORT  void 	xgraph_tri(FILE*,TRI*,const COORDINATE_PLANE);
IMPORT  void 	xgraph_tris_list(FILE*,TRI**,int,const COORDINATE_PLANE);
IMPORT  bool    point_within_range(POINT*,float*,float,int);

#if defined(USE_OVERTURE)
/*   intfc_amr.c  */

IMPORT void set_amr_intfc_tol(INTERFACE*,float); 

#endif /* if defined(USE_OVERTURE) */

#endif /* !defined(_IPROTOS_H) */
