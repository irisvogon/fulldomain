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
*				ilocprotos.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#if !defined(_ILOCPROTOS_H)
#define _ILOCPROTOS_H


 /* libintfc.a EXPORTED Function Declarations Internal to Interface*/

/*	comp.c*/

/*	iecomps.c*/
IMPORT	POINTER	copy_e_comps(INTERFACE*);

/*	iblkb.c*/
IMPORT	int construct_bdry_blk(BLK_CRX*,BLK_TRI*);

/*	igrid.c*/

/*	intfc.c*/
IMPORT	CURVE	*i_read_curve(INTERFACE*,int);
IMPORT	POINT	*i_read_print_point(INTERFACE*,const IO_TYPE*,bool);
IMPORT	NODE	*i_read_node(INTERFACE*,int);
IMPORT	POINT	*i_read_point(INTERFACE*,int);
IMPORT	bool	next_hypersurface2d(INTERFACE*,HYPER_SURF**);
IMPORT  bool next_point1d(INTERFACE*,POINT**,
			  HYPER_SURF_ELEMENT**,HYPER_SURF**);
IMPORT  bool next_point2d(INTERFACE*,POINT**,
		       HYPER_SURF_ELEMENT**,HYPER_SURF**);
IMPORT	int	i_read_print_intfc_rect_grids(const IO_TYPE*,INTERFACE*,REMAP*);
IMPORT	int	winding_number(float*,float*,int);
IMPORT	struct	Chunk	*new_chunk(struct Table*);
IMPORT	void	i_fprint_boundary_type(FILE*,const char*,int,const char*,
                                       INTERFACE*);
IMPORT	void	i_fprint_curve(FILE*,CURVE*);
IMPORT	void	i_fprint_intfc_rect_grids(FILE*,INTERFACE*);
IMPORT	void	i_fprint_point(FILE*,POINT*);
IMPORT	void	i_fprint_node(FILE*,NODE*);
IMPORT	void	print_table(struct Table*);
IMPORT	void	print_table_list(void);

/*	irefl.c*/
IMPORT	void	i_reflect_interface1d(INTERFACE*,float*,float*);
IMPORT	void	i_reflect_interface2d(INTERFACE*,float*,float*);
IMPORT	void	i_reflect_interface3d(INTERFACE*,float*,float*);

/*	iscatter.c*/

/*	isub.c*/
IMPORT	int	i_print_number_of_tangles(const char*,INTERFACE*,CROSS*);
IMPORT	void	set_tri_storage_type(TRI_STORAGE_TYPE type);

/*	iuserintfc.c*/

/*	ppcopy.c*/
IMPORT	void	i_reconstruct_surface_pointers(SURFACE*,INTERFACE*,INTERFACE*,
					       POINTER*,POINTER*,int);
IMPORT	void	i_reconstruct_tri_pointers(TRI*,INTERFACE*,INTERFACE*,
					   POINTER*,POINTER*,int);

/*	shift.c*/

/*	top.c*/
IMPORT	CURVE	*adjacent_curve_along_surface(CURVE*,ORIENTATION,
					      ANGLE_DIRECTION,
					      ORIENTATION*,float*,int);
IMPORT	bool	i_is_subdomain_boundary(HYPER_SURF*);
IMPORT	bool	i_is_subdomain_node(NODE*);
IMPORT	bool	i_is_virtual_fixed_node(NODE*);
IMPORT	void	fprint_bond(FILE*,BOND*);
IMPORT	void	set_curves_at_onode(O_NODE*);

/*	zoom.c*/

/*	comp1d.c*/
IMPORT	COMPONENT	component1d(float*,INTERFACE*);
IMPORT	bool	i_intersections1d(INTERFACE*,CROSS**,const bool);
IMPORT	bool	long_nearest_interface_point1d(float*,COMPONENT,INTERFACE*,
					       USE_BOUNDARIES,HYPER_SURF*,
					       float*,float*,
					       HYPER_SURF_ELEMENT**,
					       HYPER_SURF**);
IMPORT	bool	nearest_interface_point1d(float*,COMPONENT,INTERFACE*,
					  USE_BOUNDARIES,HYPER_SURF*,float*,
					  float*,HYPER_SURF_ELEMENT**,
					  HYPER_SURF**);
IMPORT	void	i_print_intersections1d(CROSS*,INTERFACE*);
IMPORT void	show_COMP_1d(FILE*,INTERFACE*);
IMPORT	void	show_point_comp_lists(INTERFACE*);


/*	comp2d.c*/
IMPORT	COMPONENT	component2d(float*,INTERFACE*);
IMPORT	COMPONENT	long_component2d(float*,INTERFACE*);
IMPORT	bool	long_nearest_interface_point2d(float*,COMPONENT,INTERFACE*,
					       USE_BOUNDARIES,HYPER_SURF*,
					       float*,float*,
					       HYPER_SURF_ELEMENT**,
					       HYPER_SURF**);
IMPORT	bool	long_nearest_similar_interface_point2d(float*,COMPONENT,
				                       COMPONENT,INTERFACE*,
						       USE_BOUNDARIES,
						       HYPER_SURF*,float*,
						       float*,
						       HYPER_SURF_ELEMENT**,
						       HYPER_SURF**);
IMPORT	bool	nearest_interface_point2d(float*,COMPONENT,INTERFACE*,
					  USE_BOUNDARIES,HYPER_SURF*,float*,
					  float*,HYPER_SURF_ELEMENT**,
					  HYPER_SURF**);
IMPORT	bool	nearest_similar_interface_point2d(float*,COMPONENT,COMPONENT,
						  INTERFACE*,USE_BOUNDARIES,
						  HYPER_SURF*,float*,float*,
						  HYPER_SURF_ELEMENT**,
						  HYPER_SURF**);
IMPORT	bool	nearest_interface_point2d(float*,COMPONENT,INTERFACE*,
					  USE_BOUNDARIES,HYPER_SURF*,float*,
					  float*,HYPER_SURF_ELEMENT**,
					  HYPER_SURF**);
IMPORT	void	show_COMP_2d(FILE*,INTERFACE*);

/*	cross2d.c*/

/*	isect2d.c*/
IMPORT	bool	i_intersections2d(INTERFACE*,CROSS**,const bool);
IMPORT	void	i_print_crossing_elements2d(CROSS*,INTERFACE*);
IMPORT	void	i_print_intersections2d(CROSS*,INTERFACE*);

/*	setb2d.c*/


/*	comp3d.c*/
IMPORT	COMPONENT	component3d(float*,INTERFACE*);
IMPORT	COMPONENT	long_component3d(float*,INTERFACE*);
IMPORT	COMPONENT	dir_long_component3d(float*,INTERFACE*);
IMPORT	bool	make_tri_lists(INTERFACE*);
IMPORT	bool	long_nearest_interface_point3d(float*,COMPONENT,INTERFACE*,
					       USE_BOUNDARIES,HYPER_SURF*,
					       float*,float*,
					       HYPER_SURF_ELEMENT**,
					       HYPER_SURF**);
IMPORT	bool	long_nearest_similar_interface_point3d(float*,COMPONENT,
						       COMPONENT,INTERFACE*,
						       USE_BOUNDARIES,
						       HYPER_SURF*,float*,
						       float*,
						       HYPER_SURF_ELEMENT**,
						       HYPER_SURF**);
IMPORT	bool	nearest_interface_point3d(float*,COMPONENT,INTERFACE*,
					  USE_BOUNDARIES,HYPER_SURF*,float*,
					  float*,HYPER_SURF_ELEMENT**,
					  HYPER_SURF**);
IMPORT	bool	nearest_similar_interface_point3d(float*,COMPONENT,COMPONENT,
						  INTERFACE*,USE_BOUNDARIES,
						  HYPER_SURF*,float*,float*,
						  HYPER_SURF_ELEMENT**,
						  HYPER_SURF**);
IMPORT	void	show_COMP_3d(FILE*,INTERFACE*);

/*	int3d.c*/
IMPORT	SURFACE *read_print_surface(INTERFACE*,const IO_TYPE*,int,
                                    INTERFACE_ADDRESSES*,bool);
IMPORT	bool	next_hypersurface3d(INTERFACE*,HYPER_SURF**);
IMPORT	bool	next_point3d(INTERFACE*,POINT**,HYPER_SURF_ELEMENT**,
			     HYPER_SURF**);
IMPORT	void	copy_all_surfaces(INTERFACE*,INTERFACE*);
IMPORT	void	fprint_tri(FILE*,TRI*,INTERFACE*);
IMPORT	void	fprint_tris_on_curve(FILE*,CURVE*);
IMPORT	void	i_fprint_surface(FILE*,SURFACE*);
IMPORT	void	print_tris_on_surface(SURFACE*);
IMPORT	void	read_print_tris_on_curve(FILE*,CURVE*,int,
					 INTERFACE_ADDRESSES*);
IMPORT	void	set_tri_array_numbers(INTERFACE*,int);

/*	iprt3d.c*/
IMPORT	void	threed_interface_plot(INTERFACE*,RECT_GRID*);

/*	isect3d.c*/
IMPORT	bool	i_intersections3d(INTERFACE*,CROSS**,const bool);
IMPORT	void	i_print_crossing_elements3d(CROSS*,INTERFACE*);

/*	map.c*/
IMPORT	bool	is_interior_vertex(TRI*,POINT*,TRI**,INTERFACE*);

/*	setb3d.c*/


/*	icheck3d.c*/
IMPORT void	check_tri_and_bond(TRI*,BOND*,const char*,INTERFACE*);

/*	ineighborhood.c */

#endif /* !defined(_ILOCPROTOS_H) */
