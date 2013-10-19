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
*				girgb.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#include <ginit/ginit.h>

 
EXPORT	void 	init_fluid_rigid_body(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	LAYER_FLAG      flag;
	LAYER_SYS       *layer_sys;
	LAYER           **layer;
	Front           *front = ip->root->front;
	RECT_GRID	*gr = front->rect_grid;
	INTERFACE	*intfc = front->interf;
	COMPONENT	left_c,right_c;
	RIGID_BODY_PARAMS rgb_params;
	CURVE		*curve;
	POINTER		func_params;
	float		(*func)(POINTER,float*);
	int 		i, dim = gr->dim;;
	int		num_rigid_bodys;	
	char		mesg[100];
	char		s[Gets_BUF_SIZE];

	debug_print("rigid_body","Entered init_fluid_rigid_body()\n");

	scalar(&layer_sys,sizeof(LAYER_SYS));

	screen("Enter yes to initialize a shock in the region: ");
	(void) Gets(s);
	if (s[0] != 'y' || s[0] != 'Y')
	    layer_sys->num_layers = 2;
	else
	    layer_sys->num_layers = 1;
	uni_array(&layer,5*layer_sys->num_layers+1,sizeof(LAYER *));
	layer_sys->layer = layer;
	layer_sys->front = front;
	CLEAR_LAYER_FLAG(layer_sys->flag);

	for (i = 1; i <= layer_sys->num_layers; ++i)
	{
	    layer[i] = prompt_for_each_layer(i,layer_sys->num_layers,
	    			&layer_sys->flag,front,layer,ip);
	}
	if (layer_sys->num_layers == 2)
	{
	    layer[1]->upper_surf->r_comp = layer[1]->comp;
	    layer[1]->upper_surf->l_comp = layer[2]->comp;
	    make_layer_surf(front,layer[1]->upper_surf,init);
	    right_c = layer[2]->comp;
	}
	else
	    right_c = layer[1]->comp;
	left_c = COMPOBST;
	set_obstacle_comp_type(comp_type(COMPOBST),front);

	/** Start specification of rigid bodies **/

	screen("Enter number of rigid bodys: ");
	Scanf("%d\n",&num_rigid_bodys);
	screen("Available geometrical shapes of rigid body are\n");
	screen("\tEllipse (E)\n");
	for (i = 1; i <= num_rigid_bodys; ++i)
	{
	    screen("Enter the shape of rigid body %d: ");
	    (void) Gets(s);
	    switch (s[0])
	    {
	    case 'e':
	    case 'E':
	    	func_params = init_ellipse_params(gr);
		func = ellipse_func;
	    	break;
	    }
	    curve = make_level_curve(gr,intfc,left_c,right_c,func,
	    			func_params,NO);
	    if (curve == NULL)
	    {
	    	screen("Cannot make rigid body %d\n",i);
		clean_up(ERROR);
	    }
	    sprintf(mesg,"the %d-th rigid body",i);
	    prompt_for_rigid_body_params(&rgb_params,mesg,dim);

	    start_status(curve) = end_status(curve) = INCIDENT;
	    wave_type(curve) = MOVABLE_NEUMANN_BOUNDARY;
	    total_mass(curve) = rgb_params.total_mass;
	    mom_inertial(curve) = rgb_params.mom_of_inertial;
	    angular_velo(curve) = 0.0;
	    surface_tension(curve) = 0.0;
	    for (i = 0; i < dim; ++i)
	    {
	    	center_of_mass(curve)[i] = rgb_params.center_of_mass[i];
		center_of_mass_velo(curve)[i] = 0.0;
	    }
	}

	/** End specification of rigid bodies **/

	init_comp_type(init,ip,layer_sys);
	if (debugging("rigid_body"))
	{
	    gview_plot_interface("out",intfc);
	    make_interface_topology_lists(intfc);
	    show_COMP(stdout,intfc);
	}
	debug_print("rigid_body","Left init_fluid_rigid_body()\n");
}	/* end init_fluid_rigid_body */

