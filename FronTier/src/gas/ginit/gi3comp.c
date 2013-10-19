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
*				gi3comp.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Initialized an injection jet.
*/

#define DEBUG_STRING    "i3comp"
#include <ginit/ginit.h>

	/* LOCAL Function Declarations */
LOCAL  	float   jet_inner_func(POINTER,float*);
LOCAL  	float   jet_outer_func(POINTER,float*);
LOCAL   int     fuel_comp_func(POINTER, float*);
LOCAL   int     flowin_comp_func(POINTER, float*);
LOCAL   int     flowin_comp_func_x(POINTER, float*);
LOCAL   int     cyl_pert_comp_func(POINTER, float*);
LOCAL   int     plane_pert_comp_func(POINTER, float*);
LOCAL   void    init_types_from_comp(POINTER);
LOCAL   int     compare_comps(const void *, const void *);
LOCAL   int     init_surface_type_from_comp(SURFACE *);
LOCAL   int     insert_comp(int*,int,int);
LOCAL   int     init_curve_type_from_comp(SURFACE *s[]);
LOCAL   void    init_wave_curve_types(INTERFACE*,SURFACE**,CURVE**,int,int);
LOCAL   void    init_rect_boundary(Front*);
LOCAL   void    init_bdry_surfs_flag(INTERFACE*,POINTER,INIT_DATA*,INIT_PHYSICS*);

int	RI_comp_func(void*,float*);
void	RI_invicid_solution_outer(float*,float*,float*,float*,
		float,float,float,float,RI_PERT*);
void	RI_invicid_solution_inner(float*,float*,float*,float*,
		float,float,float,float,RI_PERT*);
void	RI_invicid_growth_rate(RI_PERT*,float,float);


typedef struct {
        float *cen;             // Center of buttom of the jet 
        float r1;               // Lengths of radii in upper section 
        float r2;               // Lengths of radii in lower section 
	float *h;               // height of different sections 
        float ratio;            // ratio of height from contact 
	                        // surface to top of the nozzle 
} JET_PARAMS;

typedef struct {
	JET_PARAMS  *params_in, *params_out;
	BDRY_BOX_PARAMS  *bp;
	int     comp_top, comp_obst, comp_fluid, comp_gas, comp_bot; 
	float	**file_r;
        float   hz, r, *cen;
} FUEL_PARAMS;

LOCAL	RI_PERT		ri_param;

LOCAL  int  surf_type[10][4], curv_type[10][4];
LOCAL  int  snum, cnum;

LOCAL  void  init_types_from_comp(
	POINTER func_params)
{
FUEL_PARAMS   *params = (FUEL_PARAMS*)func_params;
int  comp_top = params->comp_top, comp_bot = params->comp_bot,
     comp_obst = params->comp_obst, comp_fluid = params->comp_fluid,
     comp_gas = params->comp_gas;

#define RK_USE_FLOW 1
int  st[][4] = { {comp_top, comp_bot, PASSIVE_BOUNDARY, 1},
		 {comp_top, comp_obst, PASSIVE_BOUNDARY, 1},
		 {comp_bot, comp_obst, PASSIVE_BOUNDARY, 1},
		 {comp_top, comp_gas, DIRICHLET_BOUNDARY, 1},
		 {comp_obst,  comp_gas,   NEUMANN_BOUNDARY, 0},
#if RK_USE_FLOW
		 {comp_bot,   comp_gas,   FLOW_BOUNDARY, 1},
#else
		 {comp_bot,   comp_gas,   DIRICHLET_BOUNDARY, 1},
#endif
		 {comp_fluid, comp_gas,   CONTACT, 0} };

int  sn = 7;
int  ct[][4] = { {comp_top, comp_bot, comp_obst, FIXED_CURVE},
		 {comp_bot, comp_gas, comp_obst, FIXED_CURVE},
		 {comp_top, comp_gas, comp_obst, FIXED_CURVE} };
int  i, cn = 3;

	snum = sn;
	cnum = cn;

	//#bjet2 for lgb
	//if(debugging("jet3d_lgb"))
	    ct[1][3] = FIXED_CURVE;

	for(i=0; i<sn; i++)
	    ft_assign(surf_type[i], st[i], 4*INT);
	for(i=0; i<cn; i++)
	    ft_assign(curv_type[i], ct[i], 4*INT);
}

LOCAL  int  compare_comps(const void *a, const void *b)
{
int  *c1 = (int*)a, *c2 = (int*)b;
	return  *c1 < *c2 ? -1 : 1;
}

LOCAL  int  insert_comp(int *carr, int n, int comp)
{
int	i;

	for(i=0; i<n; i++)
	    if(carr[i] == comp)
	        return n;

	carr[i] = comp;
	return n+1;
}

LOCAL  int  init_surface_type_from_comp(SURFACE *s)
{
int    nc = negative_component(s), pc = positive_component(s);
int    i;

	for(i=0; i<snum; i++)
	{
	    if( (pc == surf_type[i][0] && nc == surf_type[i][1]) ||
	        (pc == surf_type[i][1] && nc == surf_type[i][0]) )
	    {
		//printf("sc  %d  %d    %d\n", surf_type[i][0], surf_type[i][1], i);
	        return i;
	    }
	}
	printf("ERROR init_surface_type_from_comp, no type for such comps %d %d\n", nc, pc);
	clean_up(ERROR);
}

LOCAL  int  init_curve_type_from_comp(SURFACE  *s[])
{
int    comp[6], cc[3], i, n=0;

	for(i=0; i<3; i++)
	{
	    n = insert_comp(comp, n, positive_component(s[i]));
	    n = insert_comp(comp, n, negative_component(s[i]));
	}
	if(n != 3)
	{
	    printf("ERROR init_curve_type_from_comp, comps in surfs are not adjacent.\n");
	    clean_up(ERROR);
	}
	qsort((POINTER)comp, 3, INT, compare_comps);
	for(i=0; i<cnum; i++)
	{
	    ft_assign(cc, curv_type[i], 3*INT);
	    qsort((POINTER)cc, 3, INT, compare_comps);
	    if(cc[0] == comp[0] && cc[1] == comp[1] && cc[2] == comp[2])
	    {
	        //print_int_vector("cc", cc, 3, "  ");
		//printf(" %d \n", i);
	        return i;
	    }
	}
	
	printf("ERROR init_curve_type_from_comp, no type for such comps %d %d %d\n", 
		    comp[0], comp[1], comp[2]);
	clean_up(ERROR);
}

LOCAL  void  init_wave_curve_types(
	INTERFACE    *intfc,
	SURFACE      **s,
	CURVE	     **c,
	int	     is,
	int	     ic)
{
int	  i,j, n;
BOND_TRI  *btri;
SURFACE   *sc[3];
float	  surf_ten;

	for (i = 0; i < ic; i++)
	{
	    for(j=0; j<3; j++)
	    {
	        btri = Btris(c[i]->first)[j];
	        sc[j] = btri->surface;
	    }
	    n = c[i]->number = init_curve_type_from_comp(sc);
	    curve_type(c[i]) = curv_type[n][3];
	}
	

	if(debugging("riinit"))
	{
	    surf_ten = ri_param.sigma;
	}
	else
	{
	    surf_ten =  prompt_for_surface_tension(CONTACT, "for the jet");
	}

	for (i = 0; i < is; i++)
	{
	    n = s[i]->number = init_surface_type_from_comp(s[i]);
	    wave_type(s[i]) = surf_type[n][2];
	    if(surf_type[n][3] == 1)
	        set_is_bdry(s[i]);
	    if(wave_type(s[i]) == CONTACT)
	    {
	        surface_tension(s[i]) = surf_ten;
		printf("#surf_ten = %15.8e\n", surf_ten);
	    }
	}

}
LOCAL  void   init_rect_boundary(
	Front   *front)
{
PP_GRID    *pp_grid = front->pp_grid;
INTERFACE  *intfc = front->interf;
int	   i, j, dim = intfc->dim;
int        me[3], *G;


	G = pp_grid->gmax;
	find_Cartesian_coordinates(pp_mynode(),pp_grid,me);

	//TMP fix for jet
	
	if(me[0] == G[0] - 1)
	    //rect_boundary_type(intfc,0,1) = OPEN_BOUNDARY;
	    rect_boundary_type(intfc,0,1) = REFLECTION_BOUNDARY;
	if(me[0] == 0)
	    rect_boundary_type(intfc,0,0) = MIXED_TYPE_BOUNDARY;
//	    rect_boundary_type(intfc,0,0) = REFLECTION_BOUNDARY; /* nozzle YYU */
	if(me[1] == G[1] - 1)
//	    rect_boundary_type(intfc,1,1) = REFLECTION_BOUNDARY; //TMP_TK
	    rect_boundary_type(intfc,1,1) = MIXED_TYPE_BOUNDARY;
	if(me[1] == 0)
//	    rect_boundary_type(intfc,1,0) = REFLECTION_BOUNDARY; //TMP_TK
	    rect_boundary_type(intfc,1,0) = MIXED_TYPE_BOUNDARY;
	if(me[2] == G[2] - 1)
	    //rect_boundary_type(intfc,2,1) = REFLECTION_BOUNDARY;/* 100325 YYU*/
	    rect_boundary_type(intfc,2,1) = MIXED_TYPE_BOUNDARY;
	    //rect_boundary_type(intfc,2,1) = OPEN_BOUNDARY;
	if(me[2] == 0)
	    //rect_boundary_type(intfc,2,0) = OPEN_BOUNDARY;
	    rect_boundary_type(intfc,2,0) = PASSIVE_BOUNDARY;
	
	for (i = 0; i < dim; i++)
	{
	    for (j = 0; j < 2; j++)
	    {
		if (((me[i]>0) && (j==0)) || ((me[i]<(G[i]-1)) && (j==1)))
		    rect_boundary_type(intfc,i,j) = SUBDOMAIN_BOUNDARY;
	    }
	}
}

LOCAL  void    init_bdry_surfs_flag(
	INTERFACE	  *intfc,
	POINTER           func_params,	
	INIT_DATA	  *init,
	INIT_PHYSICS	  *ip)
{            
FUEL_PARAMS   *params = (FUEL_PARAMS *)func_params;
float         *tcen = params->cen, bcen[3];
int           bst_index1, bst_index2;
SURFACE	      **s;

	bcen[0] = params->bp->L[0];
	bcen[1] = tcen[1];
	bcen[2] = params->bp->L[2];
	
	printf("Input Dirichlet boundary condtion for top fluid %f %f %f\n",tcen[0],tcen[1],tcen[2]);
	bst_index1 = g_prompt_for_boundary_state(DIRICHLET_BOUNDARY,NULL,tcen,
                        params->comp_gas,-1,NULL,init,ip);
	
	printf("Input Dirichlet boundary condtion for bottom gas %f %f %f\n",bcen[0],bcen[1],bcen[2]);
	bst_index2 = g_prompt_for_boundary_state(DIRICHLET_BOUNDARY,NULL,bcen,
                        params->comp_gas,-1,NULL,init,ip);

	for (s = intfc->surfaces; s && *s; s++)
	{
	    printf("#s %d (%d %d) %d\n", *s, negative_component(*s), 
	        positive_component(*s), (*s)->number);    
       
	    //init_types_from_comp
	    switch((*s)->number)
	    {
	    case 5:
	        printf("#bst set %d %d\n", *s,  bst_index2);
	        bstate_index(*s) = bst_index2;
		break;
	    case 3:
	        printf("#bst set %d %d\n", *s,  bst_index1);
	        bstate_index(*s) = bst_index1;
		break;
	    }
	}
}


int	RI_init_params(
	int		comp_gas,
	int		comp_fluid,
	BDRY_BOX_PARAMS *bp)
{
	float		sigma,a,lambda;
	float		*L = bp->L, *U = bp->U;
	COMP_TYPE	*ct;
	Locstate	st;
	int		i;

	//Assumption here
	//diameter of the cylinder = domain width/6
	//a = (U[0] - L[0])/3.0/2.0;
	screen("Enter the radius of the jet: ");
	(void)  Scanf("%f\n",&a);
	
	//surf ten coef = 0.01*a
	//sigma = 1.0*a;
	screen("Enter the surface tension of the jet: ");
	(void)  Scanf("%f\n",&sigma);
	
	//amplitude is 12% the radius
	//ri_param.amp = 0.12*a;
	screen("Enter the amplitude of the perturbation: ");
	(void)  Scanf("%f\n",&(ri_param.amp));
	
	//wave length = domain height
	lambda = U[2] - L[2];

	//gas speed is 0
	ri_param.U = 0.0;
	
	ri_param.a = a;
	ri_param.kx = 2.0*PI/lambda;
	ri_param.R = U[0];
	for(i=0; i<3; i++)
	    ri_param.cen[i] = 0.5*(L[i] + U[i]);

	if(debugging("quadom"))
	{
	    ri_param.cen[0] = 0.0;
	    ri_param.cen[1] = 0.0;
	}

	ct = comp_type(comp_fluid);
	st = Ambient(ct);
	ri_param.dens1 = density(st);
	
	ct = comp_type(comp_gas);
	st = Ambient(ct);
	ri_param.dens2 = density(st);
	
	ri_param.press = pressure(st);
	ri_param.sigma = sigma;

	RI_invicid_growth_rate(&ri_param, ri_param.dens1, sigma);

	{
	    printf("\n\n   Perturbed cylinder parameters:\n");
	    printf("a      =  %15.8e\n", ri_param.a);
	    printf("kx     =  %15.8e\n", ri_param.kx);
	    printf("lambda =  %15.8e\n", 2*PI/ri_param.kx);
	    printf("U      =  %15.8e\n", ri_param.U);
	    print_general_vector("cen    = ", ri_param.cen, 3, "\n");
	    printf("sigma  =  %15.8e (surface tension coef)\n", ri_param.sigma);
	    printf("dens1  =  %15.8e\n", ri_param.dens1);
	    printf("dens2  =  %15.8e\n", ri_param.dens2);
	    printf("press  =  %15.8e\n", ri_param.press);
	    printf("beta   =  %15.8e (growth rate)\n", ri_param.beta);
	    printf("\n");
	}

	return 0;
}

//assume all the params in ri_param
//st already has a state

void	RI_initializer(Locstate, float*, int);

void	RI_initializer(
	Locstate	st,
	float		*coords,
	int		comp)
{
	float	lambda, z, r, pos, u, v, p, d, pres, dens, theta;
	int	st_type;

	r = sqrt(sqr(coords[0] - ri_param.cen[0]) 
	       + sqr(coords[1] - ri_param.cen[1]));
	theta = angle(coords[0], coords[1]);
	
	lambda = 2.0*PI/ri_param.kx;
	z = coords[2] + lambda/4.0;

	pos = ri_param.a + ri_param.amp*sin(ri_param.kx*z);

	if(comp == 2)
	{
	    if(NO && r >= pos)
	    {
		printf("ERROR RI_initializer, outer point has comp 2\n");
		printf("r = %15.8e  pos = %15.8e\n", r, pos);
		clean_up(ERROR);
	    }

	    dens = ri_param.dens1;
	    pres = ri_param.press + ri_param.sigma/ri_param.a;

	    RI_invicid_solution_inner(&u,&v,&p,&d,r,z,0.0, dens, &ri_param);
	}
	else
	{
	    if(NO && r <= pos)
	    {
		printf("ERROR RI_initializer, outer point has comp 3\n");
		printf("r = %15.8e  pos = %15.8e\n", r, pos);
		clean_up(ERROR);
	    }

	    dens = ri_param.dens2;
	    pres = ri_param.press;

	    RI_invicid_solution_outer(&u,&v,&p,&d,r,z,0.0, dens, &ri_param);
	}
	
	st_type = state_type(st);
	set_state(st, TGAS_STATE, st);

	Dens(st) = dens + d;
	Press(st) = pres + p;
	Vel(st)[0] = v*cos(theta);
	Vel(st)[1] = v*sin(theta);
	Vel(st)[2] = u;

	set_state(st, st_type, st);
}

void	shift_intfc_in_grid(
		INTERFACE	*intfc,
		RECT_GRID	*gr)
{
        HYPER_SURF              *hs;
        HYPER_SURF_ELEMENT      *hse;
        POINT                   *p;
        float			zshf;

	DEBUG_ENTER(shift_intfc_in_grid)

	zshf = 0.01*gr->h[2];

	(void) next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
	    if(Boundary_point(p))
	        continue;
	    
	    if(wave_type(hs) == CONTACT) 
	    {
		Coords(p)[2] -= zshf;
	    }
	}

	reset_normal_on_intfc(intfc);

	DEBUG_LEAVE(shift_intfc_in_grid)
}

LOCAL	void	move_ref_interface(
	INTERFACE *intfc)
{
	POINT   *p;
	SURFACE **s;
	TRI     *t;
	int     i;

	(void) next_point(intfc,NULL,NULL,NULL);
	
	for (s = intfc->surfaces; s && *s; ++s)
	{
	    for (t = first_tri(*s); !at_end_of_tri_list(t,*s); t = t->next)
	    {
		for (i = 0; i < 3; ++i)
		{
		    p = Point_of_tri(t)[i];
		    if (sorted(p) == NO && 
		        Coords(p)[2] > 4.0-0.012 && Coords(p)[2] < 4.0)
		    {
	                Coords(p)[2] -= 0.012/3;
		    }
		    sorted(p) = YES;
		}
	    }
	}

	interface_reconstructed(intfc) = NO;
}



bool make_surfaces_from_files(RECT_GRID*,SURFACE**,CURVE**,int*,int*);
LOCAL   void   get_file_r(float**, int, int);

/*
*		init_fuel_injection_jet3d():
*
*	Initializes a gas chamber for the injection jet.
*
*
*       |           
*       |                                      
*       |-------------- r_n ---->|______________________________________ 
*	|    -                   |######################################|
*	|    |                   |######################################|
*	|    |          COMPB    |######################################|
*       |    |                   |######################################|
*       |   h_u                  |######################################|
*       |____|___contact_________|##################### COMPOBST #######|
*       |    |                   |######################################|
*       |    |                   |######################################|
*	|    |                   |######################################|
*       |    -                -         \###############################|
*       |                     |               \#########################|  
*	|                     |                   \#####################|  
*	|                     |                     \###################|  
*	|                    h_l                      \#################|  
*	|                     |                        |################|  
*	|                     |                         |###############|  
*	|           COMPA     |                          \##############|  
*	|                     -                          |##############|  -
*	|                                                |##############|  |
*	|                                                |##############|  |
*	|                                                |##############|  |
*	|                                                |##############|  |
*	|                                                |##############|  |
*	|                                                |##############|  h_i
*	|                                                |##############|  |
*	|                                                |##############|  |
*	|                                                |##############|  |
*	|<-------------------------r_i------------------>|##############|  |
*	|                                                |##############|  |
*	|                                                |##############|  |
*	|                                                |<-----w_w---->|  |
*	|                                                |##############|  |
*	|________________________________________________|##############|  -
*
*/

EXPORT	void init_3comp_jet3d(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	Front           *front = ip->root->front;
	INTERFACE       *intfc = front->interf;
	RECT_GRID       *gr = front->rect_grid;
	CURVE           **jet;
	SURFACE	        **surfs ;
	JET_PARAMS      jet_paramsa, jet_paramsb;
	BDRY_BOX_PARAMS bp;
	FUEL_PARAMS     fuel_params;
	Gas_param       *gas_paramsa, *gas_paramsb,  *gas_paramsc;
	float           *L = gr->GL;
	float           *U = gr->GU;
	float	        *h = gr->h;
	float	        coords[MAXD];
	char 	        s[1024];
	int             i, dim = gr->dim;
	float           r_i, w_w, h_i, r_n, h_l, h_u;
	float           cen[3], h_in[3], h_out[3], ratio;
	float           tol = MIN_SC_SEP(intfc);
        int             is,ic;

	debug_print("jet","Entered init_3comp_jet3d()\n");

	//#bjet2 TMP fix
	//make_bubbles(NULL, front);

	set_obstacle_comp_type(comp_type(COMPOBST),front);
	exclude_comp(COMPD,intfc);
	set_obstacle_comp_type(comp_type(COMPD),front);

	/* read geometrics of the jet */
	//TMP comment for out crx case
	
	screen("\nFuel injection jet enters from below.\n");
	screen("Please enter the geometric parameters of the fuel jet.\n");
	
	screen("Enter the center values of the floor of fuel container: ");
	(void) Scanf("%f %f %f\n",&cen[0],&cen[1],&cen[2]);
	if ((L[0] > cen[0]) || (U[0] < cen[0]) ||
	    (L[1] > cen[1]) || (U[1] < cen[1]) ||
	    (L[2] > cen[2]) || (U[2] < cen[2]))
	{
	    screen("ERROR in init_fuel_injection_jet3d(), "
	       "invalid center of the floor x0 = %"FFMT", "
	       "y0 = %"FFMT", z0 = %"FFMT"\n",cen[0],cen[1],cen[2]);
	    (void) printf("L[0] = %"FFMT", U[0] = %"FFMT"\n",L[0],U[0]);
	    (void) printf("L[1] = %"FFMT", U[1] = %"FFMT"\n",L[1],U[1]);
	    (void) printf("L[2] = %"FFMT", U[2] = %"FFMT"\n",L[2],U[2]);
	    //clean_up(ERROR);
	}

	screen("Enter the inner radius of the fuel container: ");
	(void) Scanf("%f\n",&r_i);
	if ((r_i <= 0.0) || (U[0] < (L[0] + r_i)))
	{
	    screen("ERROR in init_fuel_injection_jet3d(), "
	       "invalid inner radius r_i = %"FFMT"\n",r_i);
	    (void) printf("L[0] = %"FFMT", U[0] = %"FFMT"\n",L[0],U[0]);
	    //clean_up(ERROR);
	}

	screen("Enter the width of the fuel container wall: ");
        (void) Scanf("%f\n",&w_w);
	if (w_w <= 0.0)
	{
	    screen("ERROR in init_fuel_injection_jet3d(), "
		       "invalid container width = %"FFMT"\n",w_w);
	    //clean_up(ERROR);
	}
   
	screen("Enter the radium of the nozzle: ");
	(void) Scanf("%f\n",&r_n);
	if ((r_n <= 0.0) || (r_i < r_n) || (U[0] < (L[0] + r_n)))
	{
	    screen("ERROR in init_fuel_injection_jet3d(), "
	       "invalid nozzle radium r_n = %"FFMT"\n",r_n);
	    (void) printf("L[0] = %"FFMT", U[0] = %"FFMT"\n",L[0],U[0]);
	    //clean_up(ERROR);
	}
   
	h_i = 0.0;
	screen("Enter the height of the interior vertical wall\n"
	   "\tof the fuel chamber (default = %g): ",h_i);
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    sscan_float(s,&h_i);
	    if ((h_i < 0.0) || (U[2] < (L[2] + h_i)))
	    {
	        screen("ERROR in init_fuel_injection_jet3d(), "
		           "invalid interior wall height h_i = %"FFMT"\n",h_i);
	        (void) printf("L[2] = %"FFMT", U[2] = %"FFMT"\n",L[2],U[2]);
	        //clean_up(ERROR);
	    }
	}

	h_l = 0.0;
	screen("Enter the height of the curve section of the nozzle\n"
	       "\t(default or zero = no lower section): ");
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    if ((sscan_float(s,&h_l) != 1) || (h_l < 0.0) ||
	        (U[2] < (L[2] + h_i + h_l) ))
	    {
	        screen("ERROR in init_fuel_injection_jet3d(), "
	           "invalid curve section nozzle "
		   "height nh_l = %"FFMT"\n",h_l);
	        (void) printf("L[2] = %"FFMT", U[2] = %"FFMT"\n",L[2],U[2]);
	        //clean_up(ERROR);
	    }
	}

	h_u = 0.0;
	screen("Enter the height of the cylinder section of the nozzle\n"
	           "\t(default or zero = no cylinder section): ");
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    if ((sscan_float(s,&h_u) != 1) || (h_u < 0.0) ||
	        (U[2] < (L[2] + h_i + h_l + h_u) ))
	    {
	        screen("ERROR in init_fuel_injection_jet3d(), "
	           "invalid cylinder section nozzle "
		   "height nh_u = %"FFMT"\n",
	           h_u);
	        (void) printf("L[2] = %"FFMT", U[2] = %"FFMT"\n",L[2],U[2]);
	        //clean_up(ERROR);
	    }
	}
	
	ratio = 0.5;
	screen("Enter the ratio of the height of the contace section\n"
	       "\t to the top of the nozzle (default = %f): ",ratio);
	(void) Scanf("%f\n",&ratio);
	if ((ratio > 1.0))
	{
	    screen("ERROR in init_fuel_injection_jet3d(), "
	       "invalid ratio of contace in nozzle = %"FFMT"\n",ratio);
	    clean_up(ERROR);
	}
	//TMP end comment
	
	g_prompt_for_dynamic_bubble_insertion3d(ip->root);

	(void) prompt_for_eos_params(init,ip,YES,"");
	gas_paramsb = prompt_for_eos_params(init,ip,YES,
			"\n\tfor the ambient gas");
	gas_paramsa = prompt_for_eos_params(init,ip,YES,
			"\n\tfor the injected gas");
	gas_paramsc = prompt_for_eos_params(init,ip,YES,
			"\n\tfor the bubbles in the injected gas");

	prompt_for_ambient_state(comp_type(COMPB),gas_paramsb,
                                 " for the ambient gas",front,init);
	prompt_for_ambient_state(comp_type(COMPA),gas_paramsa,
                                 " for the injected gas",front,init);
	prompt_for_ambient_state(comp_type(COMPC),gas_paramsc,
		                 " for the bubbles in the injected gas",front,init);

	/* start make surface of jet */
	h_in[0] = h_i;
	h_in[1] = h_l;
	h_in[2] = h_u;

	h_out[0] = h_i;
	h_out[1] = h_l;
	h_out[2] = h_u;

	//printf("#h %24.16e  %24.16e  %24.16e\n", h_i, h_l, h_u);
	
	jet_paramsa.cen = cen;
	jet_paramsa.r1 = r_i;
	jet_paramsa.r2 = r_n;
	jet_paramsa.h = h_in;
	jet_paramsa.ratio = ratio;
	
	jet_paramsb.cen = cen;
	jet_paramsb.r1 = r_i + w_w;
	jet_paramsb.r2 = r_n;
	jet_paramsb.h = h_out;
	jet_paramsb.ratio = ratio;

	bp.L = gr->GL;         bp.U = gr->GU;
   	
	fuel_params.bp = &bp;
	fuel_params.params_in = &jet_paramsa;
	fuel_params.params_out = &jet_paramsb;
	
	fuel_params.comp_top = exterior_component(front->interf);
	fuel_params.comp_obst = COMPOBST;
	fuel_params.comp_fluid = COMPA;
	fuel_params.comp_gas = COMPB;
	fuel_params.comp_bot = COMPD;
	fuel_params.cen = cen;

	uni_array(&surfs,10,sizeof(SURFACE*));
	uni_array(&jet,10,sizeof(CURVE*));

	is = ic = 0;
	make_surfaces_from_comp(gr, flowin_comp_func_x, (POINTER)&fuel_params, 
			surfs, jet, &is, &ic);

	//make_surfaces_from_files(gr, surfs, jet, &is, &ic); 
	
	printf("End make_jet_surfaces,is = %d,ic = %d\n", is, ic);

	order_interface(intfc);

	//printf("Af init_fuel\n");
	//print_interface(intfc);
	
	//tmp fixed layer
	//move_ref_interface(intfc);

	null_sides_are_consistent();
	check_print_intfc("After init jet", "init_jet", 'f', 
	       intfc, 1, 1, NO);

	init_types_from_comp((POINTER)&fuel_params);
	//do not change this one.
	init_wave_curve_types(intfc, surfs, jet, is, ic);
	init_rect_boundary(front);
	init_bdry_surfs_flag(intfc, (POINTER)&fuel_params, init, ip);

	//print_interface(intfc);

	debug_print("jet","Left init_3comp_jet3d()\n");
}	/* end init_fuel_injection_jet3d */

void	pert_interface_from_file(
		INTERFACE	*intfc)
{
	HYPER_SURF              *hs;
	HYPER_SURF_ELEMENT      *hse;
	POINT                   *p;
	float			**file_r, *pt;
	float			r, theta, z, t_frac, z_frac;
	int			zh, th, th1;

	bi_array(&file_r, 100, 360, FLOAT);
	get_file_r(file_r, 100, 360);

        next_point(intfc,NULL,NULL,NULL);
        while(next_point(intfc,&p,&hse,&hs))
	{
	    pt = Coords(p);

	    theta = 180.0*angle(pt[0], pt[1])/PI;
	    th = (int)theta;
	    th = th >= 360 ? 359 : th;
	    t_frac = theta - th;
	    th1 = th+1 >= 360 ? 0 : th+1;
	    z = (pt[2] + 2.2)*20;
	    zh = (int)z;
	    z_frac = z - zh;

	    r = (1.0-z_frac)*((1.0-t_frac)*file_r[zh][th] + t_frac*file_r[zh][th1])
	    	+z_frac*((1.0-t_frac)*file_r[zh+1][th] + t_frac*file_r[zh+1][th1]);

	    pt[0] = r*cos(theta);
	    pt[1] = r*sin(theta);
	}

	free(file_r);
}

EXPORT	void init_3comp_jet3d_prev(
	INIT_DATA	*init,
	INIT_PHYSICS	*ip)
{
	Front           *front = ip->root->front;
	INTERFACE       *intfc = front->interf;
	RECT_GRID       *gr = front->rect_grid;
	CURVE           **jet;
	SURFACE	        **surfs ;
	JET_PARAMS      jet_paramsa, jet_paramsb;
	BDRY_BOX_PARAMS bp;
	FUEL_PARAMS     fuel_params;
	Gas_param       *gas_paramsa, *gas_paramsb,  *gas_paramsc;
	float           *L = gr->GL;
	float           *U = gr->GU;
	float	        *h = gr->h;
	float	        coords[MAXD];
	char 	        s[1024];
	int             i, dim = gr->dim;
	float           r_i, w_w, h_i, r_n, h_l, h_u;
	float           cen[3], h_in[3], h_out[3], ratio;
	float           tol = MIN_SC_SEP(intfc);
        int             is,ic;

	debug_print("jet","Entered init_fuel_injection_jet3d()\n");

	//#bjet2 TMP fix
	//make_bubbles(NULL, front);

	set_obstacle_comp_type(comp_type(COMPOBST),front);
	exclude_comp(COMPD,intfc);
	set_obstacle_comp_type(comp_type(COMPD),front);

	/* read geometrics of the jet */
	screen("\nFuel injection jet enters from below.\n");
	screen("Please enter the geometric parameters of the fuel jet.\n");
	
	screen("Enter the center values of the floor of fuel container: ");
	(void) Scanf("%f %f %f\n",&cen[0],&cen[1],&cen[2]);
	if ((L[0] > cen[0]) || (U[0] < cen[0]) ||
	    (L[1] > cen[1]) || (U[1] < cen[1]) ||
	    (L[2] > cen[2]) || (U[2] < cen[2]))
	{
	    screen("ERROR in init_fuel_injection_jet3d(), "
	       "invalid center of the floor x0 = %"FFMT", "
	       "y0 = %"FFMT", z0 = %"FFMT"\n",cen[0],cen[1],cen[2]);
	    (void) printf("L[0] = %"FFMT", U[0] = %"FFMT"\n",L[0],U[0]);
	    (void) printf("L[1] = %"FFMT", U[1] = %"FFMT"\n",L[1],U[1]);
	    (void) printf("L[2] = %"FFMT", U[2] = %"FFMT"\n",L[2],U[2]);
	    clean_up(ERROR);
	}

	screen("Enter the inner radius of the fuel container: ");
	(void) Scanf("%f\n",&r_i);
	if ((r_i <= 0.0) || (U[0] < (L[0] + r_i)))
	{
	    screen("ERROR in init_fuel_injection_jet3d(), "
	       "invalid inner radius r_i = %"FFMT"\n",r_i);
	    (void) printf("L[0] = %"FFMT", U[0] = %"FFMT"\n",L[0],U[0]);
	    clean_up(ERROR);
	}

	screen("Enter the width of the fuel container wall: ");
        (void) Scanf("%f\n",&w_w);
	if (w_w <= 0.0)
	{
	    screen("ERROR in init_fuel_injection_jet3d(), "
		       "invalid container width = %"FFMT"\n",w_w);
	    clean_up(ERROR);
	}
   
	screen("Enter the radium of the nozzle: ");
	(void) Scanf("%f\n",&r_n);
	if ((r_n <= 0.0) || (r_i < r_n) || (U[0] < (L[0] + r_n)))
	{
	    screen("ERROR in init_fuel_injection_jet3d(), "
	       "invalid nozzle radium r_n = %"FFMT"\n",r_n);
	    (void) printf("L[0] = %"FFMT", U[0] = %"FFMT"\n",L[0],U[0]);
	    clean_up(ERROR);
	}
   
	h_i = 0.0;
	screen("Enter the height of the interior vertical wall\n"
	   "\tof the fuel chamber (default = %g): ",h_i);
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    sscan_float(s,&h_i);
	    if ((h_i < 0.0) || (U[2] < (L[2] + h_i)))
	    {
	        screen("ERROR in init_fuel_injection_jet3d(), "
		           "invalid interior wall height h_i = %"FFMT"\n",h_i);
	        (void) printf("L[2] = %"FFMT", U[2] = %"FFMT"\n",L[2],U[2]);
	        clean_up(ERROR);
	    }
	}

	h_l = 0.0;
	screen("Enter the height of the curve section of the nozzle\n"
	       "\t(default or zero = no lower section): ");
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    if ((sscan_float(s,&h_l) != 1) || (h_l < 0.0) ||
	        (U[2] < (L[2] + h_i + h_l) ))
	    {
	        screen("ERROR in init_fuel_injection_jet3d(), "
	           "invalid curve section nozzle "
		   "height nh_l = %"FFMT"\n",h_l);
	        (void) printf("L[2] = %"FFMT", U[2] = %"FFMT"\n",L[2],U[2]);
	        clean_up(ERROR);
	    }
	}

	h_u = 0.0;
	screen("Enter the height of the cylinder section of the nozzle\n"
	           "\t(default or zero = no cylinder section): ");
	(void) Gets(s);
	if (s[0] != '\0')
	{
	    if ((sscan_float(s,&h_u) != 1) || (h_u < 0.0) ||
	        (U[2] < (L[2] + h_i + h_l + h_u) ))
	    {
	        screen("ERROR in init_fuel_injection_jet3d(), "
	           "invalid cylinder section nozzle "
		   "height nh_u = %"FFMT"\n",
	           h_u);
	        (void) printf("L[2] = %"FFMT", U[2] = %"FFMT"\n",L[2],U[2]);
	        clean_up(ERROR);
	    }
	}
	
	ratio = 0.5;
	screen("Enter the ratio of the height of the contace section\n"
	       "\t to the top of the nozzle (default = %f): ",ratio);
	(void) Scanf("%f\n",&ratio);
	if ((ratio > 1.0))
	{
	    screen("ERROR in init_fuel_injection_jet3d(), "
	       "invalid ratio of contace in nozzle = %"FFMT"\n",ratio);
	    clean_up(ERROR);
	}

	g_prompt_for_dynamic_bubble_insertion3d(ip->root);

	(void) prompt_for_eos_params(init,ip,YES,"");
	gas_paramsb = prompt_for_eos_params(init,ip,YES,
			"\n\tfor the ambient gas");
	gas_paramsa = prompt_for_eos_params(init,ip,YES,
			"\n\tfor the injected gas");
	gas_paramsc = prompt_for_eos_params(init,ip,YES,
			"\n\tfor the bubbles in the injected gas");

	prompt_for_ambient_state(comp_type(COMPB),gas_paramsb,
                                 " for the ambient gas",front,init);
	prompt_for_ambient_state(comp_type(COMPA),gas_paramsa,
                                 " for the injected gas",front,init);
	prompt_for_ambient_state(comp_type(COMPC),gas_paramsc,
		                 " for the bubbles in the injected gas",front,init);

	/* start make surface of jet */
	h_in[0] = h_i;
	h_in[1] = h_l;
	h_in[2] = h_u;

	h_out[0] = h_i;
	h_out[1] = h_l;
	h_out[2] = h_u;

	//printf("#h %24.16e  %24.16e  %24.16e\n", h_i, h_l, h_u);
	
	jet_paramsa.cen = cen;
	jet_paramsa.r1 = r_i;
	jet_paramsa.r2 = r_n;
	jet_paramsa.h = h_in;
	jet_paramsa.ratio = ratio;
	
	jet_paramsb.cen = cen;
	jet_paramsb.r1 = r_i + w_w;
	jet_paramsb.r2 = r_n;
	jet_paramsb.h = h_out;
	jet_paramsb.ratio = ratio;

	bp.L = gr->GL;         bp.U = gr->GU;
   	
	fuel_params.bp = &bp;
	fuel_params.params_in = &jet_paramsa;
	fuel_params.params_out = &jet_paramsb;
	
	fuel_params.comp_top = exterior_component(front->interf);
	fuel_params.comp_obst = COMPOBST;
	fuel_params.comp_fluid = COMPA;
	fuel_params.comp_gas = COMPB;
	fuel_params.comp_bot = COMPD;
	
	uni_array(&surfs,10,sizeof(SURFACE*));
	uni_array(&jet,10,sizeof(CURVE*));

	is = ic = 0;
	make_surfaces_from_comp(gr, fuel_comp_func, (POINTER)&fuel_params, 
			surfs, jet, &is, &ic);
         
	printf("End make_jet_surfaces,is = %d,ic = %d\n", is, ic);

	order_interface(intfc);

	//printf("Af init_fuel\n");
	//print_interface(intfc);
	null_sides_are_consistent();
	check_print_intfc("After init jet", "init_jet", 'f', 
	       intfc, 1, -1, NO);

	init_types_from_comp((POINTER)&fuel_params);
	//do not change this one.
	init_wave_curve_types(intfc, surfs, jet, is, ic);
	init_rect_boundary(front);
	init_bdry_surfs_flag(intfc, (POINTER)&fuel_params, init, ip);

	//print_interface(intfc);
	//clean_up(0);

	debug_print("jet","Left init_fuel_injection_jet3d()\n");
}	/* end init_fuel_injection_jet3d */

/*
LOCAL   float jet_inner_func(
        POINTER func_params,
        float *coords)
{
JET_PARAMS	*params = (JET_PARAMS *)func_params;
float 		*cen,r0,r1,r2,h1,h2,h3,ratio,arg;

        cen = params->cen;
        r1 = params->r1;
        r2 = params->r2;
        ratio = params->ratio;

        h1 = cen[2] - params->h[0];
        h2 = h1 - params->h[1];
        h3 = h2 - params->h[2]*(1.0-ratio);
	//h3 = cen[2] - h[0] - h[1] - h[2] + h[2]*ratio, z coord of contact

        arg = 1.0;
        if (coords[2] >= h1)
            arg = sqr(coords[0]-cen[0])+sqr(coords[1]-cen[1])-sqr(r1);
        else if (coords[2] >= h3)
        {
            r0 = sqrt(fabs(sqr(r1-r2)-sqr(coords[2]-h1))) + r2;
            if (coords[2] >= h2)
                arg = sqrt(sqr(coords[0]-cen[0])+sqr(coords[1]-cen[1])) - r0;
            else
                arg = sqr(coords[0]-cen[0])+sqr(coords[1]-cen[1])-sqr(r2);
        }
        
        return arg;
}       // end jet_inner_func 

LOCAL   float jet_outer_func(
        POINTER func_params,
        float *coords)
{
JET_PARAMS	*params = (JET_PARAMS *)func_params;
float		dist, h1, h2, *cen;

	cen = params->cen;
	h1 = cen[2] - params->h[0] - params->h[1] - params->h[2];
	h2 = h1 + params->h[2]*params->ratio;   //z coord of contact

	if(coords[2] >= h1)
	{
	    dist = sqrt(sqr(coords[0]-cen[0])+sqr(coords[1]-cen[1]));

	    if(dist > params->r1)
	        return 1.0;
	    if(coords[2] > h2)
	        return -1.0;
	    if(dist < params->r2)
	        return 1.0;
	    return -1.0;
	}

        return 1.0;
}      // end jet_outer_func
*/
// straight jet
LOCAL   float jet_inner_func(
        POINTER func_params,
        float *coords)
{
JET_PARAMS	*params = (JET_PARAMS *)func_params;
float 		*cen,r2,h1,dist;

        cen = params->cen;
        r2 = params->r2;

	h1 = cen[2] - params->h[2]*(1.0-params->ratio);
	dist = sqrt(sqr(coords[0]-cen[0])+sqr(coords[1]-cen[1]));
        
        if (dist < r2 && coords[2] > h1)
	    return -1.0;
        
        return 1.0;
}       // end jet_inner_func 

LOCAL   float jet_outer_func(
        POINTER func_params,
        float *coords)
{
JET_PARAMS	*params = (JET_PARAMS *)func_params;
float		dist, h1, h2, r1, r2, *cen;

	cen = params->cen;
	r1 = params->r1;
	r2 = params->r2;
	h1 = cen[2] - params->h[2]*(1.0-params->ratio);
	h2 = cen[2] - params->h[2];

	dist = sqrt(sqr(coords[0]-cen[0])+sqr(coords[1]-cen[1]));
	
	//printf("dist=%24.16e  r2=%24.16e h2=%24.16e h1=%24.16e\n", dist, r2, h2, h1);
	
	if((dist < r2 && coords[2] > h1) || 
	   (dist >= r2 && dist < r1 && coords[2] > h2))
	    return -1.0;
	
	return 1.0;

}      // end jet_outer_func

/*
LOCAL   int   fuel_comp_func(
        POINTER func_params,
        float *coords)
{
FUEL_PARAMS      *params = (FUEL_PARAMS *) func_params;
JET_PARAMS       *params_in = params->params_in;
JET_PARAMS       *params_out = params->params_out;
BDRY_BOX_PARAMS  *bp = params->bp;
float		 fin, fout;

	if(coords[2] > bp->U[2])
	    return  params->comp_top;
	//if(coords[2] < bp->L[2])
	//    return  params->comp_bot;
	
	fin = jet_inner_func((POINTER)params_in, coords);
	if(fin < 0)
	    return params->comp_fluid;
	fout = jet_outer_func((POINTER)params_out, coords);
	if(fout < 0)
	    return params->comp_obst;

	return params->comp_gas;
}
*/

LOCAL   int   fuel_comp_func(
        POINTER func_params,
        float *coords)
{
FUEL_PARAMS      *params = (FUEL_PARAMS *) func_params;
JET_PARAMS       *params_in = params->params_in;
JET_PARAMS       *params_out = params->params_out;
BDRY_BOX_PARAMS  *bp = params->bp;
int		 ind, indth;
float		 h, r, len, theta, z, L, cl, dl, z1, coef, r_ellip, phi;


	r = sqrt(sqr(coords[0]/0.3) + sqr(coords[1]/0.84));
	
	//r = sqrt(sqr(coords[0]/0.5012) + sqr(coords[1]/0.5012) + sqr(coords[2]/0.5012));
	if(r < 1.0)
	    return  params->comp_fluid;
	return  params->comp_gas;

	r = sqrt(sqr(coords[0]) + sqr(coords[1]));

	/*
	z = coords[2] + 1.6;
	theta = angle(coords[0], coords[1]);
	
	//pert from a jet simulation
	//0.8/161 will give error so that scatter front fails.
	
	h = 0.8/161.0;
	ind = (int)fabs(z/h);
	ind = (ind/80)%2 == 0 ? ind%80 : 79-ind%80;

	//r_ellip = sqrt(sqr(0.2835*cos(theta)) + sqr(0.7155*sin(theta)));
	//r_ellip = 0.2835;
	
	//ellip radius in direction.
	phi = angle(coords[0]/0.2835, coords[1]/0.7155);
	r_ellip = sqrt(sqr(0.2835*cos(phi)) + sqr(0.7155*sin(phi)));
	
	//indth = (int)(theta/PI*180-0.001);
	indth = (int)(phi/PI*180-0.001);
	len = (params->file_r[ind][indth]/0.0091-1.0)*0.4*r_ellip + r_ellip;
	
	z1 = fabs(coords[2]);
	if(z1>1.43)
	    if(z1<1.55)
		coef = (z1-1.43)/(1.55-1.43);
	    else
		coef = 1.0;
	else
	    coef = 0.0;
	
	len = len*(1.0-coef) + r_ellip*coef;

	if(debugging("scatcrx"))
	{
	    print_general_vector("crds ", coords, 3, "  ");
	    printf(" len %15.8e %d %d\n", len, ind, indth);
	}
	*/

	len = 0.5021 + 0.07*exp(-32.0*sqr(coords[2]));

	if(len < r)
	//if(r > 0.2835)
	    return params->comp_gas;

	return params->comp_fluid;
}

/*
LOCAL   int   fuel_comp_func(
        POINTER func_params,
        float *coords)
{
FUEL_PARAMS      *params = (FUEL_PARAMS *) func_params;
JET_PARAMS       *params_in = params->params_in;
JET_PARAMS       *params_out = params->params_out;
BDRY_BOX_PARAMS  *bp = params->bp;
float		 h, r, len, theta, z, L, cl, dl;

	//r = sqr(coords[0])/sqr(0.483) + sqr(coords[1])/sqr(0.187);
	r = sqr(coords[0])/sqr(0.83) + sqr(coords[1])/sqr(0.3);

	if(r < 1.0)
	    return params->comp_fluid;
	else
	    return params->comp_gas;

}
*/

LOCAL   int   flowin_comp_func_prev(
        POINTER func_params,
        float *coords)
{
FUEL_PARAMS      *params = (FUEL_PARAMS *) func_params;
JET_PARAMS       *params_in = params->params_in;
JET_PARAMS       *params_out = params->params_out;
BDRY_BOX_PARAMS  *bp = params->bp;
float		 r2, *cen;
        
	cen = params_in->cen;
	r2 = params_in->r2;

	/* 100325 change YYU */
//	if(coords[2] > bp->U[2])
//	    return  params->comp_top;
	
	if(sqr(coords[0]-cen[0]) + sqr(coords[1]-cen[1]) < sqr(r2) &&
	       coords[2] > bp->U[2] - 0.012*2)
	    return params->comp_fluid;

        /* Crosslfow */
        if(coords[0] < bp->L[0])
            return params->comp_bot;

	return params->comp_gas;
}


LOCAL   int   flowin_comp_func(
        POINTER func_params,
        float *coords)
{
FUEL_PARAMS      *params = (FUEL_PARAMS *) func_params;
JET_PARAMS       *params_in = params->params_in;
JET_PARAMS       *params_out = params->params_out;
BDRY_BOX_PARAMS  *bp = params->bp;
float		 r2, *cen, r1, wall;
        
	cen = params_in->cen;
	r2 = params_in->r2;
	
	if(coords[2] > bp->U[2])
	    return  params->comp_top;
	
	//xiaoxue
//	if(coords[0] < bp->L[0] + r2/3.77
//	    return params->comp_bot;
	if(coords[0] < bp->L[0]) //TK
	    return params->comp_bot;
/*
	if(coords[2] < bp->U[2]-0.15)
	    return params->comp_gas;
*/
      if(coords[1] < 0.0 || coords[1] > 3.75) return params->comp_obst;
	
        if(coords[2] > bp->U[2]-0.1 && coords[2] <= bp->U[2])
	{
	    if(sqr(coords[0]-cen[0]) + sqr(coords[1]-cen[1]) < sqr(r2))
	    	return params->comp_gas;
	}

	/* Scranjet configuration */

	if(coords[0] > bp->L[0] && coords[0] <= -7)
	    wall = bp->L[2]+0.1;
	if(coords[0] > -7 && coords[0] <= -7+0.8/tan(3.1415926*10/180))
	    wall = bp->L[2]+0.1+(coords[0]+7)*tan(3.1415926*10/180);
	if(coords[0] > -7+0.8/tan(3.1415926*10/180))
	    wall = bp->L[2]+0.1+0.8;
//	printf("%lf, %lf, %lf, r1=%lf\n",coords[0],coords[1],coords[2],wall);
	
	if(coords[2] < bp->U[2]-0.1 && coords[2] > wall)
	    return params->comp_gas;

/*
	if(coords[2] > bp->U[2]-1.0 && coords[2] < bp->U[2]-1.0+0.095)
	    r1 = r2;
	if(coords[2] >= bp->U[2]-1.0+0.095 && coords[2] < bp->U[2]-1.0+0.095+0.317)
	    r1 = r2+0.556-sqrt(sqr(0.556)-sqr(coords[2]-(bp->U[2]-1.0+0.095)));
	if(coords[2] >= bp->U[2]-1.0+0.095+0.317 && coords[2] < bp->U[2]-1.0+0.095+0.317+0.061)
	    r1 = r2+0.556-sqrt(sqr(0.556)-sqr(0.317))+(coords[2]-(bp->U[2]-1.0+0.095+0.317))*tan(3.14159*67.5/180);
	if(coords[2] >= bp->U[2]-1.0+0.095+0.317+0.061 && coords[2] <= bp->U[2])
	    r1 = r2+0.556-sqrt(sqr(0.556)-sqr(0.317))+0.061*tan(3.14159*67.5/180);
*/
//	printf("%lf, %lf, %lf, r1=%lf\n",coords[0],coords[1],coords[2],r1);
	
//	if(sqr(coords[0]-cen[0]) + sqr(coords[1]-cen[1]) < sqr(r2))
//	    return params->comp_gas;

        return params->comp_obst;
}


LOCAL   int   flowin_comp_func_x(
        POINTER func_params,
        float *coords)
{
FUEL_PARAMS      *params = (FUEL_PARAMS *) func_params;
    return 3;

/*JET_PARAMS       *params_in = params->params_in;
JET_PARAMS       *params_out = params->params_out;
BDRY_BOX_PARAMS  *bp = params->bp;
float		 r2, *cen, r1, wall;
        
	cen = params_in->cen;
	r2 = params_in->r2;
	
	if(coords[2] > bp->U[2])
	    return  params->comp_top;
	
	if(coords[0] < bp->L[0]) //TK
	    return params->comp_bot;
      
        if(coords[1] < 0.0 || coords[1] > 3.75) return params->comp_obst;
	
        if(coords[2] > bp->U[2]-0.1 && coords[2] <= bp->U[2])
	{
	    if(sqr(coords[0]-cen[0]) + sqr(coords[1]-cen[1]) < sqr(r2))
	    	return params->comp_gas;
	}


	if(coords[0] > bp->L[0] && coords[0] <= -7)
	    wall = bp->L[2]+0.1;
	if(coords[0] > -7 && coords[0] <= -7+0.8/tan(3.1415926*10/180))
	    wall = bp->L[2]+0.1+(coords[0]+7)*tan(3.1415926*10/180);
	if(coords[0] > -7+0.8/tan(3.1415926*10/180))
	    wall = bp->L[2]+0.1+0.8;
	
	if(coords[2] < bp->U[2]-0.1 && coords[2] > wall)
	    return params->comp_gas;

        return params->comp_obst;*/
}

LOCAL   int   cyl_pert_comp_func(
        POINTER  func_params,
        float 	 *coords)
{
	FUEL_PARAMS      *params = (FUEL_PARAMS *) func_params;
	BDRY_BOX_PARAMS  *bp = params->bp;
	float		 lambda, amp, r, z, *cen, kx, theta;

	cen = params->cen;
	lambda = (bp->U[2] - bp->L[2])/1.0;
	amp = 0.05111*lambda;
	kx = 2*PI/lambda;
	r = sqrt(sqr(coords[0]-cen[0]) + sqr(coords[1]-cen[1]));
	z = coords[2];
	theta = angle(coords[0], coords[1]);
	
	//if(r < params->r + amp*sin(kx*z)*sin(8.0*theta))
	if(r < params->r + amp*sin(kx*z))
	    return 2;
	return 3;
}

LOCAL   int   plane_pert_comp_func(
        POINTER  func_params,
        float 	 *coords)
{
	FUEL_PARAMS      *params = (FUEL_PARAMS *) func_params;
	BDRY_BOX_PARAMS  *bp = params->bp;
	float		 lambda, lambda1, amp, r, z, *cen, kx, ky, theta;

	cen = params->cen;
	lambda = (bp->U[2] - bp->L[2])/1.0;
	lambda1 = (bp->U[0] - bp->L[0])/1.0;
	
	amp = 0.05111*lambda;
	kx = 2*PI/lambda;
	ky = 2*PI/lambda1;
	
	if(coords[1] < params->r + amp*sin(kx*coords[2])*sin(ky*(coords[0]-lambda1/4.0)))
	    return 2;
	return 3;
}


LOCAL   void   get_file_r(
	float	**file_r, 
	int	kz, 
	int	kt)
{
	int		i, j;
	FILE		*fp;
	char		fname[200];
		
	//sprintf(fname, "%ssurf_polar", get_directory());
	fp = fopen(fname, "r");
	for(i=0; i<kz; i++)
	    for(j=0; j<kt; j++)
	    {
		fscanf(fp, "%lf", &file_r[i][j]);
	    }
	fclose(fp);

}
