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
*				example.c:
*
*		User initialization example for Front Package:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*/

#include <FronTier.h>

	/*  Function Declarations */
static void test_propagate(Front*);
static float test_ellipoid_func(POINTER,float*);
static int test_circular_rotation_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,HYPER_SURF*,float*);


char *in_name,*restart_state_name,*restart_name,*out_name;
bool RestartRun;
int RestartStep;
bool binary = YES;

/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/
// Slotted sphere in 3D evolves under circular rotation around z-axis
//
//
typedef struct {
        float cen[3];
	float rad[3];
	float h;
	float w;
} TEST_ELLIPSOID_PARAMS;

typedef struct
{       
       int dim;
       float x_0,y_0;
       float omega_0;
       float grad;
} TEST_CIRCULAR_ROTATION_PARAMS;

/********************************************************************
 *	Velocity function parameters for the front	 	    *
 ********************************************************************/


int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	TEST_ELLIPSOID_PARAMS s_params;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	TEST_CIRCULAR_ROTATION_PARAMS cr_params; /* velocity function parameters */
	Locstate  sl;

	f_basic.dim = 3;	
	FrontInitStandardIO(argc,argv,&f_basic);

	/* Initialize basic computational data */

	f_basic.L[0] = 0.0;	f_basic.L[1] = 0.0; 	f_basic.L[2] = 0.0;
	f_basic.U[0] = 1.0;	f_basic.U[1] = 1.0; 	f_basic.U[2] = 1.0;
	f_basic.gmax[0] = 100;	f_basic.gmax[1] = 100; f_basic.gmax[2] = 100;
	f_basic.boundary[0][0] = f_basic.boundary[0][1] = PERIODIC_BOUNDARY;
	f_basic.boundary[1][0] = f_basic.boundary[1][1] = PERIODIC_BOUNDARY;
	f_basic.boundary[2][0] = f_basic.boundary[2][1] = PERIODIC_BOUNDARY;
	f_basic.size_of_intfc_state = 0;

        in_name                 = f_basic.in_name;
        restart_state_name      = f_basic.restart_state_name;
        out_name                = f_basic.out_name;
        restart_name            = f_basic.restart_name;
        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;

        sprintf(restart_name,"%s.ts%s",restart_name,right_flush(RestartStep,7));
#if defined(__MPI__)
        sprintf(restart_name,"%s-nd%s",restart_name,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */

	FrontStartUp(&front,&f_basic);

	if (!RestartRun)
	{
	    /* Initialize interface through level function */
	    s_params.cen[0] = 0.5;
	    s_params.cen[1] = 0.5;
	    s_params.cen[2] = 0.5;
	    s_params.rad[0] = 0.3;
	    s_params.rad[1] = 0.3;
	    s_params.rad[2] = 0.3;
	    s_params.h = 0.4;
	    s_params.w = 0.025;

	    level_func_pack.neg_component = 1;
	    level_func_pack.pos_component = 2;
	    level_func_pack.func_params = (POINTER)&s_params;
	    level_func_pack.func = test_ellipoid_func;

	    FrontInitIntfc(&front,&level_func_pack);
	}

	/* Initialize velocity field function */

	cr_params.dim = 3;
	cr_params.x_0 = 0.5;
	cr_params.y_0 = 0.5;
	cr_params.omega_0 = 1.0;
	cr_params.grad = 0.0;

	velo_func_pack.func_params = (POINTER)&cr_params;
	velo_func_pack.func = test_circular_rotation_vel;

	FrontInitVelo(&front,&velo_func_pack);

	/* For geometry-dependent velocity, use first 
	* order point propagation function, higher order
	* propagation requires surface propagate, currently
	* in writing, not yet in use. The following override
	* the assigned fourth_order_point_propagate.
	*/

	/* Propagate the front */

	test_propagate(&front);

	clean_up(0);
	return 0;
}

static  void test_propagate(
        Front *front)
{
        int ip,im,status,count;
        Front *newfront;
        float dt,dt_frac,CFL;
        bool is_print_time, is_movie_time, time_limit_reached;
        char s[10];
        float fcrds[MAXD];
        int  dim = front->rect_grid->dim;

	front->max_time = 25; 
	front->max_step = 400;
	front->print_time_interval = 2;
	front->movie_frame_interval = 0.1;

        CFL = Time_step_factor(front) = 0.5;

        printf("dim = %d\n", dim);
	printf("CFL = %f\n",CFL);
	printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
		Frequency_of_redistribution(front,GENERAL_WAVE));

	if (!RestartRun)
	{
            redistribute(front,YES,NO);

            front->time = 0.0;
            front->dt = 0.0;
	    front->step = 0;

	    // Always output the initial interface.
	    FrontPrintOut(front,out_name);
            FrontMovieFrame(front,out_name,binary);
            ip = im = 1;

	    // This is a virtual propagation to get maximum front 
	    // speed to determine the first time step.

            status = FrontAdvance(front->dt,&dt_frac,front,&newfront,
                                (POINTER)NULL);
            front->dt = CFL*FrontHypTimeStep(front); 
	}
	else
	{
	    ip = (int)(front->time/front->print_time_interval + 1.0);
            im = (int)(front->time/front->movie_frame_interval + 1.0);
	}

	front->dt = FrontOutputTimeControl(front,
			&is_movie_time,&is_print_time,
			&time_limit_reached,&im,&ip);

        for (;;)
        {
	    /* Propagating interface for time step dt */

            status = FrontAdvance(front->dt,&dt_frac,front,&newfront,
                                (POINTER)NULL);
            assign_interface_and_free_front(front,newfront);

            ++front->step;
            front->time += front->dt;

	    //Next time step determined by maximum speed of previous
	    //step, assuming the propagation is hyperbolic and
	    //is not dependent on second order derivatives of
	    //the interface such as curvature, and etc.

            front->dt = CFL*FrontHypTimeStep(front); 

            printf("\ntime = %f   step = %5d   next dt = %f\n",
                        front->time,front->step,front->dt);
            fflush(stdout);

            if (is_print_time || time_limit_reached)
                print_front_output(front,out_name);
            if (is_movie_time || time_limit_reached)
                FrontMovieFrame(front,out_name,binary);

            if (time_limit_reached)
                    break;

	    /* Output section, next dt may be modified */

	    front->dt = FrontOutputTimeControl(front,
			&is_movie_time,
			&is_print_time,
			&time_limit_reached,
			&im,&ip);
        }
        (void) delete_interface(front->interf);
}       /* end test_propagate */

/********************************************************************
 *	Sample (dummbell 3D) level function for the initial interface    *
 ********************************************************************/

static float test_ellipoid_func(
        POINTER func_params,
        float *coords)
{
        TEST_ELLIPSOID_PARAMS *e_params = (TEST_ELLIPSOID_PARAMS*)func_params;
        const float *cen,*rad;         
	float arg;          
	float x0,y0,z0,h,w;
	float dist1,dist2,dist3;

	x0 = e_params->cen[0]; y0 = e_params->cen[1]; z0 = e_params->cen[2];
	rad = e_params->rad;
	h = e_params->h;
	w = e_params->w;

        dist1 = 1.0 -
                sqr(coords[0] - x0)/sqr(rad[0]) -
                sqr(coords[1] - y0)/sqr(rad[1]) -
                sqr(coords[2] - z0)/sqr(rad[2]);

	dist2 = fabs(coords[0]-x0)-w;
	dist3 = coords[1] - y0 -(rad[1]-h);
	if (dist1 < 0)
	    arg = dist1;
	else
	{
	    if (dist2 <= 0)
	    {
	    	if(dist3>=0)
		    arg = -1*min(fabs(dist2),fabs(dist3));
		else
		    arg = min(fabs(dist1),fabs(dist3)); 
	    }
	    else
	    {
	    	if(dist3>=0)
		    arg = min(fabs(dist1), fabs(dist2));
		else
		    arg = dist1;
	    }
	}
        return -arg;
}       /* end test_ellipoid_func */


static int test_circular_rotation_vel(
	POINTER params,
	Front *front,
	POINT *p,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	float *vel)
{
	TEST_CIRCULAR_ROTATION_PARAMS *c_params;
	float x0,y0;
	float omega_0,domega_dr;
	float rad,V;
	float xcomp,ycomp;
	float *coords = Coords(p);

	c_params = (TEST_CIRCULAR_ROTATION_PARAMS*)params;

	x0 = c_params->x_0;
	y0 = c_params->y_0;
	omega_0 = c_params->omega_0;
	domega_dr = c_params->grad;

	rad = sqrt(sqr(coords[0] - x0) + sqr(coords[1] - y0));
	if (rad == 0.0)
	{
	    vel[0] = vel[1] = vel[2] = 0.0;
	    return 1;
	}
	xcomp = fabs(coords[1]-y0)/rad;
	ycomp = fabs(coords[0]-x0)/rad;
	V = rad*(omega_0 + domega_dr*rad);

	vel[2] = 0.0;
	if (coords[0]-x0 >= 0.0 && coords[1]-y0 >= 0.0) /*1st quadrant*/ 
	{
	    vel[0] = -V*xcomp;
	    vel[1] =  V*ycomp;
	}
	else if (coords[0]-x0 <= 0.0 && coords[1]-y0 >= 0.0) /*2nd quadrant*/ 
	{
	    vel[0] = -V*xcomp;
	    vel[1] = -V*ycomp;
	}
	else if (coords[0]-x0 <= 0.0 && coords[1]-y0 <= 0.0) /*3rd quadrant*/ 
	{
	    vel[0] =  V*xcomp;
	    vel[1] = -V*ycomp;
	}
	else if (coords[0]-x0 >= 0.0 && coords[1]-y0 <= 0.0) /*4th quadrant*/ 
	{
	    vel[0] =  V*xcomp;
	    vel[1] =  V*ycomp;
	}
	vel[2] = 0;
}	/* end shear_motion_vel */
