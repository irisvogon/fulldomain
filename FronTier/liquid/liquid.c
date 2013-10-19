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
*				test_incompressible.c
* This program is modified from example0.c for solving incompressible flow.
* 
* The solver is define in lcartsn.h/c.
*
*
*				example0.c:
*
*		User initialization example for Front Package:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*	
*	This example shows a circle in a double vortex field. It demonstrates
*	the resolution of the front tracking method.
*
*/

#include "liquid.h"

	/*  Function Declarations */
static void init_io( int,char**);
static void fluid_driver(Front*,CARTESIAN&);
static double level_circle_func(POINTER,double*);
static double sine_curve(POINTER,double*);
static int test_vortex_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
	                       HYPER_SURF*,double*);
static int cartesian_vel(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
                           HYPER_SURF*,double*);

char *in_name,*restart_state_name,*restart_name,*out_name;
bool RestartRun;
bool ReadFromInput;
int RestartStep;
bool binary = NO;


/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/

typedef struct {
	        /* equation for line is x^2/a^2 + y^2/b^2 = 1 */
        double x0;
        double y0;        
        double R; 
} CIRCLE_PARAMS;


/********************************************************************
 *	Velocity function parameters for the front	 	    *
 ********************************************************************/

typedef struct {
	double i1,i2;
        double cen1[2],cen2[2];
} DOUBLE_VORTEX_PARAMS;

	/* It is very strange that after using pbs instead of qsub -I 
	 * to login a node, I can make and after I make the file using 
	 * starzero, test_incompressible doesn't converge. 
	 */

int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	CIRCLE_PARAMS circle_params;	/* level function parameters */

	CARTESIAN	cartesian(front);


	/* Initialize basic computational data */

	f_basic.dim = 2;	// default dimension	
	FrontInitStandardIO(argc,argv,&f_basic);

	f_basic.L[0] = 0.0;	f_basic.L[1] = 0.0;
	f_basic.U[0] = 1.0;	f_basic.U[1] = 4.0;
	
	f_basic.gmax[0] = 256;	f_basic.gmax[1] = 1024;
	//f_basic.gmax[0] = 25;	f_basic.gmax[1] = 100;
	f_basic.boundary[0][0] = f_basic.boundary[0][1] = PERIODIC_BOUNDARY;
	//f_basic.boundary[1][0] = f_basic.boundary[1][1] = PERIODIC_BOUNDARY;
	f_basic.boundary[1][0] = f_basic.boundary[1][1] = NEUMANN_BOUNDARY;
	f_basic.size_of_intfc_state = 0;

        in_name                 = f_basic.in_name;
        restart_state_name      = f_basic.restart_state_name;
        out_name                = f_basic.out_name;
        restart_name            = f_basic.restart_name;
        RestartRun              = f_basic.RestartRun;
        RestartStep             = f_basic.RestartStep;

        sprintf(restart_state_name,"%s-state.ts%s",restart_name,
                        right_flush(RestartStep,7));
        sprintf(restart_name,"%s.ts%s",restart_name,right_flush(RestartStep,7));
        sprintf(restart_name,"%s-nd%s",restart_name,right_flush(pp_mynode(),4));
        sprintf(restart_state_name,"%s-nd%s",restart_state_name,
                        right_flush(pp_mynode(),4));

	FrontStartUp(&front,&f_basic);
        PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

	/* Initialize interface through level function */

	if (!RestartRun)
	{
	    circle_params.x0 = 0.5;
	    circle_params.y0 = 3.0;
	    circle_params.R = 0.3;

	    level_func_pack.neg_component = 0;
	    level_func_pack.pos_component = 1;
	    level_func_pack.func_params = (POINTER)&circle_params;
	    //level_func_pack.func = sine_curve;
	    level_func_pack.func = level_circle_func;

	    FrontInitIntfc(&front,&level_func_pack);
	}

	/* Initialize velocity field function */

	velo_func_pack.func_params = (POINTER)&cartesian;
	velo_func_pack.func = cartesian_vel;
	velo_func_pack.point_propagate = first_order_point_propagate;

	cartesian.initRectGrid(&front);
	cartesian.initMesh();

	FrontInitVelo(&front,&velo_func_pack);

	/* Propagate the front */

	fluid_driver(&front, cartesian);

	PetscFinalize();
	clean_up(0);
}

static  void fluid_driver(
        Front *front,
	CARTESIAN &cartesian)
{
        int ip,im;
        Front *newfront;
        double dt_frac,CFL,tmp;
	bool is_print_time,is_movie_time,time_limit_reached;
	int i,dim = front->rect_grid->dim;

	Curve_redistribution_function(front) = expansion_redistribute;

	if (ReadFromInput)
        {
	    /*TODO*/
	    ;
	}
	else
	{
	    front->max_time = 0.4;
	    front->max_step = 100000;
	    front->print_time_interval  = 0.1;
	    front->movie_frame_interval = 0.01;
            CFL = Time_step_factor(front) = 0.1;
	    Frequency_of_redistribution(front,GENERAL_WAVE) = 5;
	}

	if (RestartRun)
	{
	    FrontScatIntfcBuffer(front);
	}
	else
	{
	    FrontRedistMesh(front);
            front->time = 0.0;
            front->dt = 0.0;
            front->step = 0; 
	}

	if (RestartRun)
	    cartesian.readMeshStates(restart_state_name);
	else
	    cartesian.setInitialCondition();
	cartesian.interpolateEdgeVelocity();

        if (!RestartRun)
        {
            // Front standard output
            FrontPrintOut(front,out_name);
            FrontMovieFrame(front,out_name,binary);
	    cartesian.printMeshStates(out_name);

            // Problem specific output
	    // TODO

            FrontAdvance(front->dt,&dt_frac,front,&newfront,
                                (POINTER)NULL);
            ip = im = 1;
            front->dt = CFL*FrontHypTimeStep(front);
	    for (i = 0; i < dim; ++i)
	    	front->dt = min(front->dt,CFL*cartesian.top_h[i]);
        }
        else
        {
            ip = (int)(front->time/front->print_time_interval + 1.0);
            im = (int)(front->time/front->movie_frame_interval + 1.0);
        }

        front->dt = FrontOutputTimeControl(front,
                        &is_movie_time,&is_print_time,
                        &time_limit_reached,&im,&ip);

	printf("CFL = %f\n",CFL);
	printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
			Frequency_of_redistribution(front,GENERAL_WAVE));

        for (;;)
        {
            /* Propagating interface for time step dt */

            FrontAdvance(front->dt,&dt_frac,front,&newfront,
                                (POINTER)NULL);

	    cartesian.solve(front->dt);
            assign_interface_and_free_front(front,newfront);

            ++front->step;
            front->time += front->dt;
				
            //Next time step determined by maximum speed of previous
            //step, assuming the propagation is hyperbolic and
            //is not dependent on second order derivatives of
            //the interface such as curvature, and etc.

            front->dt = CFL*FrontHypTimeStep(front);
	    for (i = 0; i < dim; ++i)
	    	front->dt = min(front->dt,CFL*cartesian.top_h[i]);
	
            /* Output section */

            printf("\ntime = %f   step = %5d   next dt = %f\n",
                        front->time,front->step,front->dt);
            fflush(stdout);

            if (is_print_time || time_limit_reached)
	    {
            	FrontPrintOut(front,out_name);
		cartesian.printMeshStates(out_name);
	    }
            if (is_movie_time || time_limit_reached)
            	FrontMovieFrame(front,out_name,binary);

            if (time_limit_reached)
                    break;

            front->dt = FrontOutputTimeControl(front,
                        &is_movie_time,
                        &is_print_time,
                        &time_limit_reached,
                        &im,&ip);
        }
}       /* end fluid_driver */

/********************************************************************
 *	Sample (circle) level function for the initial interface    *
 ********************************************************************/

static double level_circle_func(
        POINTER func_params,
        double *coords)
{
	CIRCLE_PARAMS *circle_params = (CIRCLE_PARAMS*)func_params;
	double x0,y0,R,dist;

	x0 = circle_params->x0;
	y0 = circle_params->y0;
	R  = circle_params->R;

	dist = sqrt(sqr(coords[0] - x0) + sqr(coords[1] - y0)) - R;
	return dist;
}	/* end level_circle_func */

LOCAL double sine_curve(
        POINTER func_params,
        double *coords)
{
        double dist, theta;
	double y_intfc = 3 + 0.1*sin(2.0*PI*(coords[0]+0.25));

        dist =   coords[1] - y_intfc;
        return dist;
}       /* end sine_curve */


static int cartesian_vel(
	POINTER params,
	Front *front,
	POINT *p,
	HYPER_SURF_ELEMENT *hse,
	HYPER_SURF *hs,
	double *vel)
{
	double *coords = Coords(p);
	((CARTESIAN_EB*)params)->getVelocity(coords, vel);
}
