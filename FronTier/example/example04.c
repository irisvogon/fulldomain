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
*				example4.c:
*
*		User initialization example for Front Package:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	This example shows the test of Dirichlet boudary.
*
*/

#include <FronTier.h>

	/*  Function Declarations */
static void test_propagate(Front*);
static float multi_circle_func(POINTER,float*);
static int trans_vel_func(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
	                       HYPER_SURF*,float*);

char *in_name,*restart_state_name,*restart_name,*out_name;
bool RestartRun;
int RestartStep;
bool binary = YES;
/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/

typedef struct {
        int num_cir;
        float **cen;
        float *rad;
} TMC_PARAMS;


/********************************************************************
 *	Velocity function parameters for the front	 	    *
 ********************************************************************/

struct _TRANSV_PARAMS
{
        int dim;
	float vx,vy;
};
typedef struct _TRANSV_PARAMS TRANSV_PARAMS;

int main(int argc, char **argv)
{
	static Front front;
	static RECT_GRID comp_grid;
	static F_BASIC_DATA f_basic;
	static LEVEL_FUNC_PACK level_func_pack;
	static VELO_FUNC_PACK velo_func_pack;
	TRANSV_PARAMS trans_params; /* velocity function parameters */
	TMC_PARAMS mc_params;
	Locstate  sl;

	f_basic.dim = 2;
	FrontInitStandardIO(argc,argv,&f_basic);

	/* Initialize basic computational data */

	f_basic.L[0] = 0.0;	f_basic.L[1] = 0.0;
	f_basic.U[0] = 1.0;	f_basic.U[1] = 1.0;
	f_basic.gmax[0] = 100;	f_basic.gmax[1] = 100;
	f_basic.boundary[0][0] = f_basic.boundary[0][1] = DIRICHLET_BOUNDARY;
	f_basic.boundary[1][0] = f_basic.boundary[1][1] = DIRICHLET_BOUNDARY;
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

	    mc_params.num_cir = 3;
            vector(&mc_params.rad,mc_params.num_cir,FLOAT);
            matrix(&mc_params.cen,mc_params.num_cir,2,FLOAT);
	    mc_params.cen[0][0] = 0.3;
	    mc_params.cen[0][1] = 0.3;
	    mc_params.cen[1][0] = 0.7;
	    mc_params.cen[1][1] = 0.3;
	    mc_params.cen[2][0] = 0.5;
	    mc_params.cen[2][1] = 0.7;
	    mc_params.rad[0] = 0.1;
	    mc_params.rad[1] = 0.1;
	    mc_params.rad[2] = 0.1;

	    level_func_pack.neg_component = 1;
	    level_func_pack.pos_component = 2;
	    level_func_pack.func_params = (POINTER)&mc_params;
	    level_func_pack.func = multi_circle_func;

	    FrontInitIntfc(&front,&level_func_pack);
	}

	/* Initialize velocity field function */

	trans_params.dim = 2;
	trans_params.vx = 0.4;
	trans_params.vy = 0.23463;

	velo_func_pack.func_params = (POINTER)&trans_params;
	velo_func_pack.func = trans_vel_func;
	velo_func_pack.point_propagate = fourth_order_point_propagate;

	FrontInitVelo(&front,&velo_func_pack);

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

	front->max_time = 2.1;
	front->max_step = 10000;
	front->print_time_interval = 1.0;
	front->movie_frame_interval = 0.02;

        CFL = Time_step_factor(front);

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

	    /* Output section */

            printf("\ntime = %f   step = %5d   next dt = %f\n",
                        front->time,front->step,front->dt);
            fflush(stdout);

            if (is_print_time || time_limit_reached)
                print_front_output(front,out_name);
            if (is_movie_time || time_limit_reached)
                show_front_output(front,out_name,binary);

            if (time_limit_reached)
                    break;

	    front->dt = FrontOutputTimeControl(front,
			&is_movie_time,
			&is_print_time,
			&time_limit_reached,
			&im,&ip);
        }
        (void) delete_interface(front->interf);
}       /* end test_propagate */

LOCAL int trans_vel_func(
        POINTER params,
        Front *front,
        POINT *p,
        HYPER_SURF_ELEMENT *hse,
        HYPER_SURF *hs,
        float *vel)
{
        TRANSV_PARAMS *transv_params;
                                                                                
        transv_params = (TRANSV_PARAMS*)params;
                                                                                
	vel[0] = transv_params->vx;
	vel[1] = transv_params->vy;
}       /* end transal_vel_func */

LOCAL float multi_circle_func(
        POINTER func_params,
        float *coords)
{
        TMC_PARAMS *mc_params = (TMC_PARAMS*)func_params;
        float **cen = mc_params->cen;
        float *rad = mc_params->rad;
        float num_cir = mc_params->num_cir;
        float dist,dmin;
        int   i,imin;
                                                                                
        dmin = HUGE;
        for (i = 0; i < num_cir; ++i)
        {
            dist = sqrt(sqr(coords[0] - cen[i][0]) +
                        sqr(coords[1] - cen[i][1]));
            if (dist - rad[i] < dmin)
            {
                dmin = dist - rad[i];
                imin = i;
            }
        }
        dist = sqrt(sqr(coords[0] - cen[imin][0]) +
                    sqr(coords[1] - cen[imin][1])) - rad[imin];
        return dist;
}       /* end ellipse_func */
