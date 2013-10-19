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

#include <FronTier.h>

	/*  Function Declarations */
int     FrontInitStandardIO( int,char**);
static void test_propagate(Front*);
static int trans_vel_func(POINTER,Front*,POINT*,HYPER_SURF_ELEMENT*,
                               HYPER_SURF*,float*);

char in_name[100],sc_name[100],out_name[100];
#if defined(__MPI__)
int subdomains[MAXD];
#endif /* defined(__MPI__) */

/********************************************************************
 *	Level function parameters for the initial interface 	    *
 ********************************************************************/

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

	FrontInitStandardIO(argc,argv,argc,argvf_basic);

	/* Initialize basic computational data */

	f_basic.dim = 2;	
	/*
	f_basic.L[0] = 0.0;	f_basic.L[1] = 0.0;
	f_basic.U[0] = 1.0;	f_basic.U[1] = 1.0;
	f_basic.gmax[0] = 128;	f_basic.gmax[1] = 128;
	f_basic.boundary[0][0] = f_basic.boundary[0][1] = PERIODIC_BOUNDARY;
	f_basic.boundary[1][0] = f_basic.boundary[1][1] = PERIODIC_BOUNDARY;
	*/
	f_basic.size_of_intfc_state = 0;
#if defined(__MPI__)
        f_basic.subdomains[0] = subdomains[0];
        f_basic.subdomains[1] = subdomains[1];
#endif /* defined(__MPI__) */

	f_basic.RestartRun = YES;
	f_basic.restart_name = in_name;
	printf("in_name = %s",in_name);
	FrontStartUp(&front,&f_basic);

	/* Start testing */

	/* Initialize velocity field function */

        trans_params.dim = 2;
        trans_params.vx = 0.4;
        trans_params.vy = 0.23463;

        velo_func_pack.func_params = (POINTER)&trans_params;
        velo_func_pack.func = trans_vel_func;

	FrontInitVelo(&front,&velo_func_pack);

	/* Propagate the front */

	test_propagate(&front);

	clean_up(0);
	return 0;
}


/*ARGSUSED*/
	int	FrontInitStandardIO(
	int		argc,
	char		**argv)
{
#if defined(__MPI__)
        pp_init(&argc,&argv);
#endif /* defined(__MPI__) */
	int step;
	argc--;
	argv++;
	strcpy(out_name,"intfc");
	while (argc >= 1)
	{
	    if (argv[0][0] != '-')
	    {
		printf("Usage: example -o output\n");
		exit(1);
	    }
	    switch(argv[0][1]) {
	    case 'i':
		zero_scalar(in_name,100);
		strcpy(in_name,argv[1]);
		argc -= 2;
		argv += 2;
	    	break;
	    case 't':
                step = atoi(argv[1]);
		argc -= 2;
		argv += 2;
	    	break;
	    case 'o':
		zero_scalar(out_name,100);
		strcpy(out_name,argv[1]);
#if defined(__MPI__)
                sprintf(out_name,"%s.%d",out_name,pp_mynode());
#endif /* defined(__MPI__) */
		freopen(out_name,"w",stdout);
		zero_scalar(out_name,100);
		strcpy(out_name,argv[1]);
		argc -= 2;
		argv += 2;
		break;
#if defined(__MPI__)
            case 'p':
            case 'P':
                int i,total_num_proc = 1;
                for (i = 0; i < MAXD; ++i)
                {
                    if (argc < 2 || argv[1][0] == '-') break;
                    argc -= 1;
                    argv += 1;
                    subdomains[i] = atoi(argv[0]);
                    total_num_proc *= subdomains[i];
                }
                argc -= 1;
                argv += 1;
                if (total_num_proc != pp_numnodes())
                {
                    printf("total number of processors for the partition %d "
                           "does not equal to requested np %d\n",
                           total_num_proc,pp_numnodes());
                    clean_up(ERROR);
                }
#endif /* defined(__MPI__) */
	    }
	}
	sprintf(in_name,"%s.ts%s",in_name,right_flush(step,7));
#if defined(__MPI__)
        sprintf(in_name,"%s-nd%s",in_name,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */
}

static  void test_propagate(
        Front *front)
{
        int ip,im,status,count;
        Front *newfront;
        float dt,dt_frac,CFL;
        float print_time_interval,movie_frame_interval;
        int max_step;
        bool is_print_time, is_movie_time, time_limit_reached;
        char s[10];
        float fcrds[MAXD];
        int  dim = front->rect_grid->dim;

        front->max_time = 2.2;
        front->max_step = 1000;
        front->print_time_interval = 1.0;
        front->movie_frame_interval = 0.02;

        CFL = Time_step_factor(front);

        printf("dim = %d\n", dim);

        print_front_output(front,out_name);
        show_front_output(front,out_name,binary);

	printf("original dt = %20.16f\n",front->dt);
	printf("original time = %20.16f\n",front->time);
	printf("CFL = %f\n",CFL);
	printf("Frequency_of_redistribution(front,GENERAL_WAVE) = %d\n",
		Frequency_of_redistribution(front,GENERAL_WAVE));
        time_limit_reached = NO;

        ip = (int)(front->time/print_time_interval + 1.0);
	im = (int)(front->time/movie_frame_interval + 1.0);
	printf("ip = %d  im = %d\n",ip,im);

        for (;;)
        {
	    /* Propagating interface for time step dt */

            status = FrontAdvance(front->dt,&dt_frac,front,&newfront,
                                (POINTER)NULL);
            assign_interface_and_free_front(front,newfront);

            ++front->step;
            front->time += front->dt;

	    /* Output section */

            printf("\ntime = %f   step = %5d   dt = %f\n",
                        front->time,front->step,front->dt);
            fflush(stdout);

            if (is_print_time || time_limit_reached)
                print_front_output(front,out_name);
            if (is_movie_time || time_limit_reached)
                show_front_output(front,out_name,binary);

            if (time_limit_reached)
                    break;

	    /* Time and step control section */

	    //Next time step determined by maximum speed of previous
	    //step, assuming the propagation is hyperbolic and
	    //is not dependent on second order derivatives of
	    //the interface such as curvature, and etc.

            front->dt = CFL*FrontHypTimeStep(front); 

	    front->dt = FrontOutputTimeControl(front,
	    		movie_frame_interval,
			print_time_interval,
			max_time,
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

