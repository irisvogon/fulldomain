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
*				gweno.c
*
*/


#include <ghyp/ghyp.h>
#include <gdecs/vecdecs.h>
#include "gchem.c"
#include <gdecs/geos.h>
#include <gdecs/gstate.h>
#define   v_ref                 0.001716
#define   S                     110.4
#define   density_0             0.098
#define   WALLMODEL_KAPPA       0.41
#define   WALLMODEL_APLUS      17.0
#define   SMALL                 0.00000000001
#define   MAXSLOOPSTEP           100000
#define   WALLMODEL_STRETCHING   1.1
#define   cp                    1.42
#define   Prandle               0.72
#define   Rinv                  0.0240558095// M/R = 2/8.314
#define   WALLMODEL_PRT         0.9
#define   itmp                  100000
#define   N                     1000
#define   Wall_Temp             300
#define   tk			2
#define	  tan10			0.176327	
#define	  cos10			0.984807
#define	  sin10			0.173648	
#define	  ramp_x0		 -7.0
#define   ramp_x1 		( -7 + 0.8/tan10)
#define   ramp_z0 		1.6
#define   ramp_z1 		2.4
LOCAL double wall_pressure_avg;


	/* LOCAL Function Declarations */
LOCAL	Vec_Muscl *g_load_muscl_state_data(Muscl_Opts*,int,int*,Front*,
                                                 Wave*,Stencil*,
						 Tan_stencil*,
						 Vec_Gas*,Vec_Src*,int,int,
						 int,float,float);
LOCAL	Vec_Muscl *g_muscl_alloc_wk_space(Muscl_Opts*,int,int,Vec_Muscl*);
LOCAL	bool	g_muscl_load_source_vectors(int,int*,Tan_stencil*,float*,
                                            Vec_Gas*,Vec_Src*,int,int,
					    RECT_GRID*);
LOCAL	float	*g_set_vmuscl_vec_gas_storage(float**,Vec_Gas*,Vec_Muscl*);
LOCAL	void	apply_limiter(float**,int,int,Vec_Muscl*);
LOCAL	void	g_add_art_visc1(int,int,float**,Vec_Muscl*);
LOCAL	void	g_add_art_visc2(int,int,Vec_Muscl*);
LOCAL	void	g_bct_linear_reconstructor(int,int,Vec_Muscl*);
LOCAL	void	g_compute_art_visc_coefs(int,int,struct _Vec_Muscl*);
LOCAL	void	g_compute_eigens(int,int,struct _Vec_Muscl*);
LOCAL	void	g_compute_lapidus_art_visc_coefs(int,int,Vec_Muscl*);
LOCAL	void	g_cons_src(int,int,int,int*,Tan_stencil*,struct _Vec_Muscl*);
LOCAL	void	g_cons_to_uncons(int,int,int,float**,float**);
LOCAL	void	g_eigenstate_linear_reconstructor(int,int,Vec_Muscl*);
LOCAL	void	g_first_order_godunov_reconstructor(int,int,Vec_Muscl*);
LOCAL	void	g_first_order_godunov_half_step(int,int,float,float,Vec_Muscl*);
LOCAL	void	g_half_step(int,int,float,float,Vec_Muscl*);
LOCAL	void	g_lin_rsoln(int,int,float**,float**,float**,Vec_Muscl*);
LOCAL	void	g_muscl_flux_vectors(int,int,float**,float*,MUSCL_FLUX*,
                                     Vec_Muscl*);
LOCAL	void	g_muscl_lin_approx_rsolver(int,int,float**,Vec_Gas*,float**,
                                           Vec_Gas*,float**,Vec_Gas*,
					   MUSCL_FLUX*,Vec_Muscl*);
LOCAL	void	g_state_linear_reconstructor(int,int,Vec_Muscl*);
LOCAL	void	g_strong_wave_half_step(int,int,float,float,Vec_Muscl*);
LOCAL	void	g_uncons_to_cons(int,int,int,float**,float**);
LOCAL	void	left_multiply_state_by_matrix(float**,float***,float**,
					      int,int,int);
LOCAL	void    print_eigen(int,int,Vec_Eigen*,int);
LOCAL	void    print_linear_reconstruction(int,int,Vec_Muscl*);
LOCAL	void	limit_strong_wave_in_cell(Vec_Muscl*,int,int);

/* turbulent boundary layer model */
LOCAL   void    oned_turbulence_boundarylayer_solver(int,float,float,float,float,float,float*,Locstate);
LOCAL   void compute_square_circle(double, double, double, double, double, double, double, int*,double *);

LOCAL	void	load_state_fluxes( float**, int, int, float**, float*);
LOCAL	void 	compute_weno_weights(float**, float*);
LOCAL	void	compute_weno_flux ( Front*, Vec_Gas*, int, int, int, int,int, int, int*, float**, float*, float*, float*, float**);
LOCAL	void	compute_IS ( float*, float*); 
LOCAL	void	update_pressure_soundspeed_gamma(Front*, Vec_Gas*,int,int, int, int, int, float**, float*, float*,float*);
LOCAL	void	inner_prod ( float*, float*, float*, int);
LOCAL	void	matrix_times_vector ( float*, float**, float*, int );
LOCAL	void	g_update_point_ghost_cell(bool,bool,int, int*, bool,int,float**, float*, float*,float*,int,int,int,Vec_Gas*, Front *fr);
LOCAL	void	Update_pressure_gam_c(int, int, float**,float*,float*,float*,Vec_Gas*);

LOCAL	void	load_random_nums(float**,int,int,float,float,unsigned short int[3]);

#if defined(DEBUG_MUSCL)
LOCAL	const char **set_muscl_debugging(int,int*,int*,Tan_stencil*,Front*,int);
LOCAL	const char **toggle_muscl_debugging(const char*,int,int*,int*,
					    Front*,int);
LOCAL	void	reset_muscl_debugging(const char**);
#endif /* defined(DEBUG_MUSCL) */


/*
*			oned_weno()
*
* Reference:  "Efficient Implementation of Weighted ENO Schemes"
* 		by Guang-Shan Jiang and Chi-Wang Shu
*/

#if defined(__cplusplus)
extern "C" {
#endif /* defined(__cplusplus) */
    LOCAL   LSODE_FUNC  ChemHeatSources;
//    int tag;
#if defined(__cplusplus)
}
#endif /* defined(__cplusplus) */

void    get_constant_state(Locstate s1,int comp,float *coords,float t);

/*ARGSUSED*/
EXPORT void oned_WENO(
	int		swp_num,
	int		*iperm,
	int		*icoords,
	Wave		*wave,
	Wave            *newwave,
	Front		*fr,
	Front           *newfr,
	Stencil         *sten,
	Tan_stencil	*tsten,
	int		offset,
	int		vsize,
	Vec_Gas		*vst,
	Vec_Src		*src,
	float		dt,
	float		dn,
	int		dim)
{
	int	   sten_rad;
	int        start, end;
	int        i, j, k, kmax;
	int	   num_comps;
	float 	   dtdni = dt/dn;
	
	sten_rad = 2;
	kmax = dim+2;

	/*debug lines added by xiaoxue*/
	float rhom, pm, rhoem, um,vm,wm;

        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* Compute differences for mass fraction */
            Gas_param *params = Params(vst->state[offset]);
            if((num_comps = params->n_comps) != 1)
                kmax += num_comps;
        }
	

	int   idirs[3];
	if (iperm == NULL || swp_num == -1)//This function is called by oblique_FD
		return;

	if( iperm!= NULL && swp_num != -1)
	{
	    for (i = 0; i < dim; ++i)
		idirs[i] = iperm[(i+swp_num)%dim];
	}

	if(sten == NULL) //oned_weno is called by vector_FD
	{
	    /* Load state data into vst data structure*/
	    g_load_VGas_state_vectors(offset,vsize,vst,dim);
	
	    /* States used in RK4 */
	    start = sten_rad;
	    end = vsize - sten_rad;

	    /*print out the coodinate*/
	    float *coordsxx;
	    int icoordsxx[3];
	    icoordsxx[0] = icoords[0];
	    icoordsxx[1] = icoords[1];
	    icoordsxx[2] = icoords[2];
	    coordsxx = Rect_coords(icoordsxx,wave);
//	    printf("coordsxx = [%f  %f  %f], icoordsxx = [%d  %d  %d]\n", coordsxx[0],coordsxx[1],coordsxx[2], icoordsxx[0],icoordsxx[1], icoordsxx[2]);
	
	float **u;
	float **uu1, **uu2, **F_mid;
	float *p, *GAM, *c;
	bi_array( &u , kmax, vsize+2, sizeof(float) );  // Two is added for two additional  cells
	bi_array( &uu1 , kmax, vsize+2, sizeof(float) );
	bi_array( &uu2 , kmax, vsize+2, sizeof(float) );
	bi_array( &F_mid , vsize+3, kmax, sizeof(float) );
	uni_array( &p, vsize+2, sizeof(float) );
	uni_array( &GAM, vsize+2, sizeof(float) );
	uni_array( &c, vsize+2, sizeof(float) );
	
	/* Copy states from vst to local arrays */
	for ( j = 0; j < vsize; j++){
		u[0][j+1] = vst->rho[j+offset]; 
		u[1][j+1] = vst->m[0][j+offset]; 
		u[2][j+1] = vst->m[1][j+offset]; 
		u[3][j+1] = vst->m[2][j+offset]; 
		u[4][j+1] = vst->en_den[j+offset]; 
		
		for ( i = 5; i < kmax; i++)
			u[i][j+1] = vst->rho0[i-5][j+offset]*vst->rho[j+offset];

		p[j+1]  = vst->p[j+offset];
		GAM[j+1]  = vst->GAM[j+offset];
		c[j+1]  = vst->c[j+offset];
		
	}
	
	/*Add two extra ghost states at the end */
	for ( i = 0; i < kmax; i++){	
		u[i][0] = uu1[i][0] = uu2[i][0] = u[i][1];
		u[i][vsize+1] = uu1[i][vsize+1] = uu2[i][vsize+1] = u[i][vsize];
		
		uu1[i][1] = uu2[i][1] = u[i][1];
		uu1[i][2] = uu2[i][2] = u[i][2];

		uu1[i][vsize] = uu2[i][vsize] = u[i][vsize];
		uu1[i][vsize-1] = uu2[i][vsize-1] = u[i][vsize-1];
	}
	p[0] = p[1]; p[vsize+1] = p[vsize];
	GAM[0] = GAM[1]; GAM[vsize+1] = GAM[vsize];
	c[0] = c[1]; c[vsize+1] = c[vsize];

	start = sten_rad + 1;
	end = vsize - 1;

	/* Step 1 of RK3 */
	update_pressure_soundspeed_gamma(fr, vst, 0, 0, vsize + 1, kmax, offset, u, p, c, GAM);
	compute_weno_flux (fr, vst, start, end, sten_rad, vsize, num_comps, offset,idirs, u, p, GAM, c, F_mid);
	
	for ( j = start; j < end; j++ )
		for ( i = 0; i < kmax; i++ )
			uu1[i][j] = u[i][j] + 0.5*dtdni*(F_mid[j][i] - F_mid[j+1][i]);
	for( j = start; j < end; j++)
	{
		if ( isnan(uu1[4][j]) ){
			for ( i = 0; i < kmax; i++)
				printf("u[%d][%d] = %e  Fmid[%d][%d] = %e  Fmid[%d][%d] = %e\n",i,j,u[i][j],j,i,F_mid[j][i],j+1,i,F_mid[j+1][i]);
			printf("\n");
		}
	}

	/* Step 2 of RK3 */
	update_pressure_soundspeed_gamma(fr, vst, 0, 0, vsize + 1,  kmax, offset, uu1, p, c, GAM);
	compute_weno_flux (fr, vst, start, end, sten_rad, vsize, num_comps, offset,idirs, uu1, p, GAM, c, F_mid);
	for ( j = start; j < end; j++ )
		for ( i = 0; i < kmax; i++ )
			uu2[i][j] = 0.75*u[i][j] + 0.25*uu1[i][j] + 0.25*dtdni*(F_mid[j][i] - F_mid[j+1][i]);
	for( j = start; j < end; j++)
	{
		if ( isnan(uu2[4][j]) ){
			for ( i = 0; i < kmax; i++)
				printf("uu1[%d][%d] = %e  Fmid[%d][%d] = %e  Fmid[%d][%d] = %e\n",i,j,uu1[i][j],j,i,F_mid[j][i],j+1,i,F_mid[j+1][i]);
			printf("\n");
		}
	}
	free(uu1);

	/* Step 3 of RK3 */
	update_pressure_soundspeed_gamma(fr, vst, 0, 0, vsize + 1, kmax, offset, uu2, p, c, GAM);
	compute_weno_flux (fr, vst, start, end, sten_rad, vsize, num_comps, offset,idirs,uu2, p, GAM, c, F_mid);


	for ( j = start; j < end; j++ )
		for ( i = 0; i < kmax; i++ )
			u[i][j] = (1.0/3.0)*u[i][j] + (2.0/3.0)*uu2[i][j] + (2.0/3.0)*dtdni*(F_mid[j][i] - F_mid[j+1][i]);

	for( j = start; j < end; j++)
	{
		if ( isnan(u[4][j]) ){
			for ( i = 0; i < kmax; i++)
				printf("uu2[%d][%d] = %e  Fmid[%d][%d] = %e  Fmid[%d][%d] = %e\n",i,j,uu2[i][j],j,i,F_mid[j][i],j+1,i,F_mid[j+1][i]);
			printf("\n");
		}
	}	

free(uu2);

	/* Copy state values back to vst */
	start = sten_rad;
	end = vsize - sten_rad;
	for ( j = start; j < end; j++){
		vst->rho[j+offset] = u[0][j+1]; 
		vst->m[0][j+offset] = u[1][j+1]; 
		vst->m[1][j+offset] = u[2][j+1]; 
		vst->m[2][j+offset] = u[3][j+1]; 
		vst->en_den[j+offset] = u[4][j+1];
		
		for ( i = 0; i < num_comps; i++ )
			vst->rho0[i][j+offset] = u[5+i][j+1];
	}

	    free(u);
	    free(F_mid );
	    free(p);
	    free(GAM);
	    free(c);
	}
	else //oned_WENO is called by point_FD
	{
//	    printf("point_FD!\n");
//	    printf("coords = [%f %f %f], icoords = [%f  %f %f]\n",Coords(sten->pstore[2])[0], Coords(sten->pstore[0])[1],Coords(sten->pstore[0])[2],  sten->icoords[0][0],sten->icoords[0][1],sten->icoords[0][2]);
	    /* Whether we are calling point_FD at the upboundary*/
	    bool upboundary;
	    int **icoords = sten->icoords;
	    int icoordsp[3];
	    icoordsp[0] = icoords[0][0];
	    icoordsp[1] = icoords[0][1];
	    icoordsp[2] = icoords[0][2];
	    upboundary = YES;
	    int comp = Rect_comp(icoordsp,wave);
	    int sten_comp[5];
	    for( i = 0; i < 5; i++)
	    {
		icoordsp[idirs[0]] = icoords[i-2][idirs[0]];
		sten_comp[i] = Rect_comp(icoordsp,wave);
	    }
	    if( sten_comp[3] == sten_comp[4] && sten_comp[3] + sten_comp[4] == 6)
		upboundary = 0;
	    else if( sten_comp[0] == sten_comp[1] && sten_comp[0] + sten_comp[1] == 6)
		upboundary = 1;
	    else	
	    printf(" the componet of the point_FD! comp = [%d  %d  %d  %d  %d]\n", sten_comp[0], sten_comp[1], sten_comp[2], sten_comp[3],sten_comp[4]);
	
	    bool inflow = 0;
	    bool nozzle = 0;
	    int nozzle_indx;
	    int inflow_indx;
	    int ubuff = 3;
	    int stencil_length = vsize + 3*ubuff - 2;
	    if(sten_comp[0] == 7)
	    {
		inflow = 1;
		inflow_indx = 1;
	    }
	    else if(sten_comp[1] == 7)
	    {
		inflow = 1;
		inflow_indx = 2;
	    }  

	    if(vst->coords[stencil_length - 2][2] > 3.91 ||vst->coords[stencil_length - 3][2] > 3.91)
	    {
		nozzle = YES;
		if(vst->coords[stencil_length - 3][2] > 3.91) nozzle_indx = 2;
		else nozzle_indx = 3;
		//printf("nozzle_indx = %d\n",nozzle_indx);
	    }
	      
	    float *coordsxx = vst->coords[9];
	    
	    /*Load state varibles into uu */
	  //  int  stencil_length; // num of state in vst
	    int  u_length; //num of state in u
	  //  int  ubuff = 3;
	  //  stencil_length = vsize + 3*ubuff - 2;

	    /* Load state data into vst data structure*/
	    g_load_VGas_state_vectors(offset,stencil_length,vst,dim);

	    /*Assign the space of the 3-order weno stencil*/
	    int num_ghost_cell = 3; //num of ghost cells in u,u1,u2;
	    int ifwall; // index of the first cell off the wall in vst;
	    int num_cells_in_wall = 2;//num of the ghost cells in vst; 
	    bool boundary_cell = 0;
	    int	 Fmid_indx;
	    if(inflow)
	    {
		num_cells_in_wall = inflow_indx;
		ifwall = inflow_indx;
	    }
	    else if(nozzle)
	    {
		num_cells_in_wall = 4 - nozzle_indx;
		ifwall = 7 + nozzle_indx - 1; 
	//	printf("num_cells_in_wall = %d\n",num_cells_in_wall);
	    }
	    else if(!upboundary)
	    {
		if (Reflect_wall(vst->state[1]) || Ramp_cellz(vst->state[1]))//point_FD update the second cell off the wall
		{ 
		    num_cells_in_wall = 1;
		    ifwall = 1;
		    if(Reflect_wall(vst->state[1]))
		    boundary_cell = 1;
		    Fmid_indx = ifwall;
		}
		else if(Reflect_wall(vst->state[2]) || Ramp_cellz(vst->state[2]))//point_FD update the first cell off the wall
		{
		    num_cells_in_wall = 2;		    
		    ifwall = 2;
		    if(Reflect_wall(vst->state[2]))  boundary_cell = 1;
		    Fmid_indx = ifwall;
		}
		else
		{
		    printf("error! in point_FD in vst stencil in lowwwall\n");
		    /*for( i = 0; i < 5; i++)
			printf("coords[%d] = [%f  %f  %f]\n", i, vst->coords[i][0], vst->coords[i][1],vst->coords[i][2]);
		    for( i = 0; i < stencil_length; i++)printf("Nozzle[%d] = %d  ", i, Nozzle(vst->state[i]));printf("\n");
		    printf(" Reflect state = [%d  %d  %d  %d  %d] Ramp_cellz = [%d  %d  %d  %d  %d]\n",Reflect_wall(vst->state[0]),Reflect_wall(vst->state[1]),Reflect_wall(vst->state[2]),Reflect_wall(vst->state[3]),Reflect_wall(vst->state[4]),Ramp_cellz(vst->state[0]),Ramp_cellz(vst->state[1]),Ramp_cellz(vst->state[2]),Ramp_cellz(vst->state[3]),Ramp_cellz(vst->state[4]));*/
		} 
	    }
	    else
	    {
		if (Reflect_wall(vst->state[stencil_length - 1 - 1]) || Ramp_cellx(vst->state[stencil_length - 1 - 1]))//point_FD update the second cell off the wall
		{
		    num_cells_in_wall = 1;
		    ifwall = stencil_length - 1 - 1;
		    if(Reflect_wall(vst->state[stencil_length - 1 - 1]))  boundary_cell = 1;
		    Fmid_indx = ifwall + 1;
		}
		else if(Reflect_wall(vst->state[stencil_length - 1 - 2]) || Ramp_cellx(vst->state[stencil_length - 1 - 2]))//point_FD update the first cell off the wall
		{
		    num_cells_in_wall = 2;
		    ifwall = stencil_length - 1 - 2;
		    if(Reflect_wall(vst->state[stencil_length - 1 - 2]))boundary_cell = 1;
		    Fmid_indx = ifwall + 1;
		}
		else
		{
		    printf("error! in point_FD in vst stencil in upwall\n");
		    //for( i = 0; i < stencil_length; i++)printf("Nozzle[%d] = %d  ", i, Nozzle(vst->state[i]));printf("\n");
		    for( i = 0; i < 5; i++)
			 printf("coords[%d] = [%f  %f  %f]\n", i, vst->coords[stencil_length - 1 - i][0], vst->coords[stencil_length - 1 - i][1],vst->coords[stencil_length - 1 - i][2]);
		    printf(" Reflect state = [%d  %d  %d  %d  %d] Ramp_cellx = [%d  %d  %d  %d  %d]\n",Reflect_wall(vst->state[stencil_length - 1]),Reflect_wall(vst->state[stencil_length - 1 -1 ]),Reflect_wall(vst->state[stencil_length - 1 - 2]),Reflect_wall(vst->state[stencil_length - 1 - 3]),Reflect_wall(vst->state[stencil_length - 1 - 4]),Ramp_cellx(vst->state[stencil_length - 1]),Ramp_cellx(vst->state[stencil_length - 1 - 1]),Ramp_cellx(vst->state[stencil_length - 1 - 2]),Ramp_cellx(vst->state[stencil_length - 1 - 3]),Ramp_cellx(vst->state[stencil_length - 1 - 4]));
		} 
	    }
	    
	    u_length = vsize + 3*ubuff - 2 - num_cells_in_wall + num_ghost_cell;
//	    printf(" u_length =  %d \n", u_length);
		
	    float **u, **uu1, **uu2, **F_mid;
	    float *p, *GAM, *c;
	    bi_array( &u ,   kmax, u_length, sizeof(float));
	    bi_array( &uu1 , kmax, u_length, sizeof(float));
	    bi_array( &uu2 , kmax, u_length, sizeof(float));
	    bi_array( &F_mid ,u_length + 1, kmax, sizeof(float));
	    uni_array( &p, u_length, sizeof(float));
	    uni_array( &GAM, u_length, sizeof(float));
	    uni_array( &c,u_length, sizeof(float));

	    int  startx[3],endx[3], indx_ans_u, indx_ans_vst,fcell,lcell;
	    int  startp[3],endp[3];
	    
	    if(!nozzle && ( inflow || !upboundary))
	    {
		fcell = num_ghost_cell;
		lcell = u_length -1;
		/*Load state data from vst into u*/
	    	for( j = 0 ; j < stencil_length - num_cells_in_wall; j++)
	    	{
		    u[0][j + num_ghost_cell] = vst->rho[j + num_cells_in_wall];
		    u[1][j + num_ghost_cell] = vst->m[0][j + num_cells_in_wall];
		    u[2][j + num_ghost_cell] = vst->m[1][j + num_cells_in_wall];
		    u[3][j + num_ghost_cell] = vst->m[2][j + num_cells_in_wall];
		    u[4][j + num_ghost_cell] = vst->en_den[j + num_cells_in_wall];
		    for ( i = 5; i < kmax; i++)
			u[i][j + num_ghost_cell] = vst->rho0[i-5][j + num_cells_in_wall]*vst->rho[j + num_cells_in_wall];
		    p[j + num_ghost_cell] = vst->p[j + num_cells_in_wall];
		    GAM[j + num_ghost_cell] = vst->GAM[j + num_cells_in_wall] ;
		    c[j + num_ghost_cell] = vst->c[j + num_cells_in_wall];
	    	} 

		startx[0] = num_ghost_cell;  endx[0] = u_length - ubuff;
		startx[1] = num_ghost_cell;  endx[1] = u_length - 2*ubuff;
		startx[2] = num_ghost_cell;  endx[2] = u_length - 3*ubuff;

		startp[0] = 0; endp[0] = u_length;
		startp[1] = 0; endp[1] = u_length - ubuff;
		startp[2] = 0; endp[2] = u_length - 2*ubuff;

		indx_ans_u =  num_ghost_cell + 2 - num_cells_in_wall;
		indx_ans_vst = 2;
	    }
	    else 
	    {
		fcell = 0;
		lcell = u_length - num_ghost_cell - 1;
		/*Load state data into u*/
	    	for( j = 0 ;j < stencil_length - num_cells_in_wall; j++)
	    	{
		    u[0][j] = vst->rho[j];
		    u[1][j] = vst->m[0][j];
		    u[2][j] = vst->m[1][j];
		    u[3][j] = vst->m[2][j];
		    u[4][j] = vst->en_den[j];
		    for ( i = 5; i < kmax; i++)
			u[i][j] = vst->rho0[i-5][j]*vst->rho[j];
		    p[j] = vst->p[j];
		    GAM[j] = vst->GAM[j] ;
		    c[j] = vst->c[j];
	    	}

		startx[0] = ubuff;    endx[0] = u_length - num_ghost_cell;
		startx[1] = 2*ubuff;  endx[1] = u_length - num_ghost_cell;
		startx[2] = 3*ubuff;  endx[2] = u_length - num_ghost_cell;
	
		startp[0] = 0;		endp[0] = u_length;
		startp[1] = ubuff; 	endp[1] = u_length;
		startp[2] = 2*ubuff; 	endp[2] = u_length;

		indx_ans_u =  3*ubuff;
		indx_ans_vst = 3*ubuff;
	    }

		/* Step 1 of RK3 */
		g_update_point_ghost_cell(nozzle,inflow, swp_num,iperm,upboundary, kmax,u, p, GAM, c, num_ghost_cell,fcell,lcell,vst,fr);
		update_pressure_soundspeed_gamma(fr, vst, indx_ans_vst, startp[0], endp[0], kmax, offset, u, p, c, GAM);
		for( i = startp[0]; i < endp[0];i++) 
		{
		    float sum_frac = 0.0;
		    for( j = 5; j < kmax; j++) sum_frac += u[j][i];
		    if( fabs(sum_frac - u[0][i]) > 100*MACH_EPS) 
		    {
			    printf("PARTIAL DENS not consistent1!\n");
			    printf(" sum_frac = %f u[0][%d] = %f  ", sum_frac, i, u[0][i]);
			    for( j = 5; j < kmax; j++) printf("%g  ", u[j][i]);
			    printf("\n");
		    }
		}
		compute_weno_flux(fr,vst,startx[0],endx[0],sten_rad,vsize,num_comps,offset,idirs, u, p, GAM, c, F_mid);

		if(boundary_cell){for(i = 0; i < kmax; i++){if( i!= 1) F_mid[Fmid_indx][i] = 0;}}

		for( j = startx[0]; j < endx[0] + 1; j++)
		{
			float sum_flux = 0.0;
			for( i = 5; i < kmax; i++) sum_flux += F_mid[j][i];
			if( fabs( sum_flux - F_mid[j][0]) > 1000*MACH_EPS)
			{
				printf("partial flux is not consistant!1\n");
				printf(" sum_flux = %f F_mid[j][0] = %f  ", sum_flux, F_mid[j][0]);
				for( i = 0; i < kmax; i++) printf("%g  ", F_mid[j][i]);
					printf("\n");
			}
		}

	for( j = startx[0]; j < endx[0]+1; j++)
	{
		rhom = 0.5*(u[0][j-1] + u[0][j]);
		pm = 0.5*(p[j-1] + p[j]);
		rhoem = 0.5*(u[4][j-1] + u[4][j]);
		um = 0.5*(u[1][j-1]/u[0][j-1] + u[1][j]/u[0][j]);
		vm = 0.5*(u[2][j-1]/u[0][j-1] + u[2][j]/u[0][j]);
		wm = 0.5*(u[3][j-1]/u[0][j-1] + u[3][j]/u[0][j]);
	}
		for ( j = startx[0]; j < endx[0]; j++)
		    for ( i = 0; i < kmax; i++)
			uu1[i][j] = u[i][j] + 0.5*dtdni*(F_mid[j][i] - F_mid[j+1][i]);

		/* Step 2 of RK3 */
	   	g_update_point_ghost_cell(nozzle,inflow,swp_num,iperm,upboundary,kmax,uu1,p,GAM,c,num_ghost_cell, fcell,lcell,vst,fr);
		update_pressure_soundspeed_gamma(fr, vst, indx_ans_vst,startp[1],endp[1], kmax, offset, uu1, p, c, GAM);
		for( i = startp[1]; i < endp[1];i++) 
		{
		    float sum_frac = 0.0;
		    for( j = 5; j < kmax; j++) sum_frac += uu1[j][i];
		    if( fabs(sum_frac - uu1[0][i]) > 100*MACH_EPS) 
		    {
			    printf("PARTIAL DENS not consistent2!\n");
			    printf(" sum_frac = %f uu1[0][%d] = %f  ", sum_frac, i, uu1[0][i]);
			    for( j = 5; j < kmax; j++) printf("  %g", uu1[j][i]);
			    printf("\n");
		    }
		}


		compute_weno_flux(fr, vst,startx[1],endx[1], sten_rad, vsize, num_comps, offset, idirs, uu1, p, GAM, c, F_mid);
		if(boundary_cell){for(i = 0; i < kmax; i++){if( i!= 1) F_mid[Fmid_indx][i] = 0;}}
		for( j = startx[1]; j < endx[1] + 1; j++)
		{
			float sum_flux = 0.0;
			for( i = 5; i < kmax; i++) sum_flux += F_mid[j][i];
			if( fabs( sum_flux - F_mid[j][0]) > 1000*MACH_EPS)
			{
				printf("partial flux is not consistant!2\n");
				printf(" sum_flux = %f F_mid[j][0] = %f ", sum_flux, F_mid[j][0]);
				for( i = 0; i < kmax; i++) printf(" %g ", F_mid[j][i]);
				printf("\n");
			}
		}

	    	for( j = startx[0]; j < endx[1]+1; j++)
		{
			rhom = 0.5*(uu1[0][j-1] + uu1[0][j]);
			pm = 0.5*(p[j-1] + p[j]);
			rhoem = 0.5*(uu1[4][j-1] + uu1[4][j]);
			um = 0.5*(uu1[1][j-1]/uu1[0][j-1] + uu1[1][j]/uu1[0][j]);
			vm = 0.5*(uu1[2][j-1]/uu1[0][j-1] + uu1[2][j]/uu1[0][j]);
			wm = 0.5*(uu1[3][j-1]/uu1[0][j-1] + uu1[3][j]/uu1[0][j]);
		}

		for ( j = startx[1]; j < endx[1]; j++)
		    for ( i = 0; i < kmax; i++)
				uu2[i][j] = 0.75*u[i][j] + 0.25*uu1[i][j] + 0.25*dtdni*(F_mid[j][i] - F_mid[j+1][i]);

		/* Step 3 of RK3 */
	   	g_update_point_ghost_cell(nozzle,inflow,swp_num, iperm, upboundary, kmax,uu2, p, GAM, c, num_ghost_cell, fcell,lcell, vst,fr);
			update_pressure_soundspeed_gamma(fr, vst, indx_ans_vst, startp[2], endp[2], kmax, offset, uu2, p, c, GAM);
	    	for( i = startp[2]; i < endp[2];i++) 
		{
		    float sum_frac = 0.0;
		    for( j = 5; j < kmax; j++) sum_frac += uu2[j][i];
		    if( fabs(sum_frac - uu2[0][i]) > 100*MACH_EPS) 
		    {
			    printf("PARTIAL DENS not consistent3!\n");
			    printf(" sum_frac = %f uu2[0][%d] = %f  ", sum_frac, i, uu2[0][i]);
			    for( j = 5; j < kmax; j++) printf("  %g", uu2[j][i]);
			    printf("\n");
		    }
		}
		compute_weno_flux (fr, vst,startx[2],endx[2], sten_rad, vsize, num_comps, offset, idirs, uu2, p, GAM, c, F_mid);
		if(boundary_cell){for(i = 0; i < kmax; i++){if( i!= 1) F_mid[Fmid_indx][i] = 0;}}
		for( j = startx[2]; j < endx[2] + 1; j++)
		{
			float sum_flux = 0.0;
			for( i = 5; i < kmax; i++) sum_flux += F_mid[j][i];
			if( fabs( sum_flux - F_mid[j][0]) > 1000*MACH_EPS)
			{
				printf("partial flux is not consistant!2\n");
				printf(" sum_flux = %f F_mid[j][0] = %f ", sum_flux, F_mid[j][0]);
				for( i = 0; i < kmax; i++) printf("  %g", F_mid[j][i]);
				printf("\n");
			}
		}

	    	for( j = startx[0]; j < endx[2]+1; j++)
		{
			rhom = 0.5*(uu2[0][j-1] + uu2[0][j]);
			pm = 0.5*(p[j-1] + p[j]);
			rhoem = 0.5*(uu2[4][j-1] + uu2[4][j]);
			um = 0.5*(uu2[1][j-1]/uu2[0][j-1] + uu2[1][j]/uu2[0][j]);
			vm = 0.5*(uu2[2][j-1]/uu2[0][j-1] + uu2[2][j]/uu2[0][j]);
			wm = 0.5*(uu2[3][j-1]/uu2[0][j-1] + uu2[3][j]/uu2[0][j]);
		}	

		for ( j = startx[2]; j < endx[2]; j++)
		    for ( i = 0; i < kmax; i++)
			u[i][j] = (1.0/3.0)*u[i][j] + (2.0/3.0)*uu2[i][j] + (2.0/3.0)*dtdni*(F_mid[j][i] - F_mid[j+1][i]);

/*		for( i = startx[2]; i < endx[2];i++) 
		{
		    float sum_frac = 0.0;
		    for( j = 5; j < kmax; j++) sum_frac += u[j][i];
		    if( fabs(sum_frac - 1.0) > 100*MACH_EPS) 
		    {
			    printf("PARTIAL DENS not consistent4!\n");
			    printf(" sum_frac = %f ", sum_frac);
			    for( j = 5; j < kmax; j++) printf("%g", u[j][i]);
			    printf("\n");
		    }
		}
*/
		
		free(uu1);
		free(uu2);

//		printf("sweep direction:%d\n", idirs[0]);
//		printf("old: in_u = %d in_ans_vst = %d rho = %f u = [%f  %f  %f] Energy = %f\n",indx_ans_u, indx_ans_vst, vst->rho[indx_ans_vst],vst->m[0][indx_ans_vst]/vst->rho[indx_ans_vst], vst->m[1][indx_ans_vst]/vst->rho[indx_ans_vst] ,vst->m[2][indx_ans_vst]/vst->rho[indx_ans_vst] ,vst->en_den[indx_ans_vst]);
//		printf("changes: in_u = %d in_ans_vst = %d rho = %f u = [%f  %f  %f] Energy = %f\n",indx_ans_u, indx_ans_vst, u[0][indx_ans_u] - vst->rho[indx_ans_vst],(u[1][indx_ans_u] -vst->m[0][indx_ans_vst])/vst->rho[indx_ans_vst], (u[2][indx_ans_u] - vst->m[1][indx_ans_vst])/vst->rho[indx_ans_vst] ,(u[3][indx_ans_u] -vst->m[2][indx_ans_vst])/vst->rho[indx_ans_vst] ,(u[4][indx_ans_u] - vst->en_den[indx_ans_vst]));
		vst->rho[indx_ans_vst]= u[0][indx_ans_u] ; 
		vst->m[0][indx_ans_vst] = u[1][indx_ans_u] ;
		vst->m[1][indx_ans_vst] = u[2][indx_ans_u] ;
		vst->m[2][indx_ans_vst] = u[3][indx_ans_u];
		vst->en_den[indx_ans_vst]  = u[4][indx_ans_u];
		for ( i = 5; i < kmax; i++)
		    vst->rho0[i-5][indx_ans_vst] = u[i][indx_ans_u];
		start = indx_ans_vst;
		end = indx_ans_vst + 1;
//		printf("new: in_u = %d in_ans_vst = %d rho = %f u = [%f  %f  %f] Energy = %f\n",indx_ans_u, indx_ans_vst, vst->rho[indx_ans_vst],vst->m[0][indx_ans_vst]/vst->rho[indx_ans_vst], vst->m[1][indx_ans_vst]/vst->rho[indx_ans_vst] ,vst->m[2][indx_ans_vst]/vst->rho[indx_ans_vst] ,vst->en_den[indx_ans_vst]);
		free(u);
		free(F_mid );
		free(c);
		free(p);
		free(GAM);
	}	

	bool *Flag;
	uni_array(&Flag, kmax + 1, sizeof(bool));
	for (i = 0; i < kmax; ++i)
		Flag[i] = YES;
//	    Flag[i] = NO;

	if((swp_num == 0) || (swp_num == 1) || (swp_num == 2))
	    for (j = start; j < end; ++j)
	    {
		for (i = 0; i < kmax; ++i)
		{
		    //if (vmuscl->du[i][j] != 0)
		    //	Flag[i] = YES;
			Flag[i] = NO;
		    Lflag(vst->state[j])[i] = Flag[i];
		}
	    }
	free(Flag);
	
	//xiaoxue
	/*
	int   idirs[3];
	if((iperm != NULL && swp_num == -1) && (iperm == NULL && swp_num != -1)) 
		printf("error in iperm.");

	if( iperm!= NULL && swp_num != -1){
	     for (i = 0; i < dim; ++i)
	     	idirs[i] = iperm[(i+swp_num)%dim];
	}
	else{
	    for (i = 0; i < dim; ++i) idirs[i] = i;
	}

	for ( j = start; j < end; j++){
		if( Merge_cell(vst->state[j])  &&fr->turb_boundary_layer == YES && iperm != NULL && swp_num != -1 ){
	    	i = 2*iperm[swp_num%3]; //sweep direction.
	    	//density flux
	    	Flux(vst->state[j])[i][0] = F_mid[j][0];
	    	Flux(vst->state[j])[i+1][0] = F_mid[j+1][0];
	    	flux_flag(vst->state[j]) = 1;
	    	//energy flux
	    	Flux(vst->state[j])[i][1] = F_mid[j][4];
	    	Flux(vst->state[j])[i+1][1] = F_mid[j+1][4];
	    	//momentum flux
	    	for(k = 0; k < 3; k++){
	    		Flux(vst->state[j])[i][idirs[k]+2]   = F_mid[j][k+1];
	    		Flux(vst->state[j])[i+1][idirs[k]+2] = F_mid[j+1][k+1];
	    	}
	    	//partial density flux
	    	for(k = 5; k < kmax; k++){
			Flux(vst->state[j])[i][k]   = F_mid[j][k];
			Flux(vst->state[j])[i+1][k] = F_mid[j+1][k];
	    	}
		}
		else{

			for ( k = 0; k < kmax; k++ ){
				u[k][j] += dtdni * ( F_mid[j][k] - F_mid[j+1][k]);
			}
		}
	}*/
	
	/*  Solve for and add chemistry source term */
        if( (swp_num == 2) /*&& (sten==NULL)*/)
	{
	    double *Y, *yp;
	    float *coords;
	    int istate;
	    int ic[3];
	    ic[idirs[0]] = icoords[idirs[0]];
	    ic[idirs[1]] = icoords[idirs[1]];
	    ic[idirs[2]] = icoords[idirs[2]];
	    for(j = start; j < end; ++j)
	    {
		ic[idirs[0]] = j;
		coords = Rect_coords(ic,wave);
		float T = temperature(vst->state[j]);
		if(is_obstacle_state(vst->state[j])) printf("ERROR ! coords = [%f  %f  %f]\n",coords[0],coords[1],coords[2]);

		if (T>=200 && T <=3500 && !is_obstacle_state(vst->state[j]))
		{
                	uni_array(&Y,num_comps+5,sizeof(double));

			Y[num_comps+1] = vst->en_den[j];  //set initial energy for ODE
			for(k = 0; k < num_comps; k++)
				Y[1+k] = vst->rho0[k][j];
			Y[num_comps+2] = vst->rho[j];  //set initial density for ODE
			Y[num_comps+3] = T;  //set initial temperature for ODE
//			Y[num_comps+4] = vst->p[j];  //set initial pressure for ODE
			Y[num_comps+4] = pressure(vst->state[j]);	
	                istate = 1;
			
			/*if ( (fabs(coords[0]-20.277993)< 1e-3) &&
			     (fabs(coords[1]+0.781250)<1e-3) &&
			     (fabs(coords[2]-3.427083)<1e-3))
			{
				for ( k =0; k < 15; k++) printf("VM: Y[%d] = %e\n",k,Y[k]);
				printf("c = %e\n",vst->c[j]);
				printf("\n");
			}*/
			
                	if (!ode_solver(ChemHeatSources,Y,0.0,dt,1,1.e-8,
                	                 10*MACH_EPS,
                	                 1.e-4,&istate,NULL))
                	{
			    printf("Before ODE crashes, energy = %E\n", vst->en_den[j]);
			    for(k=0;k<10;k++)
				printf("partial density %d = %E\n",k,pdens(vst->state[j])[k]);
			    if ( sten == NULL ) printf("Called by vector_FD\n");
			    printf("T = %f\n",temperature(vst->state[j]));
	 		    printf("P = %f\n",pressure(vst->state[j]));
			    printf("After ODE crashes, energy = %E\n", Y[NUM_OF_SPECIES+1]);
			    printf("Temperature: %E\n", Y[num_comps+3]);
			    printf("Pressure: %E\n", Y[num_comps+4]);
			    printf("IE = %f\n",internal_energy(vst->state[j]));
		            printf("iCoords: (%d, %d, %d)\n", icoords[0], icoords[1], icoords[2]);
		    	    printf("Coords: (%f, %f, %f)\n", coords[0], coords[1], coords[2]);
		    	    printf("j = %d\n",j);
		    	    printf("start = %d\n",start);
		    	    printf("end = %d\n",end);
		   	    screen("ERROR in oned_WENO(), "
                           	"ode_solver() failed for ChemHeatSources\n");
                            clean_up(ERROR);
                	}

			vst->en_den[j] = Y[num_comps + 1];
                	vst->rho[j]= Y[num_comps+2];   
                	for(k = 0; k < num_comps; k++)
                	    vst->rho0[k][j] = Y[k+1];

                	free(Y);
		}
            } 
        }
}		/*end oned_weno*/

LOCAL	void	g_update_point_ghost_cell(
		bool nozzle,
		bool inflow,
		int swp_num,
		int *iperm,
		bool upboundary,
		int kmax,
		float	**u,
		float	*p,
		float	*GAM,
		float	*c,
		int num_ghost_cell,
		int start,//index of the fisrt interior point in **u 
		int end,//index of the last interior point in **u 
		Vec_Gas *vst,
		Front	*fr)
{
//	printf(" in update ghost cell inflow = %d nozzle = %d upboundary = %d num_ghost_cell = %d start = %d end = %d\n", inflow, nozzle, upboundary,  num_ghost_cell, start, end);
	int i,j,k;
	float	**coords;
	coords = vst->coords;
	float	*coordsx;
	int dim = 3;
	int idirs[3];
	for( i = 0; i < 3; i++)
	    idirs[i] = iperm[(i + swp_num)%dim];
	if(inflow && !upboundary)
	{
		for( i = 0; i < num_ghost_cell; i++)
		{
		    for( j = 0 ; j < kmax; j++)
		    	u[j][i] = u[j][num_ghost_cell];
		}
	        return;
	}

	if(nozzle)
	{
		//printf("nozzle! end + 1 = %d end + num_ghost_cell = %d\n", end + 1, end + num_ghost_cell);
		Locstate sl = NULL;
		if ( sl == NULL )  g_alloc_state(&sl, fr->sizest);
		set_type_of_state(sl,TGAS_STATE);
		Press(sl) = 12.5;
		Dens(sl) = 12.3e-4;
		Vel(sl)[0] = 0.0;
		Vel(sl)[1] = 0.0;
		Vel(sl)[2] = -123.0;
		for(i = 0; i < 10; i++)
			pdens(sl)[i] = 0.0;
		pdens(sl)[4] = Dens(sl);
		Set_params(sl,vst->state[9]);	
		int st_type = state_type(vst->state[9]);
		set_state(sl, st_type, sl);

		for( i = end + 1; i < end + num_ghost_cell + 1; i++)
		{
		      u[0][i] = Dens(sl);
		      for( j = 1; j < 4; j++)u[j][i] = Mom(sl)[idirs[j-1]];
		      u[4][i] = Energy(sl);
		      for( j = 5; j < kmax; j++) u[j][i] = pdens(sl)[j-5];
//		      printf(" i = %d  u = %f  %f  %f  %f  %f\n", i, u[0][i], u[1][i], u[2][i], u[3][i], u[4][i]);
		}
		return;
	}
	if(inflow && upboundary) printf("ERROR!\n");
	    float nvec[3], nnvec[3];
	    float u_times_nnvec ;
//	    printf("idirs = [%d  %d  %d]  iperm = [%d  %d  %d] swp_num = %d \n", idirs[0], idirs[1], idirs[2], iperm[0], iperm[1], iperm[2], swp_num);
	if(idirs[0] == 0 && upboundary) // We need to assign the ghost cell inside the ramp for the x direction sweep
	{
	    if(coords[0][2] > 2.5) printf("Error in g_update_point_ghost_cell!\n");
	    nvec[0] = -sin10;
	    nvec[1] = 0.0;
	    nvec[2] = cos10;
//	    printf("ramp upboundary!\n");
	    for( j = 0; j < dim; j++)	nnvec[j] = nvec[idirs[j]];
	    for( i = 0; i < num_ghost_cell; i++)
	    {
	//	printf(" 1.end + 1 + i = %d \n", end + 1 + i);
		u[0][end + 1 + i] = u[0][end - i]; //rho
		u[4][end + 1 + i] = u[4][end - i]; //en_den
		u_times_nnvec = 0;
		for( j = 0; j < 3; j++)
		    u_times_nnvec += u[1 + j][end - i]*nnvec[j];
		for( j = 0; j < 3; j++)
		    u[1 + j][end + 1+ i] = u[1 + j][end - i] - 2 * nnvec[j]*u_times_nnvec;
		for( j = 5; j < kmax;j++)    u[j][end + 1 + i] = u[j][end - i];
		p[end + 1 + i] = p[end - i];GAM[end +1 + i] = GAM[end - i];c[end + 1 + i] = GAM[end - i]; 
//		for( j = 0; j < 3; j++)
//	    	printf("u[%d][%d] = %f  u[%d][%d] = %f\n", 1 + j ,end + 1+ i, u[1 + j][end + 1+ i], 1 + j, end - i ,u[1 + j][end - i]);
//		printf("nvec = [%f  %f  %f] nnvec = [%f  %f  %f] u_times_nnvec = %f \n", nvec[0], nvec[1], nvec[2], nnvec[0], nnvec[1], nnvec[2],u_times_nnvec);
	    }	   
	}
	else if(idirs[0] == 2 && !upboundary && ( coords[0][0] < ramp_x0 || coords[0][0] > ramp_x1))
	{
	    nvec[0] = 0;
	    nvec[1] = 0;
	    nvec[2] = 1;
	
	    for( j = 0; j < dim; j++)	nnvec[j] = nvec[idirs[j]];
	    for( i = 0; i < num_ghost_cell; i++)
	    {
//		 printf(" lowboundary! \n");
		u[0][num_ghost_cell - 1 - i] = u[0][num_ghost_cell + i]; //rho
		u[4][num_ghost_cell - 1 - i] = u[4][num_ghost_cell + i]; //en_den
//		printf("2. num_ghost_cell - 1 - i  = %d nnvec = [%f  %f  %f]\n", num_ghost_cell - 1 - i, nnvec[0],nnvec[1],nnvec[2]);
		u_times_nnvec = 0.0;
		for( j = 0; j < 3; j++)
		    u_times_nnvec += u[1 + j][num_ghost_cell + i]*nnvec[j];
		for( j = 0; j < 3; j++)
		    u[1 + j][num_ghost_cell - 1 - i] = u[1 + j][num_ghost_cell + i] - 2 * nnvec[j]*u_times_nnvec;
//		for( j = 0; j < 3; j++)
//			printf("u[%d][%d] = %f  u[%d][%d] = %f\n", 1 + j ,num_ghost_cell - 1 - i, u[1 + j][num_ghost_cell - 1 - i], 1 + j, num_ghost_cell + i,u[1 + j][num_ghost_cell + i]);
//		printf("nvec = [%f  %f  %f] nnvec = [%f  %f  %f] u_times_nnvec = %f \n",nvec[0], nvec[1], nvec[2],  nnvec[0], nnvec[1], nnvec[2],u_times_nnvec);
		for( j = 5; j < kmax;j++)    u[j][num_ghost_cell - 1 - i] = u[j][num_ghost_cell + i];
		p[num_ghost_cell - 1 - i] = p[num_ghost_cell + i];
		GAM[num_ghost_cell - 1 - i] = GAM[num_ghost_cell + i];
		c[num_ghost_cell - 1 - i] = c[num_ghost_cell + i];
	    }
	    
	}
	else if(idirs[0] == 2 && !upboundary && coords[0][0] > ramp_x0 && coords[0][0] < ramp_x1)
	{
	    nvec[0] = -sin10;
	    nvec[1] = 0.0;
	    nvec[2] = cos10;
	    for( j = 0; j < dim; j++)	nnvec[j] = nvec[idirs[j]];
	    for( i = 0; i < num_ghost_cell; i++)
	    {
		u[0][num_ghost_cell - 1 - i] = u[0][num_ghost_cell + i]; //rho
		u[4][num_ghost_cell - 1 - i] = u[4][num_ghost_cell + i]; //en_den
		u_times_nnvec = 0;
		for( j = 0; j < 3; j++)
		    u_times_nnvec += u[1 + j][num_ghost_cell + i]*nnvec[j];
		for( j = 0; j < 3; j++)
		    u[1 + j][num_ghost_cell - 1 - i] = u[1 + j][num_ghost_cell + i] - 2 * nnvec[j]*u_times_nnvec;
		for( j = 5; j < kmax;j++)    u[j][num_ghost_cell - 1 - i] = u[j][num_ghost_cell + i];
		p[num_ghost_cell - 1 - i] = p[num_ghost_cell + i];
		GAM[num_ghost_cell - 1 - i] = GAM[num_ghost_cell + i];
		c[num_ghost_cell - 1 - i] = c[num_ghost_cell + i];
	    }
	}
	else if(idirs[0] == 2 && upboundary) // sweep in z direction
	{
	    nvec[0] = 0.0;
	    nvec[1] = 0.0;
	    nvec[2] = 1.0;
	    for( j = 0; j < dim; j++)	nnvec[j] = nvec[idirs[j]];
	    for( i = 0; i < num_ghost_cell; i++)
	    {
		u[0][end + 1 + i] = u[0][end - i]; //rho
		u[4][end+ 1 + i] = u[4][end - i]; //en_den
		u_times_nnvec = 0;
		for( j = 0; j < 3; j++)
		    u_times_nnvec += u[1 + j][end - i]*nnvec[j];
		for( j = 0; j < 3; j++)
		    u[1 + j][end + 1 + i] = u[1 + j][end - i] - 2 * nnvec[j]*u_times_nnvec;
		for( j = 5; j < kmax;j++)    u[j][end + 1 + i] = u[j][end - i];
		p[end + 1 + i] = p[end - i];
		GAM[end + 1 + i] = GAM[end - i];
		c[end + 1 + i] = c[end - i];
	    }	    
	}
	else
	    printf("ERROR! in g_update_point_ghost_cell!\n");
}
		
LOCAL	void	Update_pressure_gam_c(int start, int end, float** u, float* p, float* GAM, float* c,Vec_Gas *vst)
{
	float      *mom0, *mom1, *mom2, *rho, *en_den;
	float	   v0,v1,v2,re,e,gam,gm,c2, pressure;
	int i, j;
	rho = u[0];
	mom0 = u[1];
	mom1 = u[2];
	mom2 = u[3];
	en_den = u[4];

	for( i = start; i < end; i++)
	{
	    v0 = mom0[j]/rho[i];
	    v1 = mom1[j]/rho[i];
	    v2 = mom2[j]/rho[i];
	    re = en_den[i] - 0.5*rho[i]*(v0*v0 + v1*v1 + v2*v2);
	    e = re/rho[i];
	    gm = MGamma(vst->state[0]);
	    gam = gm - 1;
	    p[i] = gam*re;
	    GAM[i] = gam;
	    c2 = gm*p[i]/rho[i];
	    c[i] = sqrt(c2);
	}

	return;
}
	

LOCAL	void	compute_weno_flux ( Front* fr, Vec_Gas* vst, int start, int end, int sten_rad,int vsize, int num_comps, int offset, int *idirs,  float** u, float* p, float* GAM, float* c, float** F_mid)
{
	float	   **L, **R;		// Left and Right eigenvector matrices.
	
	float	   **u_tmp;		// Vector of state variables at a given point.
					// 	Used in an intermediate step for the 
					//	calculation of the WENO flux.
	float	   **char_u;		// u_tmp mapped to characteristic field.		
	
	float      **F_old;		// Flux at cell centers at current time step.
	float      **F_tmp;
	float	   **F_plus;		// Positive part of flux used in flux splitting.
	float	   **F_minus;
	float	   **char_F;		// F_tmp mapped to characteristic field ( L . F)
	float	   *Fm_plus;
	float	   *Fm_minus;
	
	float	   *maxeig;		// Maximum eigenvalues ( used in flux splitting (2.5) )
	
	float	   GAMm;		// Gruneisen coefficient at midstate.
	float	   *mom0, *mom1, *mom2, *rho, *en_den;  // State variables at cell centers
	float	   vm0, vm1, vm2, cm, rhom, pm, hm;	// State variables at midstate.
	float	   rhom0[num_comps];				
	
	int	   s;			// Index in characteristic space.


	/* Set pointers */
	rho = u[0];
	mom0 = u[1];	
	mom1 = u[2];	
	mom2 = u[3];	
	en_den = u[4];

	/* Allocate memory for variables above */
	uni_array( &maxeig, 5, sizeof(float));
	uni_array( &Fm_plus, 5, sizeof(float));
	uni_array( &Fm_minus, 5, sizeof(float));
	bi_array( &u_tmp, 6, 5, sizeof(float));
	bi_array( &char_u, 6, 5, sizeof(float));
	bi_array( &char_F, 6, 5, sizeof(float));
	bi_array( &F_plus, 5, 5, sizeof(float));
	bi_array( &F_minus, 5, 5, sizeof(float));
	bi_array( &F_old, end + 3, 5, sizeof(float));
	bi_array( &F_tmp, 6, 5, sizeof(float));
	bi_array( &L, 5, 5, sizeof(float));
	bi_array( &R, 5, 5, sizeof(float));
	
	/* Compute fluxes at cell centers and store in F_old*/
	load_state_fluxes(F_old, start- 3, end + 3, u, p);
	
	/* Compute maximum eigenvalues */
	for ( j = start; j < end; j++ ){
		maxeig[0] = max( maxeig[0], fabs( mom0[j]/rho[j] - c[j] ));
		maxeig[1] = max( maxeig[1], fabs( mom0[j]/rho[j] ));
		maxeig[4] = max( maxeig[4], fabs( mom0[j]/rho[j] + c[j] ));
	}
	maxeig[2] = maxeig[3] = maxeig[1];

	/* Main loop to compute weno fluxes used to update states 
	   for the next time step. */
	for ( j = start; j < end + 1; j++){
		/* Compute midstate variables, Equation (4.1),
		   by averaging two adjacent states */
		rhom = 0.5 * (rho[j-1] + rho[j]);
		pm = 0.5 * ( p[j-1] + p[j] );
		vm0 = 0.5 * ( mom0[j-1] + mom0[j] ) / rhom;
		vm1 = 0.5 * ( mom1[j-1] + mom1[j] ) / rhom;
		vm2 = 0.5 * ( mom2[j-1] + mom2[j] ) / rhom;
		GAMm = 0.5 * ( GAM[j-1] + GAM[j] );
		cm = sqrt( (GAMm+1.0)*pm/rhom );//0.5*( c[j-1] + c[j] );
		hm = 0.5*(vm0*vm0 + vm1*vm1 + vm2*vm2) + (GAMm+1.0)*pm/(GAMm*rhom);
		
		/* Load right eigenvector matrix */
	    	R[0][0] = 1.0;
	    	R[0][1] = 1.0;
	    	R[0][2] = 0.0;
	    	R[0][3] = 0.0;
	    	R[0][4] = 1.0;
		
		R[1][0] = vm0 - cm;
	    	R[1][1] = vm0;
	    	R[1][2] = 0.0;
	    	R[1][3] = 0.0;
	    	R[1][4] = vm0 + cm;
	    	
		R[2][0] = vm1;
	    	R[2][1] = vm1;
	    	R[2][2] = 1.0;
	    	R[2][3] = 0.0;
	    	R[2][4] = vm1;
		
		R[3][0] = vm2;
		R[3][1] = vm2;
		R[3][2] = 0.0;
		R[3][3] = 1.0;
		R[3][4] = vm2;
		
		R[4][0] = hm - vm0*cm;
		R[4][1] = .5*(vm0*vm0 + vm1*vm1 + vm2*vm2);
		R[4][2] = vm1;
		R[4][3] = vm2;
		R[4][4] = hm + vm0*cm;
		
		/* Load left eigenvector matrix */
	    	L[0][0] = hm + cm*(vm0-cm)/GAMm;
	    	L[0][1] = -(vm0 + cm/GAMm);
	    	L[0][2] = -vm1;
	    	L[0][3] = -vm2;
	    	L[0][4] = 1.0;
	    	
		L[1][0] = -2.0*hm + 4*cm*cm/GAMm;
	    	L[1][1] = 2.0*vm0;
	    	L[1][2] = 2.0*vm1;
	    	L[1][3] = 2.0*vm2;
	    	L[1][4] = -2.0;
	    
		L[2][0] = -2*vm1*cm*cm/GAMm;
	    	L[2][1] = 0.0;
	    	L[2][2] = 2*cm*cm/GAMm;
	    	L[2][3] = 0.0;
	    	L[2][4] = 0.0;
		
		L[3][0] = -2.0*vm2*cm*cm/GAMm;
		L[3][1] = 0.0;
		L[3][2] = 0.0;
		L[3][3] = 2.0*cm*cm/GAMm;
		L[3][4] = 0.0;
		
		L[4][0] = hm - cm*(vm0+cm)/GAMm;
		L[4][1] = -vm0 + cm/GAMm;
		L[4][2] = -vm1;
		L[4][3] = -vm2;
		L[4][4] = 1.0;	

		for ( i = 0; i < 5; i++ )
			for ( k = 0; k < 5; k++ )
				L[i][k] *= GAMm / (2.0*cm*cm) ;

		/* Set mid state mass fractions */
		for ( i = 0; i < num_comps; i++)
			rhom0[i] = (u[5+i][j-1] + u[5+i][j])/(u[0][j-1] + u[0][j]);

		/* Flux Splitting */
		for ( i = 0; i < 6; i++){
			F_tmp[i] = F_old[j-3+i];
			u_tmp[i][0] = rho[j-3+i];
			u_tmp[i][1] = mom0[j-3+i];
			u_tmp[i][2] = mom1[j-3+i];
			u_tmp[i][3] = mom2[j-3+i];
			u_tmp[i][4] = en_den[j-3+i];
		}			
		
		/* Map to characteristic field */
		for ( i = 0; i < 6; i++){
			matrix_times_vector( char_F[i], L, F_tmp[i], 5 );
			matrix_times_vector( char_u[i], L, u_tmp[i], 5 );
		}
		
		for ( i = 0; i < 5; i++){ // loop over stencil
			for ( k = 0; k < 5; k++){
				F_plus[i][k] = 0.5 * ( char_F[k][i] + maxeig[i]*char_u[k][i] );
				F_minus[i][k] = 0.5 * ( char_F[5-k][i] - maxeig[i]*char_u[5-k][i] );
			}
		}

		compute_weno_weights( F_plus, Fm_plus);
		compute_weno_weights( F_minus, Fm_minus);

		for ( i = 0; i < 5; i++)
			F_tmp[0][i] = Fm_plus[i] + Fm_minus[i];

		/* Map back to component field */
		matrix_times_vector( F_mid[j], R, F_tmp[0], 5);
		

		if (isnan(F_mid[j][0]) || isnan(F_mid[j][1]) || isnan(F_mid[j][2]) || isnan(F_mid[j][3]) || isnan(F_mid[j][4])){
			printf("In compute_weno_flux()\n");
			printf("j = %d\n",j);
			printf("vsize = %d\n",vsize);
			printf("start = %d\n",start);
			printf("end = %d\n",end);
			printf("rhom = %e\n",rhom);
			printf("(vm0, vm1, vm2) = (%e, %e, %e)\n",vm0,vm1,vm2);
			printf("pm = %e\n",pm);
			printf("cm = %e, c0 = %e, c1 = %e\n",cm,c[j-1],c[j]);
			printf("hm = %e\n",hm);
			printf("GAMm = %e\n",GAMm);
			printf("\n");
		}

		/* Compute fluxes for species */
		for ( k = 5; k < 5 + num_comps; k++)
			F_mid[j][k] = rhom0[k-5] * F_mid[j][0];
		
	}

	/* Free memory */
	free( maxeig );
	free( Fm_plus );
	free( Fm_minus );
	free( u_tmp );
	free( char_u );
	free( char_F );
	free( F_plus );
	free( F_minus );
	free( F_old );
	free( F_tmp );
	free( L );
	free( R );

}
LOCAL	void	update_pressure_soundspeed_gamma(Front* fr, Vec_Gas* vst, int indx, int start, int end, int kmax, int offset, float** u, float* p, float* c,float* GAM)
{
	Locstate tmpst = NULL;
	
	if ( tmpst == NULL )
		g_alloc_state(&tmpst, fr->sizest);

	set_type_of_state(tmpst,GAS_STATE);

	int i,j;
	for ( i = start; i < end; i++){
		Set_params(tmpst,vst->state[indx]);
		
		Dens(tmpst) = u[0][i];
		Mom(tmpst)[0] = u[1][i];
		Mom(tmpst)[1] = u[2][i];
		Mom(tmpst)[2] = u[3][i];
		Energy(tmpst) = u[4][i];
		
		for ( j = 5; j < kmax; j++ ){
			pdens(tmpst)[j-5] = u[j][i];
		}

		p[i] = pressure(tmpst);
		GAM[i] = gruneisen_gamma(tmpst);
		c[i] = sqrt( sound_speed_squared(tmpst) );
		if ( isnan(c[i]) ){ 
			printf("NAN sound speed in updates step!\n");
			printf("GAM = %e\n",GAM[i]);
			printf("rho = %e, rho-1 = %e, rho+1 =%e, rho =%e\n",Dens(tmpst),u[0][i-1],u[0][i+1],u[0][i]);
			printf("ie = %e\n",internal_energy(tmpst));
			for ( j = 0; j < kmax; j++)
				printf("u[%d][%d] = %e\n",j,i,u[j][i]);
		}
	}
	free ( tmpst );
}
LOCAL	void	compute_weno_weights ( float** LF, float* f_new)
{
	int	i, j, k, s;
	float	C[3] = {1.0/10.0, 6.0/10.0, 3.0/10.0}; 			// Table II
	float	a[3][5] = { {1.0/3.0, -7.0/6.0, 11.0/6.0, 0.0, 0.0},
			    {0.0, -1.0/6.0, 5.0/6.0, 1.0/3.0, 0.0},
			    {0.0, 0.0, 1.0/3.0, 5.0/6.0, -1.0/6.0} }; 	// Table I
	float	eps = 1e-6;
	float	*lsf;
	float	alpha[3];
	float	sum_alpha;
	float	omega[3];
	float	IS[3]; 
	float	fs; 
	//float	q[3];
	float	q[5];

	for ( s = 0; s < 5; s++){
		lsf = LF[s];

		
		compute_IS( IS, lsf );
		
		sum_alpha = 0.0;
		for ( k = 0; k < 3; k++){
			alpha[k] = C[k] / pow(eps + IS[k],2); // Equation (2.16)
			sum_alpha += alpha[k];
		}
		
		for ( k = 0; k < 3; k++)
			omega[k] = alpha[k] / sum_alpha;  // Equation (4.3)
		
		fs = 0.0;
		for ( i = 0; i < 5; i++ ){
			q[i] = 0.0;
			for ( j = 0; j < 3; j++ ){
				q[i] += omega[j] * a[j][i];
			}
		}
		for ( i =0; i < 5; i++){
			fs += q[i] * lsf[i];
		}
		/*
		if ( isnan(fs) ){
			printf("VINAY DEBUG\n");
			for ( i = 0 ; i < 5; i++)
				printf("lsf[%d] = %f\n",i,lsf[i]);
			for ( i = 0 ; i < 5; i++)
				printf("q[%d] = %f\n",i,q[i]);
			for ( i = 0; i < 3; i++)
				printf("omega[%d] = %f\n",i,omega[i]);
			clean_up(-1);
		}*/
		f_new[s] = fs;
	}
}

LOCAL	void	compute_IS ( float* IS, float* lsf ) 
{
	// Equations (3.2) - (3.4)
	IS[0] = ( ( 13.0/12.0) * pow( lsf[0] - 2*lsf[1] + lsf[2], 2 ))
		   + ((1.0/4.0) * pow(lsf[0] - 4*lsf[1] + 3*lsf[2], 2 ));

	IS[1] = ( ( 13.0/12.0) * pow( lsf[1] - 2*lsf[2] + lsf[3], 2 ))
		   + ((1.0/4.0) * pow( lsf[1] - lsf[3], 2 ));

	IS[2] = ( ( 13.0/12.0) * pow( lsf[2] - 2*lsf[3] + lsf[4], 2 ))
		   + ((1.0/4.0) * pow( 3*lsf[2] - 4*lsf[3] + lsf[4], 2 ));
}
	
/*  Computes flux vectors from cell centered values */
LOCAL	void	load_state_fluxes( float** F, int start, int end, float** u, float* p)
{
	int	i, j;
	float*	rho = u[0];
	float*	mom0 = u[1];
	float*	mom1 = u[2];
	float*	mom2 = u[3];
	float*	E = u[4];
	float	v0, v1, v2;

	/* Each weighted flux is calculated using the cell centered fluxes 
	    of a five point stencil */
	for ( j = start; j < end; j++){
		v0 = mom0[j]/rho[j];
		v1 = mom1[j]/rho[j];
		v2 = mom2[j]/rho[j];

	    	/* mass */
	        F[j][0] = rho[j]*v0;
	    	/* Sweep component of momentum */
	        F[j][1] = F[j][0]*v0 + p[j];
	    	/* Off sweep component of momentum */
	        F[j][2] = F[j][0]*v1;
	    	/* Off sweep component of momentum */
	        F[j][3] = F[j][0]*v2;
	    	/* energy */
	        F[j][4] = v0 * (E[j]+p[j] );

		if ( isnan(F[j][0]) || isnan(F[j][1]) || isnan(F[j][2]) || isnan(F[j][3]) || isnan(F[j][4])){
		printf("rho[%d] = %e\n",j,rho[j]);
		printf("E[%d] = %e\n",j,E[j]);
		printf("p[%d] = %e\n",j,p[j]);
		printf("(v0,v1,v2) = (%e, %e, %e)\n",v0,v1,v2);
		printf("\n");
	}
	}
	
}
LOCAL	void	inner_prod ( float* c, float* a, float* b, int dim)
{
	int i;
	*c = 0.0;
	for ( i = 0; i < dim; i++)
		*c += a[i]*b[i];	
}
LOCAL	void	matrix_times_vector ( float* b, float** A, float* x, int dim )
{
	int i, j;
	
	for ( i = 0; i < dim; i++){
		b[i] = 0.0;
		for ( j = 0; j < dim; j++){
			b[i] += A[i][j]*x[j];
		}
	}
}
EXPORT void oned_MUSCL(
	int		swp_num,
	int		*iperm,
	int		*icoords,
	Wave		*wave,
	Wave            *newwave,
	Front		*fr,
	Front           *newfr,
	Stencil         *sten,
	Tan_stencil	*tsten,
	int		offset,
	int		vsize,
	Vec_Gas		*vst,
	Vec_Src		*src,
	float		dt,
	float		dn,
	int		dim)
{
	Vec_Muscl  *vmuscl;
	float	   dtdni = dt/dn;
	float	   **ucon, **F, **source;
	float	   *uconk, *Fk, *sourcek;
	int        start, end;
	int        i, j, k, kmax;
	int        sten_rad;
//#define	DEBUG_MUSCL
#if defined(DEBUG_MUSCL)
	const char	**debug_strings;
#endif /* defined(DEBUG_MUSCL) */

	debug_print("MUSCL","Entered oned_MUSCL()\n");

#if defined(DEBUG_MUSCL)
	debug_strings = set_muscl_debugging(swp_num,iperm,
	        			    icoords,tsten,fr,vsize);
#endif /* defined(DEBUG_MUSCL) */

	vmuscl = load_state_data(swp_num,iperm,fr,wave,sten,tsten,
	        		 vst,src,offset,vsize,dim,dt,dn);	//g_load_muscl_state_data

	sten_rad = vmuscl->sten_rad;
	
	/* compute the eigenvalues, eigenuni_arrays and max wavespeed*/
	compute_eigens(0,vsize,vmuscl);

	start = sten_rad-1;
	end = vsize-start;
	/* compute the coefficients for Lapidus and slope limiting viscosity */
	compute_art_visc_coefs(start,end,vmuscl);

	/* compute the linear reconstruction of the state variables */
	reconstructor(start,end,vmuscl);

	kmax = dim+2;
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* Compute differences for mass fraction */
            Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
            int    num_comps;
            if((num_comps = params->n_comps) != 1)
                kmax += num_comps;
        }

	bool *Flag;
	uni_array(&Flag, kmax + 1, sizeof(bool));
	for (i = 0; i < kmax; ++i)
	    Flag[i] = NO;
	if((swp_num == 0) || (swp_num == 1) || (swp_num == 2))
	    for (j = start; j < end; ++j)
	    {
		for (i = 0; i < kmax; ++i)
		{
		    //if (vmuscl->du[i][j] != 0)
		    //	Flag[i] = YES;
			Flag[i] = NO;
		    Lflag(vmuscl->vst->state[j])[i] = Flag[i];
		}
	    }
	free(Flag);

	/* Evolve for half time step.  The returned data uL and uR
	 * are these evolved states on the two sides of each mesh edge. 
	 * They will be put in the Riemann solver. */
	start = sten_rad;
	end = vsize - sten_rad + 1;
	half_step(start,end,dt,dn,vmuscl);

	/* Solve the Riemann problem at each mesh edge with the states 
	 * uL and uR on the two sides.  The returned data uM is the solution
	 * state at the middle. */

#if defined(DEBUG_MUSCL_TIME)
	start_clock("Riemann_solver");
#endif /* defined(DEBUG_MUSCL_TIME) */

	rsolver(start,end,vmuscl->uL,&vmuscl->VLst,vmuscl->uR,&vmuscl->VRst,
	                  vmuscl->uM,&vmuscl->VMst,&vmuscl->Flux,vmuscl);
	if (vmuscl->avisc)
	    add_art_visc2(start,end,vmuscl);

#if defined(DEBUG_MUSCL)
	if (debugging("lriem"))
	{
	    int i;
	    float **uM = vmuscl->uM;
	    float *pM = vmuscl->pM;
	    g_printout_vec_data("Mid state on the edge, "
	        		"obtained by Riemann solver",
	        	        uM[0],uM[1],uM+2,dim,start,end,
	        		"muncons");
	    (void) printf("Mid state pressure\n"); 
	    (void) printf("%-4s %-14s\n","n","pM"); 
	    for (i = start; i < end; ++i)
	    	(void) printf("%-4d %-14g\n",i,pM[i]);
	    (void) printf("\n");
	}
	if (debugging("vflux")) 
	{
	    float **F = vmuscl->Flux.F;
	    g_printout_vec_data("Here are the flux uni_arrays",F[0],F[1],F+2,
	        		dim,start,end,"cons");
        }
#endif /* defined(DEBUG_MUSCL) */

#if defined(DEBUG_MUSCL_TIME)
	stop_clock("Riemann_solver");
#endif /* defined(DEBUG_MUSCL_TIME) */

	/* Set conservative source terms */
	sten_rad = vmuscl->sten_rad;
	start = sten_rad;
	end = vsize - sten_rad;
	cons_src(start,end,swp_num,iperm,tsten,vmuscl);

	/* compute the cell average of the approximate solution for the
	   next time step */
	kmax =		dim+2;
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* Compute differences for mass fraction */
            Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
            int    num_comps;
            if((num_comps = params->n_comps) != 1)
                kmax += num_comps;
        }
	ucon =		vmuscl->ucon;
	F =		vmuscl->Flux.F;
	source =	vmuscl->source;

        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* Convert mass fraction to partial density */
            Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
            int    num_comps;
            if((num_comps = params->n_comps) != 1)
            {
                for(k = 0; k < num_comps; k++)
                {
                    uconk = ucon[k+dim+2];
                    for (j = start; j < end; ++j)
                        uconk[j] *= ucon[0][j];
                }
            }
        }

	/*Y[0] to Y[num_comps-1] contain the mass fractions of the components.
	  Y[num_comps] contains energy.
	  Y[num_comps+1] and Y[num_comps+2] contain temp, and density (not part of ODE)*/
	
        double *Y, *yp;
        Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
        int    num_comps;
        num_comps = params->n_comps;
        int istate;

	/*pointer to function that computes right hand side of dY/dt = rhs(t,Y)*/
	/*Solves step 1 of operator splitting method*/
	//Only the average density, energy, momentum in normal cell are updated here. The physical variable in cutting cells are updated in function Update_merge_cell.
		
	if((iperm != NULL && swp_num == -1) && (iperm == NULL && swp_num != -1)) printf("error in iperm.");
	int   idirs[3];
	if( iperm!= NULL && swp_num != -1)
	{
	    for (i = 0; i < dim; ++i) 
		idirs[i] = iperm[(i+swp_num)%dim];
	}		
	else
	{
	    for (i = 0; i < dim; ++i) idirs[i] = i;
	}
	
	for (j = start; j < end; ++j)
	{
	    if( Merge_cell(vmuscl->vst->state[j])  &&fr->turb_boundary_layer == YES && iperm != NULL && swp_num != -1 )
	    {
		for (k = 0; k < 2; ++k)
		{
		    uconk = ucon[k];
	    	    Fk = F[k];
		    i = 2*iperm[swp_num%3]; /*sweep direction.*/
		    Flux(vmuscl->vst->state[j])[i][k] = Fk[j]; 
                    Flux(vmuscl->vst->state[j])[i+1][k] = Fk[j+1];	
		    flux_flag(vmuscl->vst->state[j]) = 1;
		}
		for (k = 2; k < 5; ++k)
		{
		    uconk = ucon[k];
	            Fk = F[k];
		    i = 2*iperm[swp_num%3];/*sweep direction.*/
		    Flux(vmuscl->vst->state[j])[i][idirs[k-2]+2]   = Fk[j]; 
                    Flux(vmuscl->vst->state[j])[i+1][idirs[k-2]+2] = Fk[j+1]; 		
		}

		for (k = 5; k < kmax; ++k)
		{
	  	    uconk = ucon[k];
	   	    Fk = F[k];
		    i = 2*iperm[swp_num%3];
		    Flux(vmuscl->vst->state[j])[i][k] = Fk[j]; //  j
                    Flux(vmuscl->vst->state[j])[i+1][k] = Fk[j+1]; //j+1		
	   	}
	    }
            else
	    {
		for (k = 0; k < kmax; ++k)
		{ 
		uconk = ucon[k];
		Fk = F[k];
		uconk[j] += dtdni*(Fk[j] - Fk[j+1]);		
		}
	    }

	}


        /*Solves step 2 of operator splitting method*/
	if(NO)
        {
	    float *coords;
	    int ix, iy, iz;
	    ix = icoords[0];
	    iy = icoords[1];
	    iz = icoords[2];
	    for(j = start; j < end; ++j)
	    {
		float T = temperature(vmuscl->vst->state[j]);
		if (T>=200 && T <=3500)
		{
                	uni_array(&Y,num_comps+5,sizeof(double));

			Y[num_comps+1] = ucon[1][j];  //set initial energy for ODE
			for(k = 0; k < params->n_comps; k++)
				Y[1+k] = ucon[k+dim+2][j];
			Y[num_comps+2] = ucon[0][j];  //set initial density for ODE
			Y[num_comps+3] = T;  //set initial temperature for ODE
			Y[num_comps+4] = pressure(vmuscl->vst->state[j]);  //set initial pressure for ODE
	
	                istate = 1;

/*			float *coords;
			coords = Rect_coords(icoords, wave);
			printf("At coords: (%f, %f, %f)\n", coords[0], coords[1], coords[2]);*/

                	if (!ode_solver(ChemHeatSources,Y,0.0,dt,1,1.e-8,
                	                 10*MACH_EPS,
                	                 1.e-4,&istate,NULL))
                	{
			    printf("Before ODE crashes, energy = %E\n", ucon[1][j]);
			    for(k=0;k<10;k++)
				printf("partial density %d = %E\n",k,pdens(vmuscl->vst->state[j])[k]);
			    printf("T = %f\n",temperature(vmuscl->vst->state[j]));
	 		    printf("P = %f\n",pressure(vmuscl->vst->state[j]));
			    printf("After ODE crashes, energy = %E\n", Y[NUM_OF_SPECIES+1]);
			    printf("Temperature: %E\n", Y[num_comps+3]);
			    printf("Pressure: %E\n", Y[num_comps+4]);
			    printf("IE = %f\n",internal_energy(vmuscl->vst->state[j]));
		            printf("iCoords: (%d, %d, %d)\n", icoords[0], icoords[1], icoords[2]);
		    	    printf("Coords: (%f, %f, %f)\n", coords[0], coords[1], coords[2]);
		    	    printf("j = %d\n",j);
		    	    printf("start = %d\n",start);
		    	    printf("end = %d\n",end);
                    
		   	    screen("ERROR in oned_MUSCL(), "
                           	"ode_solver() failed for ChemHeatSources\n");
                            clean_up(ERROR);
                	}

                	ucon[1][j] = Y[num_comps+1];   
                	for(k = 0; k < params->n_comps; k++)
                	    ucon[k+dim+2][j] = Y[k+1];

                	free(Y);
		}
            } 
        } 

	/*if (source != NULL)

	{
	    for (k = 0; k < kmax; ++k)
	    {
	    	uconk = ucon[k];
	    	Fk = F[k];
	    	sourcek = source[k];
	    	for (j = start; j < end; ++j)
		{
	    	    uconk[j] += dtdni*(Fk[j] - Fk[j+1]) + dt*sourcek[j];
		}
	    }
	}
	else
	{
	    for (k = 0; k < kmax; ++k)
	    {
	    	uconk = ucon[k];
	    	Fk = F[k];
	    	for (j = start; j < end; ++j)
		{
	    	    uconk[j] += dtdni*(Fk[j] - Fk[j+1]);
		}
	    }
	}*/

#if defined(DEBUG_MUSCL)
	if (debugging("vgdnv"))
	{
	    g_printout_vec_data("Here are the updated cell averages:",
	        	        ucon[0],ucon[1],ucon+2,dim,start,end,"cons");
	}
	if (debugging("vgdnv_ie"))
	{
	    g_print_internal_energy("Internal energy of updated cell averages:",
	        	            ucon,vmuscl,start,end);
	}
#endif /* defined(DEBUG_MUSCL) */

	debug_print("MUSCL","Left oned_MUSCL()\n");
#if defined(DEBUG_MUSCL)
	reset_muscl_debugging(debug_strings);
#endif /* defined(DEBUG_MUSCL) */
}		/*end oned_MUSCL*/

#if defined(__cplusplus)
extern "C" {
#endif /* defined(__cplusplus) */

LOCAL void ChemHeatSources(
        int             *neq,
        float           *t,
        float           *Y,
        float           *yp)
{
        //printf("In the ChemHeatSources\n");  //you can print out this for the number of subcycling.

        double temp, dens, pres;
	double TBConc;
        int i, j, k;
        float YY[NUM_OF_SPECIES];
     
        *neq = NUM_OF_SPECIES + 2;   //num_comp+1
        dens = Y[NUM_OF_SPECIES + 2];
        temp = Y[NUM_OF_SPECIES + 3];
	pres = Y[NUM_OF_SPECIES + 4];
        for (i = 0; i < NUM_OF_SPECIES; ++i)
            YY[i] = Y[i+1]/dens;       
        
	ReadThermoData();

        for(i = 0 ; i < NUM_OF_REACTIONS; ++i)
	{
	    	TBConc = 0.0;
	    	if (reaction[i].flag_tb == 1)
		{
		    if (reaction[i].tbIndex == 0)
		    {
			k = reaction[i].ThirdBody[0].index;
			TBConc = YY[k] / species[k].MolarMass;
		    }
		    else
		    {
			//TMP
			//printf("Entering third body M\n");
			for (j = 0; j < reaction[i].nofcomp; j++)
			{
			    k = reaction[i].ThirdBody[j].index;
			    TBConc += reaction[i].ThirdBody[j].coeff * YY[k] / species[k].MolarMass;
			}
			//TMP_XY hardcode
			switch(reaction[i].tbIndex)
			{
			    case 1: TBConc += YY[1] / species[1].MolarMass + YY[2] / species[2].MolarMass
				    + YY[3] / species[3].MolarMass + YY[6] / species[6].MolarMass + YY[8] / species[8].MolarMass; break;
			    case 2: TBConc += YY[0] / species[0].MolarMass + YY[1] / species[1].MolarMass
				    + YY[2] / species[2].MolarMass + YY[3] / species[3].MolarMass + YY[6] / species[6].MolarMass; break;
		            case 3: TBConc += YY[1] / species[1].MolarMass + YY[2] / species[2].MolarMass
				    + YY[3] / species[3].MolarMass + YY[6] / species[6].MolarMass + YY[7] / species[7].MolarMass; break;
			    case 4: TBConc += YY[1] / species[1].MolarMass + YY[2] / species[2].MolarMass
				    + YY[3] / species[3].MolarMass + YY[6]/ species[6].MolarMass; break;
			    case 5: TBConc += YY[0] / species[0].MolarMass + YY[1] / species[1].MolarMass + YY[2] / species[2].MolarMass
				    + YY[3] / species[3].MolarMass + YY[6] / species[6].MolarMass + YY[8] / species[8].MolarMass; break;
			    case 6: TBConc += YY[0] / species[0].MolarMass + YY[1] / species[1].MolarMass + YY[2] / species[2].MolarMass
				    + YY[3] / species[3].MolarMass + YY[6] / species[6].MolarMass + YY[8] / species[8].MolarMass; break;
			}

		    }
		    TBConc *= dens;
		}
		else
		{
		    TBConc = 1;
		}
		if (reaction[i].flag_fm == 0)
		{
		    ComputeRateCoefficients( &reaction[i].RateCoeff, reaction[i].PreExp, reaction[i].TempExp, reaction[i].ActEnergy / RGAS, temp );
		}
		else
		{
		    double Conc;
		    double k_knot, k_infty, Pr;
		    double c, n, fcent = 0.0, d = 0.14, f;
		    if (reaction[i].flag_tb == 1)
		    {
			Conc = TBConc;
		    }
		    else
		    {
			Conc = pres / (RGAS * temp);
		    }
		    ComputeRateCoefficients( &k_knot, reaction[i].PreExp, reaction[i].TempExp, reaction[i].ActEnergy / RGAS, temp);
		    ComputeRateCoefficients( &k_infty, reaction[i].PreExp_i, reaction[i].TempExp_i, reaction[i].ActEnergy_i / RGAS, temp);
		    Pr = k_knot * Conc / ((k_infty > 1.e-60) ? k_infty : 1.e-60); 
		    
		    if( reaction[i].fca || reaction[i].fcb || reaction[i].fcc )
		    {
			if( reaction[i].fca )
			{
			    fcent += reaction[i].fca * exp( -temp / reaction[i].fcta );
			}
			if( reaction[i].fcb)
			{
			    fcent += reaction[i].fcb * exp( -temp / reaction[i].fctb );
			}
			if( reaction[i].fcc)
			{
			    fcent += reaction[i].fcc * exp( -reaction[i].fctc / temp );
			}
		    }
		    n = 0.75 - 1.27 * myLog10( fcent );
		    c = -0.4 - 0.67 * myLog10( fcent );
		    f = (myLog10( Pr ) + c) / (n - d * (myLog10( Pr ) + c));
		    f = pow( fcent, 1.0 / (f * f + 1.0));
		    reaction[i].RateCoeff = k_infty * Pr * f / (1.0 + Pr);
		    //TMP
		    //printf("k_knot = %e\tk_infty = %e\tPr = %e\nfcent = %e\tn = %e\tc = %e\tf = %e\n", k_knot, k_infty, Pr, fcent, n, c, f);
		    //printf("TROE: Rate coefficient for reaction %d is: %e.\n", i, reaction[i].RateCoeff);
		}
                ComputeBackwardRateCoefficients( temp,i );
		//TMP
		//printf("Backward rate coefficient for reaction %d is: %e.\n", i, reaction[i].BackwardRateCoeff);
		if( reaction[i].PreExp_i != 0 ) TBConc = 1;
                ComputeReactionRate( dens ,YY, i, TBConc );
		//TMP_XY
//		if( i == 24 ) reaction[i].ReactionRate *= 1.e-7;
		if ( temp <= 295 )
		    reaction[24].ReactionRate = 0;
		//TMP
		//printf("Reaction rate for reaction %d = %e.\n", i, reaction[i].ReactionRate);
        }
	yp[0] = 0.0;
        for(i = 0; i < NUM_OF_SPECIES; ++i){
                species[i].ProductionRate = 0;
                ComputeProductionRate(i);
		yp[i+1] = species[i].ProductionRate;
		//TMP
		//printf("Production rate for species %d = %e.\n", i, species[i].ProductionRate);
        }
	//TMP_1109
/*///10	if( flag == 1)
	{
//09	printf("In ODE:\n");
	for( i = 0; i < species[1].NofUsedReactions; i++ )
	{
	    int ii = species[1].UsedReactions->vec[i];
	    if( reaction[ii].ReactionRate >1.e-10 )
	    {
		printf("Reaction %d: %e\n", ii, reaction[ii].ReactionRate);
	    }
	}
	}
	if ( flag == 1 )
	{
	    printf("yp pdens[H] = %e\t", yp[2]);
	    printf("yp dens[H2O] = %e\n", yp[6]);
///	printf("Forward: %e\tBackward: %e\n", reaction[15].RateCoeff, reaction[15].BackwardRateCoeff);
///	printf("Reaction 15: %e\n", reaction[15].ReactionRate);
	}*///10
        //Heat release term
	yp[NUM_OF_SPECIES + 1] = 0.0; 
        for(i = 0; i < NUM_OF_REACTIONS; ++i){
                ComputeHeatRelease(i,temp);
                yp[NUM_OF_SPECIES + 1] += reaction[i].HeatRelease;
		//TMP
		//printf("Heat release for reaction %d = %e.\n", i, reaction[i].HeatRelease);
		//yp[10] += (1.0/3.0)*reaction[i].HeatRelease;
        }
	//TMP_mark
/*	if (tag == 1)
	{
	    printf("Heat release: %E\n", yp[NUM_OF_SPECIES + 1]);
	    tag = 0;
	}*/
        return;
}

#if defined(__cplusplus)
}
#endif /* defined(__cplusplus) */


EXPORT	void	set_muscl_default_opts(
	Muscl_Opts  *mopts,
	Hyp_method  *method,
	int         dim)
{
	MUSCL_PromptType_Reconstructor *Sintrp;
	MUSCL_PromptType_Rsolver       *Rsolver;

	set_muscl_default_hooks(mopts,method,dim);

	uni_array(&Sintrp,5,sizeof(MUSCL_PromptType_Reconstructor));
	mopts->Sintrp = Sintrp;
	Sintrp[0].prompt = "Reconstruct density, energy, velocity";
	Sintrp[0].select = "d";
	Sintrp[0].reconstructor = g_state_linear_reconstructor;
	Sintrp[0].half_step = g_half_step;
	Sintrp[0].strong_wave_half_step = g_strong_wave_half_step;
	Sintrp[1].prompt = "Reconstruct eigen coordinates";
	Sintrp[1].select = "e";
	Sintrp[1].reconstructor = g_eigenstate_linear_reconstructor;
	Sintrp[1].half_step = g_half_step;
	Sintrp[1].strong_wave_half_step = g_strong_wave_half_step;
	Sintrp[2].prompt = "Bell-Colella-Trangenstein reconstruction";
	Sintrp[2].select = "b";
	Sintrp[2].reconstructor = g_bct_linear_reconstructor;
	Sintrp[2].half_step = g_half_step;
	Sintrp[2].strong_wave_half_step = g_strong_wave_half_step;
	Sintrp[3].prompt = "First order Godunov reconstruction (zero slopes)";
	Sintrp[3].select = "f";
	Sintrp[3].reconstructor = g_first_order_godunov_reconstructor;
	Sintrp[3].half_step = g_first_order_godunov_half_step;
	Sintrp[3].strong_wave_half_step = NULL;

	uni_array(&Rsolver,6,sizeof(MUSCL_PromptType_Rsolver));
	mopts->Rsolver = Rsolver;
	Rsolver[0].prompt = "Exact Riemann solver";
	Rsolver[0].select = "e";
	Rsolver[0].rsolver = g_muscl_exact_rsolver;
	Rsolver[0].rmidstate = g_exact_Riemann_midstate;
	Rsolver[1].prompt = "Linear approximate Riemann solver";
	Rsolver[1].select = "l";
	Rsolver[1].rsolver = g_muscl_lin_approx_rsolver;
	Rsolver[1].rmidstate = NULL;/*TODO*/
	Rsolver[2].prompt = "Colella-Glaz's approximate Riemann solver";
	Rsolver[2].select = "c";
	Rsolver[2].rsolver = cg_rsolve;
	Rsolver[2].rmidstate = g_cg_Riemann_midstate;
	Rsolver[3].prompt = "Linear US/UP fit (Dukowicz)";
	Rsolver[3].select = "d";
	Rsolver[3].rsolver = g_linear_us_up_rsolver;
	Rsolver[3].rmidstate = g_linear_us_up_Riemann_midstate;
	Rsolver[4].prompt = "Gamma Law fit";
	Rsolver[4].select = "g";
	Rsolver[4].rsolver = g_gamma_law_fit_rsolver;
	Rsolver[4].rmidstate = g_gamma_law_fit_Riemann_midstate;
}		/*end set_muscl_default_opts*/

EXPORT	void set_muscl_default_hooks(
	Muscl_Opts *mopts,
	Hyp_method  *method,
	int         dim)
{
	mopts->kmax = dim + 2;
	mopts->nfloats = dim+2;
	mopts->_npt_solver = point_FD;
	mopts->_one_side_npt_tang_solver = oblique_FD;
	mopts->_npt_tang_solver = g_two_side_npt_tang_solver;
	mopts->_print_internal_energy = g_print_internal_energy;
	mopts->_npts_tan_sten = method->npts_sten;
	mopts->_compute_art_visc_coefs = g_compute_art_visc_coefs;
	mopts->_load_state_vectors = g_load_state_vectors;
	mopts->_alloc_phys_vecs = g_muscl_alloc_phys_vecs;
	mopts->_free_phys_vecs = g_muscl_free_phys_vecs;
	mopts->_assign_wave_state_vectors = g_assign_wave_state_vectors;

	mopts->_reconstructor = g_state_linear_reconstructor;
	mopts->_load_state_data = g_load_muscl_state_data;
	mopts->_compute_eigens = g_compute_eigens;
	mopts->_flux_vectors = g_muscl_flux_vectors;
	mopts->_rsolver = g_muscl_exact_rsolver;
	mopts->_rmidstate = g_exact_Riemann_midstate;
	mopts->worksp_len = 5;
	mopts->_half_step = g_half_step;
	mopts->_strong_wave_half_step = g_strong_wave_half_step;
	mopts->_cons_src = g_cons_src;
	mopts->_add_art_visc1 = g_add_art_visc1;
	mopts->_add_art_visc2 = g_add_art_visc2;
	mopts->monotone_reconstruction = NO;
	mopts->link_reconstructions = NO;

	mopts->_Cg_params.pv_iterations = 4;
	mopts->_Cg_params.min_p_jump = 1e-06;          /*TOLERANCE*/
	mopts->_Cg_params.min_v_jump = 1e-06;          /*TOLERANCE*/
	mopts->_Cg_params.sw_tol = 100.0;              /*TOLERANCE*/
	mopts->_Cg_params.min_mass_flux = 1.0e-12;     /*TOLERANCE*/
	mopts->_Cg_params.sqr_min_mass_flux = 1.0e-24; /*TOLERANCE*/
}		/*end set_muscl_default_hooks*/
	

/*ARGSUSED*/
LOCAL	Vec_Muscl *g_load_muscl_state_data(
	Muscl_Opts	*mopts,
	int		swp_num,
	int		*iperm,
	Front		*fr,
	Wave		*wave,
	Stencil         *sten,
	Tan_stencil	*tsten,
	Vec_Gas		*vst,
	Vec_Src		*src,
	int		offset,
	int		vsize,
	int		dim,
	float		dt,
	float		dn)
{
	static Vec_Muscl *vmuscl = NULL;
	float		 **u, **ucon, **source;
	const float* const *Q = vst->Q;
	bool		 is_src;
	int		 i;

	if (vmuscl == NULL || vsize > vmuscl->max_vsize || dim > vmuscl->dim)
	{
	    vmuscl = g_muscl_free_wk_space(vmuscl);
	    vmuscl = g_muscl_alloc_wk_space(mopts,vsize,dim,vmuscl);
	}
	vmuscl->sizest = fr->sizest;
	vmuscl->idir = (iperm != NULL) ? iperm[swp_num] : dim;
	if(iperm!=NULL){
	for(i = 0; i< dim; ++i)
	   vmuscl->iperm[i] = iperm[i];
	}
	vmuscl->swp_num= swp_num;
	vmuscl->Q = Q;
	vmuscl->dt = dt;
	vmuscl->dn = dn;
	vmuscl->front = fr;
	vmuscl->wave = wave;
        vmuscl->sten = sten;
        vmuscl->tsten = tsten;
        vmuscl->coords = vst->coords;
        for (i = 0; i < dim; ++i)
        {
            vmuscl->GL[i] = fr->rect_grid->GL[i];
            vmuscl->GU[i] = fr->rect_grid->GU[i];
            vmuscl->dh[i] = fr->rect_grid->h[i];
        }
        vmuscl->turb_boundary_layer = fr->turb_boundary_layer;
	if (wave != NULL)
	    vmuscl->sten_rad = wave->npts_sten/2;
	else if (tsten != NULL)
	    vmuscl->sten_rad = tsten->npts/2;
	else
	{
	    vmuscl->sten_rad = 2;
	    screen("ERROR in g_load_muscl_state_data(), can't determine "
	           "stencil size\n");
	    clean_up(ERROR);
	}
	vmuscl->vst = vst;
	vmuscl->src = src;
	vmuscl->offset = offset;	vmuscl->vsize = vsize;

	g_load_VGas_state_vectors(offset,vsize,vst,dim);

	clear_Vec_Gas_set_flags(&vmuscl->VMst);
	clear_Vec_Gas_set_flags(&vmuscl->VLst);
	clear_Vec_Gas_set_flags(&vmuscl->VRst);
	vmuscl->VMst.Q = vmuscl->VLst.Q = vmuscl->VRst.Q = Q;

	/* Set params jump points for vec gas states */
	copy_vec_state_params(&vmuscl->VMst,vst,offset,vsize);
	copy_vec_state_params(&vmuscl->VRst,vst,offset,vsize);
	copy_vec_state_params(&vmuscl->VLst,vst,offset-1,vsize);

	/*put the state variables and source terms in a form more
	  easily usable by the subroutines*/

	u = vmuscl->u;
	ucon = vmuscl->ucon;
	u[vmuscl->index.density] = ucon[0] = vst->rho + offset;
	u[vmuscl->index.energy] = vst->e + offset;
	ucon[1] = vst->en_den + offset;
	for (i = 0; i < dim; ++i)    
	{
	    u[vmuscl->index.v[i]] = vst->v[i] + offset;
	    ucon[i+2] = vst->m[i] + offset;
	}
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* Include Mass fraction into vst */
            /* It is converted in g_load_VGas_state_vectors() */
            Gas_param *params = Params(vst->state[offset]);
            int    num_comps;
            if((num_comps = params->n_comps) != 1)
            {
                for(i = 0; i < num_comps; i++)
                {
                    ucon[vmuscl->index.prho[i]] = u[vmuscl->index.prho[i]]
                               = vst->rho0[i] + offset;
                }
            }
        }
	if ((src != NULL) && ((source = vmuscl->source) != NULL))
	{
	    source[0] = src->mass + offset;
	    source[1] = src->energy + offset;
	    for (i=0; i < dim; ++i)	 
	    	source[i+2] = src->mom[i] + offset;
            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
                /* Include Mass fraction into vst */
                /* It is converted in g_load_VGas_state_vectors() */
                Gas_param *params = Params(vst->state[offset]);
                int    num_comps, j;
                if((num_comps = params->n_comps) != 1)
                {
                    for(i = 0; i < num_comps; i++)
                    {
                        source[vmuscl->index.prho[i]] = src->rho0[i] + offset;
                        /* New 051005 */
                        for(j = 0; j < vsize; j++)
                            source[vmuscl->index.prho[i]][j] = 0.0;
                        /* End of New 051005 */
                    }
                }
            }
	}

	if (vmuscl->grav != NULL)
	{
	    float time = fr->time + 0.5*dt;
	    float **coords = vst->coords + offset;
	    float *g = vmuscl->grav + offset;
	    int   j;

	    if (tsten == NULL)
	    {
	        for (j = 0; j < vsize; ++j)
	            g[j] = gravity(coords[j],time)[iperm[swp_num]];
	    }
	    else
	    {
	    	for (j = 0; j < vsize; ++j)
	    	    g[j] = scalar_product(Q[0],gravity(coords[j],time),dim);
	    }
	}

	is_src = g_muscl_load_source_vectors(swp_num,iperm,tsten,vmuscl->grav,
			                     vst,src,offset,vsize,
					     fr->rect_grid);

	vmuscl->is_src = is_src;

#if defined(DEBUG_MUSCL)
	if (debugging("load_Vgas")) 
	{
	    (void) printf("All quantities before MUSCL computation:\n");

	    g_printout_vec_data("The state variables",
				u[vmuscl->index.density],
				u[vmuscl->index.energy],
				u+vmuscl->index.v[0],
				dim,0,vsize,"muncons");
	    g_printout_vec_data("The conserved variables",
				ucon[0],ucon[1],ucon+2,dim,0,vsize,"cons");
	    if (is_src)
	    {
	    	g_printout_vec_data("The source terms",
	    			    source[0],source[1],source+2,dim,1,
	    			    vsize-1,"src");
	    }
 
	    (void) printf("Miscelaneous quantities:\n");
	    (void) printf("%-4s %-14s %-14s %-14s\n","n","pressure",
	    	          "sound speed","Grun. coef");
	    for (i = offset; i < vsize; ++i) 
	        (void) printf("%-4d %-14g %-14g %-14g\n",
			      i,vst->p[i],vst->c[i],vst->GAM[i]);
	    (void) printf("\n");
	}
#endif /* defined(DEBUG_MUSCL) */
	return vmuscl;
}		/*end g_load_muscl_state_data*/


LOCAL bool g_muscl_load_source_vectors(
	int		   swp_num,
	int		   *iperm,
	Tan_stencil	   *tsten,
	float		   *grav,
	Vec_Gas		   *vst,
	Vec_Src		   *src,
	int		   offset,
	int		   vs,
	RECT_GRID	   *gr)
{
	int		   i, j;
	int		   dim = gr->dim;
	const float* const *Q = vst->Q;
	float		   *rho_src;
	float		   *e_src;
	float		   *v_src[SMAXD];
	float		   *v_src0, *v_srci;
	static	bool	   first = YES;
	float		*pr;
	float		*v0, *rho;
	static	float	alpha;

	if (src == NULL)
	    return NO;

	rho_src = src->mass + offset;
	e_src = src->energy + offset;

	for (i = 0; i < dim; ++i)
	    v_src[i] = src->mom[i] + offset;
	if (first == YES) 
	{
	    first = NO;
	    if (is_rotational_symmetry())
	    	alpha = rotational_symmetry();
	}

	for (j = 0; j < vs; ++j)
	{
	    rho_src[j] = 0.0;
	    e_src[j] = 0.0;
	}
	for (i = 0; i < dim; ++i)
	{
	    v_srci = v_src[i];
	    for (j = 0; j < vs; ++j)
		v_srci[j] = 0.0;
	}
	if (grav != NULL)
	{
	    float	*g = grav + offset;
	    v_src0 = v_src[0];
	    for (j = 0; j < vs; ++j)
	    	v_src0[j] += g[j];
	}

	if (tsten == NULL)
	{
	    /* Rect sweep */

	    if (is_rotational_symmetry() && 
		(alpha > 0.0) && (iperm[swp_num]==0))
	    {
		/* Include cylindrical source terms */
 
		float *radii = src->radii + offset;
		float rmin = src->rmin;
		float *m;

		m = vst->m[0] + offset;
		pr = vst->p + offset;
		rho = vst->rho + offset;
		v0 = vst->v[0] + offset;
		for (j = 0; j < vs; ++j)
		{
		    if (fabs(radii[j]) > fabs(rmin))
		    {
		        rho_src[j] -= alpha*m[j]/radii[j];
		        e_src[j] -= alpha*pr[j]*v0[j]/(radii[j]*rho[j]);
		    }
		}
	    }
	}
	else if (is_rotational_symmetry())
	{
	    /* Tangential sweep */

	    /* For cylindrical coordinates the mass equation for the
	    *  tangential sweep becomes
	    *
	    *	drho/dt + d(rho*vtan)/dr = - rho*vtan/r * dr/dtan
	    *
	    *  where dr/dtan is the change in radius with respect to
	    *  changes in the tangential coordinate.  If we let
	    *  (r, z) = tT + nN, i.e. represent r and z by a combination
	    *  of uni_arrays T & N tangent and normal to the curve, dr/dtan
	    *  becomes d(tT[0])/dt = T[0] ( = Q[0][0] )
	    */

	    if (alpha > 0.0)
	    {
	        float rmin, rad;
		float *mom0;
	        POINT **pt;

	        mom0 = vst->m[0];
	        pt = tsten->p;
	        pr = vst->p;
	        rho = vst->rho;
	        v0 = vst->v[0];

		rmin = src->rmin;
	        for (j = 0; j < vs; ++j)
	        {
	            rad = pos_radius(Coords(pt[j-2])[0],gr);
		    if (fabs(rad) > fabs(rmin))
		    {
	                rho_src[j] -= alpha*mom0[j]*Q[0][0]/rad;
	                e_src[j] -= alpha*Q[0][0]*v0[j]*pr[j]/(rho[j]*rad);
		    }
	        }
	    }
	}

#if defined(DEBUG_MUSCL)
	if (debugging("vsrc"))
	{
	    g_printout_vec_data("The source uni_arrays from "
	        	        "g_muscl_load_source_vectors:",
	        	        src->mass,src->energy,src->mom,dim,
	        	        offset,vs,"src");
	}
#endif /* defined(DEBUG_MUSCL) */
	return YES;
}		/*end g_muscl_load_source_vectors*/


/*
*			g_muscl_alloc_wk_space():
*
*	Allocates work space arrays for the uni_arrayized MUSCL finite
*	difference method.
*
*	The memory map is as follows:
*
*	pM = worksp[0][0]
*	backward = worksp[0]
*	uL[i] = worksp[0][i]
*
*	central = worksp[1]
*	uR[i] = worksp[1][i]
*
*	forward = worksp[2]
*	du = worksp[2]
*	right/left_wvs[i] = uM[i] = worksp[2][i]
*
*	source_sum = worksp[3]
*	q = worksp[3]
*	dulim = worksp[3]
*	F[i] = worksp[3][i]
*
*	dq[i] = worksp[4][i]
*/

LOCAL	Vec_Muscl *g_muscl_alloc_wk_space(
	Muscl_Opts *mopts,
	int	   vsize,
	int	   dim,
	Vec_Muscl  *vmuscl)
{
	float		 ***worksp;
	AVISC		 Avisc;
	Vec_Eigen	 *vegn;
	int		 nfloats = mopts->nfloats;
	int		 negn = 3;
	int		 worksp_len = mopts->worksp_len;
	int              i, nvar_u;

	vmuscl = alloc_Vec_Muscl(vmuscl);
	vegn = &vmuscl->Vegn;

	vmuscl->Opts = *mopts;
	vmuscl->dim = dim;

	nvar_u = 0;
	vmuscl->index.density = nvar_u++;
	vmuscl->index.energy = nvar_u++;
	for (i = 0; i < dim; ++i)
	    vmuscl->index.v[i] = nvar_u++;
	vmuscl->nvar_u = nvar_u;
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            int num_comp;
            num_comp = mopts->kmax - (dim+2);
            for(i = 0; i < num_comp; ++i)
                vmuscl->index.prho[i] = nvar_u++;
            vmuscl->nvar_u = nvar_u;
        }

	zero_scalar(vegn,sizeof(Vec_Eigen));
	vegn->negn = negn;
	set_no_alloc(vegn,vegn);
	MATRIX(vegn,lambda,negn,vsize,FLOAT);
	MATRIX(vegn,sgnl,negn,vsize,FLOAT);
	MATRIX(vegn,sgnr,negn,vsize,FLOAT);
	MATRIX(vegn,sgnm,negn,vsize,FLOAT);
	TRI_ARRAY(vegn,l,negn,negn,vsize,FLOAT);
	TRI_ARRAY(vegn,r,negn,negn,vsize,FLOAT);

	if (is_gravity() == YES)
	    VECTOR(vmuscl,grav,vsize,FLOAT);

	TRI_ARRAY(vmuscl,worksp,worksp_len,nfloats,vsize,FLOAT);
	worksp = vmuscl->worksp;

	VECTOR(vmuscl,u,nfloats,sizeof(float*));
	VECTOR(vmuscl,ucon,nfloats,sizeof(float*));
	if (source_terms_exist() == YES)
	    VECTOR(vmuscl,source,nfloats,sizeof(float*));

	use_artificial_dissipation(&Avisc);
	if (use_lapidus_artificial_viscosity(Avisc))
	{
	    scalar(&vmuscl->avisc,sizeof(Vec_Avisc));
	    set_alloc(vmuscl->avisc,avisc);
	    MATRIX(vmuscl->avisc,g,3,vsize,FLOAT);
	    ASSIGN_ARRAY_POINTER(vmuscl->avisc,cs_ave,worksp[0][0]);
	    ASSIGN_ARRAY_POINTER(vmuscl->avisc,c_ave,worksp[0][1]);
	    ASSIGN_ARRAY_POINTER(vmuscl->avisc,vn_ave,worksp[0][2]);
	    ASSIGN_ARRAY_POINTER(vmuscl->avisc,b,worksp[1]);
	    ASSIGN_ARRAY_POINTER(vmuscl->avisc,visc,worksp[2][0]);
	    ASSIGN_ARRAY_POINTER(vmuscl->avisc,uconM,worksp[1]);
	    set_no_alloc(vmuscl->avisc,mdlambda);
	}
	else if (use_upwind_artificial_viscosity(Avisc) ||
	         use_linear_artificial_viscosity(Avisc))
	{
	    scalar(&vmuscl->avisc,sizeof(Vec_Avisc));
	    set_alloc(vmuscl->avisc,avisc);

/*          041503, change to old code fashion
	    MATRIX(vmuscl->avisc,g,1,vsize,sizeof(float*));
*/ 
            MATRIX(vmuscl->avisc,g,1,vsize,sizeof(float));
/* end of change 041503 */

	    set_no_alloc(vmuscl->avisc,cs_ave);
	    set_no_alloc(vmuscl->avisc,c_ave);
	    set_no_alloc(vmuscl->avisc,vn_ave);
	    set_no_alloc(vmuscl->avisc,b);
	    set_no_alloc(vmuscl->avisc,visc);
	    set_no_alloc(vmuscl->avisc,uconM);
	    ASSIGN_ARRAY_POINTER(vmuscl->avisc,mdlambda,worksp[0][0]);
	}

	if (vmuscl->avisc != NULL)
	{
	    vmuscl->avisc->use_lapidus=use_lapidus_artificial_viscosity(Avisc);
	    vmuscl->avisc->use_linear = use_linear_artificial_viscosity(Avisc);
	    vmuscl->avisc->use_upwind = use_upwind_artificial_viscosity(Avisc);
	}

	if (use_muscl_slope_flattening(Avisc))
	{
	    scalar(&vmuscl->msf,sizeof(Vec_MSF));
	    set_alloc(vmuscl->msf,msf);
	    VECTOR(vmuscl->msf,chi,vsize,FLOAT);
	}
	else
	    vmuscl->msf = NULL;
	vmuscl->monotone_reconstruction = mopts->monotone_reconstruction;
	vmuscl->link_reconstructions = mopts->link_reconstructions;

	vmuscl->max_vsize = vsize;

	        /* Set Linear reconstruction data structure */
	ASSIGN_ARRAY_POINTER(vmuscl,backward,worksp[0]);
	ASSIGN_ARRAY_POINTER(vmuscl,central,worksp[1]);
	ASSIGN_ARRAY_POINTER(vmuscl,forward,worksp[2]);
	ASSIGN_ARRAY_POINTER(vmuscl,du,worksp[2]);
	ASSIGN_ARRAY_POINTER(vmuscl,q,worksp[3]);

	        /* Set half step calculation data */

	ASSIGN_ARRAY_POINTER(vmuscl,uL,worksp[0]);
	ASSIGN_ARRAY_POINTER(vmuscl,uR,worksp[1]);    
	ASSIGN_ARRAY_POINTER(vmuscl,uM,worksp[2]);
	ASSIGN_ARRAY_POINTER(vmuscl,right_wvs,worksp[2]);
	ASSIGN_ARRAY_POINTER(vmuscl,left_wvs,worksp[2]);
	ASSIGN_ARRAY_POINTER(vmuscl,Flux.F,worksp[3]);
	ASSIGN_ARRAY_POINTER(vmuscl,source_sum,worksp[3]);
	ASSIGN_ARRAY_POINTER(vmuscl,awv,worksp[3]);
	ASSIGN_ARRAY_POINTER(vmuscl,dq,worksp[4]);

	        /* Set Riemann solver date */
	vmuscl->pL = g_set_vmuscl_vec_gas_storage(vmuscl->uL,&vmuscl->VLst,
	        	                          vmuscl);
	vmuscl->pR = g_set_vmuscl_vec_gas_storage(vmuscl->uR,&vmuscl->VRst,
	        	                          vmuscl);
	vmuscl->pM = g_set_vmuscl_vec_gas_storage(vmuscl->uM,&vmuscl->VMst,
	        	                          vmuscl);
	set_no_alloc(vmuscl,A);
	set_no_alloc(vmuscl,dV);
	return vmuscl;
}		/*end g_muscl_alloc_wk_space*/


LOCAL	float *g_set_vmuscl_vec_gas_storage(
	float	  **u,
	Vec_Gas	  *vst,
	Vec_Muscl *vmuscl)
{
	int	  dim = vmuscl->dim;
	int	  vsize = vmuscl->max_vsize;
	int i;

	zero_scalar(vst,sizeof(Vec_Gas));
	ASSIGN_ARRAY_POINTER(vst,rho,u[vmuscl->index.density]);
	ASSIGN_ARRAY_POINTER(vst,e,u[vmuscl->index.energy]);
	VECTOR(vst,v,dim,sizeof(float*));
	for (i = 0; i < dim; ++i)
	    vst->v[i] = u[vmuscl->index.v[i]];
	VECTOR(vst,prms_jmp,vsize+1,INT);
	VECTOR(vst,p,vsize,FLOAT);
	VECTOR(vst,re,vsize,FLOAT);
	VECTOR(vst,GAM,vsize,FLOAT);
	VECTOR(vst,FD,vsize,FLOAT);
	VECTOR(vst,c2,vsize,FLOAT);
	VECTOR(vst,c,vsize,FLOAT);
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            int num_comp;
            num_comp = vmuscl->Opts.kmax - (dim+2);
            VECTOR(vst,rho0,num_comp,sizeof(float*));
            for(i = 0; i < num_comp; ++i)
                vst->rho0[i] = u[vmuscl->index.prho[i]];
        }
	return vst->p;
}		/*end g_set_vmuscl_vec_gas_storage*/


/*
*			g_muscl_flux_vectors();
*
*	This function computes the conservative fluxes for the system.
*/

LOCAL void g_muscl_flux_vectors(
	int	   start,
	int	   end,
	float	   **u,
	float	   *p,
	MUSCL_FLUX *Flux,
	Vec_Muscl  *vmuscl)
{
	float 	*rho = u[vmuscl->index.density];
	float 	*e = u[vmuscl->index.energy];
	float 	*v0, *v1, *v2;
	float 	*F0, *F1, *F2, *F3, *F4;
	float 	**F = Flux->F;
	int   	dim = vmuscl->dim;
	int   	i,j;

	float *coords;
    	int    *icoords;
    	float GL[3],GU[3],dh[3], uvel[3];
    	for(i = 0; i< 3; i++)
   	{
       	 GL[i] = vmuscl->GL[i];
         GU[i] = vmuscl->GU[i];
         dh[i] = vmuscl->dh[i];
        }
	float	wall_z[2];
	wall_z[0] = GL[2] + 0.9;
	wall_z[1] = GU[2] - 0.1;
        float 	dist_z[2];
        float 	T, P;
        float 	flux[4];
        int    	idir = vmuscl->idir;
        float 	delta;
	int    	iperm[3], swp_num;
	int	idirs[3];
	int	irdirs[3];

	switch (dim)
	{
	case 1:
	    v0 = u[vmuscl->index.v[0]];
	    F0  = F[0]; F1 = F[1]; F2 = F[2];
	    for (j = start; j < end; ++j)
	    {
	    	/* mass */
	        F0[j] = rho[j]*v0[j];
	    	/* energy */
	        F1[j] = v0[j]*(rho[j]*(0.5*sqr(v0[j]) + e[j]) + p[j]);
	    	/* momentum */
	        F2[j] = F0[j]*v0[j] + p[j];
                if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                {
                    /* Partial mass flux = total_mass_flux * mass_fraction*/
                    Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
                    int       i, num_comps;
                    if((num_comps = params->n_comps) != 1)
                    {
                        for(i = 0; i < num_comps; i++)
                            F[vmuscl->index.prho[i]][j] = F0[j]*u[vmuscl->index.prho[i]][j];
                    }
                }
	    }
	    break;

	case 2:
	    v0 = u[vmuscl->index.v[0]]; v1 = u[vmuscl->index.v[1]];
	    F0  = F[0]; 
	    F1 = F[1]; 
	    F2 = F[2];
	    F3 = F[3];
	    for (j = start-1; j < end; ++j)
	    {
	    	/* mass */
	        F0[j] = rho[j]*v0[j];
	    	/* energy */
	        F1[j] = v0[j]*
	            (rho[j]*(0.5*(sqr(v0[j])+sqr(v1[j]))+e[j])+p[j]);
	    	/* Sweep component of momentum */
	        F2[j] = F0[j]*v0[j] + p[j];
	    	/* Off sweep component of momentum */
	        F3[j] = F0[j]*v1[j];
                if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                {
                    /* Partial mass flux = total_mass_flux * mass_fraction*/
                    Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
                    int       i, num_comps;
                    if((num_comps = params->n_comps) != 1)
                    {
                        for(i = 0; i < num_comps; i++)
                            F[vmuscl->index.prho[i]][j] = F0[j]*u[vmuscl->index.prho[i]][j];
                    }
                }
	    }
	    break;

	case 3:
	    v0 = u[vmuscl->index.v[0]];
	    v1 = u[vmuscl->index.v[1]];
	    v2 = u[vmuscl->index.v[2]];
	    F0 = F[0]; 
	    F1 = F[1]; 
	    F2 = F[2];
	    F3 = F[3];
	    F4 = F[4];

	    for (j = start; j < end; ++j)
	    {

			/* mass */
		        F0[j] = rho[j]*v0[j];
		    	/* energy */
		        F1[j] = v0[j]*
		    	    (rho[j]*(0.5*(sqr(v0[j])+sqr(v1[j])+sqr(v2[j]))+e[j])+p[j]);
		    	/* Sweep component of momentum */
	       		F2[j] = F0[j]*v0[j] + p[j];
		    	/* Off sweep component of momentum */
		        F3[j] = F0[j]*v1[j];
		    	/* Off sweep component of momentum */
		        F4[j] = F0[j]*v2[j];
			
			coords = vmuscl->coords[j];
			if (coords[0] > ramp_x1)
                     		dist_z[0] = coords[2] - ramp_z1;
                	else if( coords[0] < ramp_x0)
                     		dist_z[0] = coords[2] - ramp_z0 ;
                	else dist_z[0] = coords[2] - (ramp_z0 + tan10*(coords[0] - ramp_x0));
			/*Find the direction after permutation.*/
                	iperm[1] = vmuscl->iperm[1];
                	iperm[0] = vmuscl->iperm[0];
                	iperm[2] = vmuscl->iperm[2];
                	swp_num = vmuscl->swp_num;
                	if(swp_num != -1)
                	{
                    	for (i = 0; i < dim; ++i) idirs[i] = iperm[(i+swp_num)%dim];
                	}

		        if(vmuscl->turb_boundary_layer == YES && idirs[0] == 2 && coords[0] > ramp_x0 && coords[0] < ramp_x1 && dist_z[0] > 0 && dist_z[0] <= dh[2] && swp_num != -1)
			{
		    		/*Set the physical flux across a wall to be zero here.*/
		    		F0[j] = 0;
	            		F1[j] = 0;
		    		F2[j] = p[j];
		    		F3[j] = 0;
		    		F4[j] = 0;
			}
       		        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
               		{
                    	/* Partial mass flux = total_mass_flux * mass_fraction*/
                    		Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
                    		int       i, num_comps;
                  		if((num_comps = params->n_comps) != 1)
                    		{
	                        for(i = 0; i < num_comps; i++)
        	                    F[vmuscl->index.prho[i]][j] = F0[j]*u[vmuscl->index.prho[i]][j];
                		 }
                	 }
	
	
              	}
	    break;

	default:
	    screen("ERROR in g_muscl_flux_vectors(), "
	           "invalid dimension %d\n",dim);
	    clean_up(ERROR);
	    break;
	}
	debug_print("mflux","Left g_muscl_flux_vectors()\n");
}		/*end g_muscl_flux_vectors*/


/*
*			g_cons_src():
*
*	Transforms source terms from nonconservative to conservative form.
*	This function computes the conservative form of the source terms
*	for the Euler equations at the half-step level.  source*dt
*	is a approximation of the integral of the source vector over the
*	space/time cell [x[j]-dx/2,x[j]+dx/2]X[t,t+dt].
*
*	The macros below define the transformation from the nonconservation
*	representation of the Euler equations using density, specific
*	internal energy, and velocity into the conservative form 
*	using density, total energy, and momentum.
*/

LOCAL void g_cons_src(
	int		start,
	int		end,
	int		swp_num,
	int		*iperm,
	Tan_stencil	*tsten,
	Vec_Muscl	*vmuscl)
{
	float		*grav = vmuscl->grav;
	const float* const *Q = vmuscl->Q;
	int		i, j,dim;
	float		*rho_src, *E_src, *m_src[MAXD], *m_srci;
	float		**uM, **source;
	float		*rho;
	float		*v0;
	static bool	first = YES;
	float		*vi;
	float		**F;
	float		*rho_flux, *E_flux;
	static	float	alpha;
        float           *prho0_src[MAX_NUM_GAS_COMPS], *prho0_flux[MAX_NUM_GAS_COMPS];
        int             n_comps = 1;

	if (((source = vmuscl->source) == NULL) || (vmuscl->is_src == NO)) 
	    return;

	dim = vmuscl->dim;
	uM = vmuscl->uM;
	rho = uM[0];
	v0 = uM[2];
	rho_src = source[0];
	E_src = source[1];
	for (i = 0; i < dim; ++i)
	    m_src[i] = source[i+2];
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* hook up partial density pointers */
            Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
            if((n_comps = params->n_comps) != 1)
            {
                for(i = 0; i < params->n_comps; i++)
                    prho0_src[i] = source[vmuscl->index.prho[i]];
            }
        }

	if (first == YES)
	{
	    first = NO;
	    if (is_rotational_symmetry())
		alpha = rotational_symmetry();
	}
	for (j = start; j <= end; ++j)
	{
	    rho_src[j] = 0.0;
	    E_src[j] = 0.0;
	}
	for (i = 0; i < dim; ++i)
	{
	    m_srci = m_src[i];
	    for (j = start; j <= end; ++j)
	    	m_srci[j] = 0.0;
	}

        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* hook up partial density pointers */
            Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
            if((n_comps = params->n_comps) != 1)
            {
		for( j = start; j<=end; ++j)
			for( i=0;i < params->n_comps; i++)
		 		prho0_src[i][j] = 0.0;
		
		/*ReadThermoData();
                for (j = start; j <= end; ++j)
                {
		    ComputeChemSourceTerm(vmuscl->vst, vmuscl->src,temperature(vmuscl->vst->state[j]),vmuscl->vst->rho[j],j,n_comps);
                    for(i = 0; i < params->n_comps; i++){
                        prho0_src[i][j] += (1.0/3)*(vmuscl->src->rho0[i][j]);
		    }
		    for(i=0; i<NUM_OF_REACTIONS; ++i)
		    {
			ComputeHeatRelease(i, temperature(vmuscl->vst->state[j]) );
			E_src[j] += (1.0/3)*reaction[i].HeatRelease;
		    }
                }*/
            }
        }
	if (tsten == NULL)
	{
	    /* Rect sweep */

	    if (is_rotational_symmetry() &&
		(alpha > 0.0) && (iperm[swp_num]==0))
	    {
	        float	rmin, *radii;
	    	/* Include cylindrical source terms */

		rmin = vmuscl->src->rmin;
	    	radii = vmuscl->src->radii;
	    	F = vmuscl->Flux.F;
	    	rho_flux = F[0];
	    	E_flux = F[1];
                if(g_composition_type() == MULTI_COMP_NON_REACTIVE &&
                   n_comps != 1)
                {
                    for(i = 0; i < n_comps; i++)
                        prho0_flux[i] = F[vmuscl->index.prho[i]];
                }
	    	for (j = start; j < end; ++j)
	    	{
		    if (fabs(radii[j]) > fabs(rmin))
		    {
	    	        rho_src[j] -= alpha*0.5*(rho_flux[j]+rho_flux[j+1])/
				      radii[j];
	    	        E_src[j] -= alpha*0.5*(E_flux[j]+E_flux[j+1])/radii[j];
		    }
	        }
                if(g_composition_type() == MULTI_COMP_NON_REACTIVE &&
                       n_comps != 1)
                {
                    for(i = 0; i < n_comps; i++)
                    {
                        for(j = start; j < end; ++j)
                        {
                            if (fabs(radii[j]) > fabs(rmin))
                            {
                                prho0_src[i][j] -= alpha*0.5*(prho0_flux[i][j]
					+prho0_flux[i][j+1])/radii[j];
                            }
                        }
                    }
                }
	    	for (i = 0; i < dim; ++i)
	    	{
	    	    m_srci = m_src[i];
	    	    vi = uM[2+i];
	    	    for (j = start; j < end; ++j)
		    {
	    	    	m_srci[j] += rho_src[j]*0.5*(vi[j]+vi[j+1]);
		    }
	    	}
	    }
	    if (grav != NULL)
	    {
	    	m_srci = m_src[0];
	    	for (j = start; j < end; ++j)
	    	{
	    	    m_srci[j] += 0.5*(rho[j] + rho[j+1])*grav[j];
	    	    E_src[j] += 0.5*(rho[j]*v0[j] + rho[j+1]*v0[j+1])*grav[j];
	    	}
	    }
	}
	else
	{
	    /* Tangential sweep */

	    if (is_rotational_symmetry() &&
		alpha > 0.0)
	    {
	    	POINT	**pt;
	        RECT_GRID *gr = vmuscl->front->rect_grid;
	        float	rmin, rad;
	        float	a = 0.5*alpha*Q[0][0];

	    	/* Include cylindrical source terms */
	    	F = vmuscl->Flux.F;
	    	rho_flux = F[0];
	    	E_flux = F[1];
                if(g_composition_type() == MULTI_COMP_NON_REACTIVE &&
                   n_comps != 1)
                {
                    for(i = 0; i < n_comps; i++)
                        prho0_flux[i] = F[vmuscl->index.prho[i]];
                }
	    	pt = tsten->p;
		rmin = pos_radius(0.0,gr);
	    	for (j = start; j < end; ++j)
	    	{
	            rad = pos_radius(Coords(pt[j-2])[0],gr);
		    if (fabs(rad) > fabs(rmin))
		    {
	                rho_src[j] -= a*(rho_flux[j]+rho_flux[j+1])/rad;
	                E_src[j] -= a*(E_flux[j]+E_flux[j+1])/rad;
		    }
	        }
                if(g_composition_type() == MULTI_COMP_NON_REACTIVE &&
                   n_comps != 1)
                {
                    for(i = 0; i < n_comps; i++)
                    {
                        for(j = start; j < end; ++j)
                        {
                            rad = pos_radius(Coords(pt[j-2])[0],gr);
                            if (fabs(rad) > fabs(rmin))
                                prho0_src[i][j] -= a*(prho0_flux[i][j]+prho0_flux[i][j+1])/rad;
                        }
                    }
                }
	        for (i = 0; i < dim; ++i)
	        {
	            m_srci = m_src[i];
	            vi = uM[2+i];
	            for (j = start; j < end; ++j)
                    {
	      	        m_srci[j] += rho_src[j]*0.5*(vi[j]+vi[j+1]);
                    }
	        }
	    }

	    if (grav != NULL)
	    {
	    	m_srci = m_src[0];
	    	for (j = start; j < end; ++j)
	        {
	            m_srci[j] += 0.5*(rho[j] + rho[j+1])*grav[j];
	            E_src[j] += 0.5*(rho[j]*v0[j] + rho[j+1]*v0[j+1])*grav[j];
	        }
	    }
	}
#if defined(DEBUG_MUSCL)
	if (debugging("csrc"))
	{
	    g_printout_vec_data("Here are the conservative source uni_arrays",
	    	                source[0],source[1],source+2,
	    	                vmuscl->dim,start,end,"src");
	}
#endif /* defined(DEBUG_MUSCL) */
}		/*end g_cons_src*/


/*
*			g_state_linear_reconstructor():
*
*	This subroutine constructs piecewise linear functions for the
*	cell averages.	The construction is processed along the direction
*	of the right eigenuni_arrays.  A 'limiter' is applied to avoid
*	oscillation.  The slopes of the piecewise linear function in the
*	directions of the eigenuni_arrays are given by:
*
*	slope[m][j] = du[][j] * l[m][][j]
*
*	du[m][j] = min( 2 * |b[m][j]|, .5 * |c[m][j]|, 2 * |f[m][j]|) *
*					sgn(c[m][j])
*			if sgn(b[m][j]) == sgn(c[m][j]) == sgn(f[m][j])
*		   = 0		otherwise.
*
*	where
*	b[m][j] = (u[m][j]   - u[m][j-1])
*	c[m][j] = (u[m][j+1] - u[m][j-1])
*	f[m][j] = (u[m][j+1] - u[m][j]).
*
*	This limiter can be found (in scalar form) in van Leer, J. Comp.
*	Physics 23, p. 289 (1977).
*
*	Once again we are not introducing any notation for the last two
*	eigenvalues or eigenuni_arrays.
*
*	When using g_load_muscl_state_data() to load the state uni_arrays,
*	we have
*           u[vmuscl->index.density] = rho,
*           u[vmuscl->index.energy] = e,
*           u[vmuscl->index.v[i]] = velocity(i), 0 <= i < dim
*/

LOCAL void g_state_linear_reconstructor(
	int		start,
	int		end,
	Vec_Muscl	*vmuscl)
{
	int		dim;
	int		j, k;	/* j is the mesh index */
	int		kmax;
	float		**u, **dq, **backward, **central, **forward, **du;
	float		***l,***r;
	float		*uk, *bk, *ck, *fk;
	int		num_comps = 0;

#if defined(DEBUG_MUSCL)
	debug_print("lcnst","Entered g_state_linear_reconstructor()\n");
#endif /* defined(DEBUG_MUSCL) */

	dim =		vmuscl->dim;
	kmax =		dim+2;
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* Compute differences for mass fraction */
            Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
            if(params->n_comps != 1)
            {
                num_comps = params->n_comps;
                kmax += num_comps;
            }
        }
	u =		vmuscl->u;
	l =		vmuscl->Vegn.l;
	r =		vmuscl->Vegn.r;
	backward =	vmuscl->backward;
	central =	vmuscl->central;
	forward =	vmuscl->forward;
	du =		vmuscl->du;
	dq =		vmuscl->dq;

	for (k = 0; k < kmax; ++k)
	{
	    bk = backward[k]; ck = central[k]; fk = forward[k];
	    uk = u[k];
	    for (j = start; j < end; ++j)
	    {
	    	bk[j] =	      (uk[j]	 - uk[j - 1]);
	    	ck[j] = 0.5 * (uk[j + 1] - uk[j - 1]);
	    	fk[j] =	      (uk[j + 1] - uk[j]);
	    }
	}

	apply_limiter(du,start,end,vmuscl);
	
	/* New 051005 change argument "dim" to "dim+num_comps" */
	left_multiply_state_by_matrix(dq,l,du,dim+num_comps,start,end);
/*
	limit_strong_wave_in_cell(vmuscl,start,end);
*/
#if defined(DEBUG_MUSCL)
	if (debugging("lcnst"))
	    print_linear_reconstruction(start, end, vmuscl);
	debug_print("lcnst","Left g_state_linear_reconstructor()\n");
#endif /* defined(DEBUG_MUSCL) */
}		/*end g_state_linear_reconstructor*/


/*
*			g_bct_linear_reconstructor():
*
*
*    This subroutine constructs piecewise linear functions for the
*    cell averages.  The construction is processed along the direction
*    of the right eigenuni_arrays.	 A 'limiter' is applied to avoid
*    oscillation.  The slopes of the piecewise linear function in the
*    directions of the eigenuni_arrays are given by:
*
*	slope[m][j] = min( 2 * |b[m][j]|, .5 * |c[m][j]|, 2 * |f[m][j]|) *
*					sgn(c[m][j])
*			if sgn(b[m][j]) == sgn(c[m][j]) == sgn(f[m][j])
*		   = 0		otherwise.
*
*	where
*	b[m][j] = (u[][j]   - u[][j-1]) dot l[m][][j]
*	c[m][j] = (u[][j+1] - u[][j-1]) dot l[m][][j]
*	f[m][j] = (u[][j+1] - u[][j]  ) dot l[m][][j].
*
*	This limiter can be found (in scalar form) in van Leer, J. Comp.
*	Physics 23, p. 289 (1977).
*
*	Once again we are not introducing any notation for the last two
*	eigenvalues or eigenuni_arrays.
*/

LOCAL void g_bct_linear_reconstructor(
	int		start,
	int		end,
	Vec_Muscl	*vmuscl)
{
	int		dim;
	int		j, k;	/* j is the mesh index */
	int		kmax;
	float		**u, **backward, **central, **forward;
	float		***l;
	float		*bk, *ck, *fk;
	float		*u0, *u1, *u2, *uk;
	float		*lk0, *lk1, *lk2;

#if defined(DEBUG_MUSCL)
	debug_print("lcnst","Entered g_bct_linear_reconstructor()\n");
#endif /* defined(DEBUG_MUSCL) */

	compute_slope_limiting_coeffs(start, end, vmuscl);

	dim =		vmuscl->dim;
	kmax =		dim+2;
	u =		vmuscl->u;
	l =		vmuscl->Vegn.l;
	backward =	vmuscl->backward;
	central =	vmuscl->central;
	forward =	vmuscl->forward;

	u0 = u[0];	u1 = u[1];	u2 = u[2];
	for (k = 0; k < 3; ++k)
	{
	    bk = backward[k]; ck = central[k]; fk = forward[k];
	    lk0 = l[k][0]; lk1 = l[k][1]; lk2 = l[k][2];
	    for (j = start; j < end; ++j)
	    {
	    	bk[j] = ((u0[j] - u0[j - 1])	* lk0[j] +
		    	 (u1[j] - u1[j - 1])	* lk1[j] +
	    		 (u2[j] - u2[j - 1])	* lk2[j]);
		ck[j] = 0.5 * (
			(u0[j + 1] - u0[j - 1]) * lk0[j] +
			(u1[j + 1] - u1[j - 1]) * lk1[j] +
			(u2[j + 1] - u2[j - 1]) * lk2[j]);
		fk[j] = (
			(u0[j + 1] - u0[j])	* lk0[j] +
			(u1[j + 1] - u1[j])	* lk1[j] +
			(u2[j + 1] - u2[j])	* lk2[j]);
	    }
	}
	for (k = 3; k < kmax; ++k)
	{
	    bk = backward[k]; ck = central[k]; fk = forward[k];
	    uk = u[k];
	    for (j = start; j < end; ++j)
	    {
	    	bk[j] =	      (uk[j]	 - uk[j - 1]);
	    	ck[j] = 0.5 * (uk[j + 1] - uk[j - 1]);
	    	fk[j] =	      (uk[j + 1] - uk[j]);
	    }
	}

        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* Compute differences for mass fraction */
            Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
            int    num_comps;
            if((num_comps = params->n_comps) != 1)
            {
                for (k = vmuscl->index.prho[0]; k < vmuscl->index.prho[0]+num_comps; ++k)
                {
                    bk = backward[k]; ck = central[k]; fk = forward[k];
                    uk = u[k];
                    for (j = start; j < end; ++j)
                    {
                        bk[j] =       (uk[j]     - uk[j - 1]);
                        ck[j] = 0.5 * (uk[j + 1] - uk[j - 1]);
                        fk[j] =       (uk[j + 1] - uk[j]);
                    }
                }
            }
        }
	apply_limiter(vmuscl->dq,start,end,vmuscl);

#if defined(DEBUG_MUSCL)
	if (debugging("lcnst"))
	    print_linear_reconstruction(start, end, vmuscl);
	debug_print("lcnst","Left g_bct_linear_reconstructor()\n");
#endif /* defined(DEBUG_MUSCL) */
}		/*end g_bct_linear_reconstructor*/

/*
*			g_eigenstate_linear_reconstructor():
*
*	This subroutine constructs piecewise linear functions for the
*	cell averages.	The construction is processed along the direction
*	of the right eigenuni_arrays.  A 'limiter' is applied to avoid
*	oscillation.  The slopes of the piecewise linear functions in the
*	directions of the eigenuni_arrays are given by:
*
*	slope[m][j] = min(2 * |b[m][j]|, .5 * |c[m][j]|, |f[m][j]|) *
*				sgn(c[m][j])
*			if sgn(b[m][j]) == sgn(c[m][j]) == sgn(f[m][j])
*		= 0	otherwise.
*
*	where
*	b[m][j] = (q[m][j]   - q[m][j-1])
*	c[m][j] = (q[m][j+1] - q[m][j-1]),
*	f[m][j] = (q[m][j+1] - q[m][j]	)
*
*	and q[m][j] = u[][j] dot l[m][][j]
*
*
*	This limiter can be found (in scalar form) in van Leer, J. Comp.
*	Physics 23, p. 289 (1977).
*
*	Once again we are not introducing any notation for the last two
*	eigenvalues or eigenuni_arrays.
*/

LOCAL void g_eigenstate_linear_reconstructor(
	int		start,
	int		end,
	Vec_Muscl	*vmuscl)
{
	int		dim;
	int		j, k;	/* j is the mesh index */
	int		kmax;
	float		**q = vmuscl->q, *qk;
	float		**u, **backward, **central, **forward;
	float		***l;
	float		*bk, *ck, *fk;

#if defined(DEBUG_MUSCL)
	debug_print("lcnst","Entered g_eigenstate_linear_reconstructor()\n");
#endif /* defined(DEBUG_MUSCL) */

	dim =		vmuscl->dim;
	kmax =		dim+2;
	u =		vmuscl->u;
	l =		vmuscl->Vegn.l;
	backward =	vmuscl->backward;
	central =	vmuscl->central;
	forward =	vmuscl->forward;

	left_multiply_state_by_matrix(q,l,u,dim,start,end);
	for (k = 0; k < kmax; ++k)
	{
	    bk = backward[k]; ck = central[k]; fk = forward[k];
	    qk = q[k];
	    for (j = start; j < end; ++j)
	    {
	    	bk[j] =	      (qk[j]	 - qk[j - 1]);
	    	ck[j] = 0.5 * (qk[j + 1] - qk[j - 1]);
	    	fk[j] =	      (qk[j + 1] - qk[j]);
	    }
	}

	apply_limiter(vmuscl->dq,start,end,vmuscl);

#if defined(DEBUG_MUSCL)
	if (debugging("lcnst"))
	    print_linear_reconstruction(start, end, vmuscl);
	debug_print("lcnst","Left g_eigenstate_linear_reconstructor()\n");
#endif /* defined(DEBUG_MUSCL) */
}		/*end g_eigenstate_linear_reconstructor*/

LOCAL	void	apply_limiter(
	float		**slope,
	int		start,
	int		end,
	Vec_Muscl	*vmuscl)
{
	int		dim;
	int		j, k;    /* j is the mesh index */
	int		kmax;
	float		**backward, **central, **forward;
	float		*bk, *ck, *fk;
	float		*slopek;
	float		sign, dl, abk, afk, ack;

	dim =		vmuscl->dim;
	kmax =		dim+2;
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* Compute differences for mass fraction */
            Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
            int    num_comps;
            if((num_comps = params->n_comps) != 1)
                kmax += num_comps;
        }
	backward =	vmuscl->backward;
	central =	vmuscl->central;
	forward =	vmuscl->forward;

#if defined(DEBUG_MUSCL)
	if (debugging("lcnst"))
	{
	    (void) printf("\n");
	    g_printout_vec_data("Backward differences",backward[0],backward[1],
	        		backward+2,dim,start,end,"muncons");
	    g_printout_vec_data("Central differences",central[0],central[1],
	        		central+2,dim,start,end,"muncons");
	    g_printout_vec_data("Forward differences",forward[0],forward[1],
	        		forward+2,dim,start,end,"muncons");
	}
#endif /* defined(DEBUG_MUSCL) */

	compute_slope_limiting_coeffs(start,end,vmuscl);

	/* Compute the standard Van Leer slope limiter */

	for (k = 0; k < kmax; ++k)
	{
	    bk = backward[k]; ck = central[k]; fk = forward[k];
	    slopek = slope[k];
	    for (j = start; j < end; ++j)
	    {
	    	if (bk[j]*fk[j] > 0.0)
	    	{
	    	    abk = fabs(bk[j]);
	    	    afk = fabs(fk[j]);
	    	    dl = 2.0*min(abk,afk);
	    	    sign = (ck[j] == 0.0) ? 0.0 : (ck[j] > 0.0) ? 1.0 : -1.0;
	    	    ack = fabs(ck[j]);
	    	    slopek[j] = sign*min(ack,dl);
	    	}
	    	else
	    	{
	    	    slopek[j] = 0.0;
	    	}
	    }
	}
	if (vmuscl->msf != NULL)
	{
	    float *chi = vmuscl->msf->chi;
	    for (k = 0; k < kmax; ++k)
	    {
	    	slopek = slope[k];
	    	for (j = start; j < end; ++j)
	    	    slopek[j] *= chi[j];
	    }
	}
	   
}		/*end apply_limiter*/
	

LOCAL void g_first_order_godunov_reconstructor(
	int		start,
	int		end,
	Vec_Muscl	*vmuscl)
{
	int		dim;
	int		j, k;	/* j is the mesh index */
	int		kmax;
	float		**dq, *dqk;

#if defined(DEBUG_MUSCL)
	debug_print("lcnst","Entered g_first_order_godunov_reconstructor()\n");
#endif /* defined(DEBUG_MUSCL) */

	dim =		vmuscl->dim;
	kmax =		dim+2;
	dq =		vmuscl->dq;

	for (k = 0; k < kmax; ++k)
	{
	    dqk = dq[k];
	    for (j = start; j < end; ++j)
	        dqk[j] = 0.0;
	}

#if defined(DEBUG_MUSCL)
	if (debugging("lcnst"))
	    print_linear_reconstruction(start, end, vmuscl);
	debug_print("lcnst","Left g_first_order_godunov_reconstructor()\n");
#endif /* defined(DEBUG_MUSCL) */
}		/*end g_first_order_godunov_reconstructor*/


/*
*			 g_compute_eigens():
*
*	Calculate the eigenvalues and eigenuni_arrays of the linearized gas
*	dynamical equations in the normal/tangential direction.
*	The state vector is u = (rho,e,v0,v1,v2), which are to be
*	interpreted as cell averages and with all
*	uni_arrays written in these coordinates.  v0 is always the normal
*	velocity.
*
*
*	The eigenvalues are:
*
*		lambda[0] = v0 - c
*		lambda[1] = lambda[3] = lambda[4] = v0
*		lambda[2] = v0 + c
*
*	c2 = sound speed squared, GAM = Gruneisen exponent.
*	The right eigenuni_arrays are (column uni_arrays):
*
*	      r[][0]	       r[][1]	          r[][2]   r[][3] r[][4]
*
*	    ( -rho/c	         GAM/c       	   rho/c      0	     0	 )
*	    (-p/(rho*c) (GAM*p/rho-c2)/(rho*c)   p/(rho*c)    0	     0	 )
*	    (	 1	         0	            1	      0	     0	 )
*	    (	 0	         0	            0	      1	     0	 )
*	    (	 0	         0	            0	      0	     1	 )
*
*	The left eigenuni_arrays are (row uni_arrays):
*
*	l[0][] ( (GAM*p/rho-c2)(2*rho*c)  -GAM/(2c)  1/2      0	     0	 )
*	l[1][] (    p/(rho*c)	          -rho/c      0	      0	     0	 )
*	l[2][] (-(GAM*p/rho-c2)/(2*rho*c)  GAM/(2c)  1/2      0	     0	 )
*	l[3][] (       0		      0	      0	      1	     0	 )
*	l[4][] (       0		      0	      0	      0	     1	 )
*
*	Since l3 and l4 (resp. r3 and r4) are perpendicular to l0, l1, l2
*	(resp. r0, r1, r2), we shall not use the third or fourth 
*	components of l0, l1, l2 (resp. r0, r1, r2). We shall also not
*	introduce the notations l3, r3, lambda3, l4, r4, or lambda4.
*
*	Returns the maximum absolute value of the eigenvalues.
*/


LOCAL void g_compute_eigens(
	int		start,
	int		end,
	Vec_Muscl	*vmuscl)
{
	Vec_Eigen  *vegn = &vmuscl->Vegn;
	Vec_Gas    *vst = vmuscl->vst;
	int        offset = vmuscl->offset;
	int        j, k;
	float      **lambda = vegn->lambda;
	float      ***l = vegn->l, ***r = vegn->r;
	float      *rho = vst->rho + offset;
	float      *v0 = vst->v[0] + offset;
	float      *c = vst->c + offset;
	float      *p = vst->p + offset;
	float      *GAM = vst->GAM + offset;
	float      *c2 = vst->c2 + offset;
	float      **sgnl = vegn->sgnl;
	float      **sgnr = vegn->sgnr;
	float      speed;
	float      *sgnlk, *sgnrk;
	float      *lambdak;
	float      *lambda0, *lambda1, *lambda2;
	float      *r00, *r01, *r02;
	float      *r10, *r11, *r12;
	float      *r20, *r21, *r22;
	float      *l00, *l01, *l02;
	float      *l10, *l11, *l12;
	float      *l20, *l21, *l22;
	Locstate   *state = vst->state + offset;

	lambda0 = lambda[0]; lambda1 = lambda[1]; lambda2 = lambda[2];
	r00 = r[0][0]; r01 = r[0][1]; r02 = r[0][2];
	r10 = r[1][0]; r11 = r[1][1]; r12 = r[1][2];
	r20 = r[2][0]; r21 = r[2][1]; r22 = r[2][2];

	l00 = l[0][0]; l01 = l[0][1]; l02 = l[0][2];
	l10 = l[1][0]; l11 = l[1][1]; l12 = l[1][2];
	l20 = l[2][0]; l21 = l[2][1]; l22 = l[2][2];
	
	for (j = start; j < end; ++j)
	{
	    lambda0[j] = v0[j] - c[j];
	    lambda1[j] = v0[j];
	    lambda2[j] = v0[j] + c[j];

	    r00[j] = -rho[j] / c[j];
	    r10[j] = -p[j] / (rho[j]*c[j]);
	    r20[j] = 1.0;

	    r01[j] = GAM[j] / c[j];
	    r11[j] = (GAM[j]*p[j]/rho[j]-c2[j]) / (rho[j]*c[j]);
	    r21[j] = 0.0;

	    r02[j] = rho[j] / c[j];
	    r12[j] = p[j] / (rho[j]*c[j]);
	    r22[j] = 1.0;

	    l00[j] = 0.5 * (GAM[j]*p[j]/rho[j]-c2[j]) / (rho[j]*c[j]);
	    l01[j] = -0.5 * GAM[j] / c[j];
	    l02[j] = 0.5;

	    l10[j] = p[j] / (rho[j]*c[j]);
	    l11[j] = -rho[j] / c[j];
	    l12[j] = 0.0;

	    l20[j] = 0.5 * (c2[j]-GAM[j]*p[j]/rho[j]) / (rho[j]*c[j]);
	    l21[j] = 0.5 * GAM[j] / c[j];
	    l22[j] = 0.5;
	}
	for (j = start; j < end; ++j)
	{
	    if (Params(state[j]) != NULL)
	    {
	        speed = Params(state[j])->avisc.sp_coef*(fabs(v0[j]) + c[j]);
	        if (vmuscl->wave != NULL)
	        {
	            set_max_wave_speed(vmuscl->idir,speed,vst->state[j+offset],
	                               vst->coords[j+offset],vmuscl->wave);
	        }
	        else
	        {
	            int    k, dim = vmuscl->dim;

	            for (k = 0; k < dim; ++k)
	            {
	                set_max_front_speed(k,vmuscl->Q[0][k]*speed,
	                                    vst->state[j+offset],
					    vst->coords[j+offset],
					    vmuscl->front);
	            }
	            set_max_front_speed(dim,speed/vmuscl->dn,
	                                vst->state[j+offset],
					vst->coords[j+offset],
					vmuscl->front);
	        }
	    }
	}

	for (k = 0; k < 3; ++k)
	{
	    sgnlk = sgnl[k]; sgnrk = sgnr[k]; lambdak = lambda[k];
	    for (j = start; j < end; ++j)
	    {
	        if (lambdak[j] > 0.0)
	        {
	            sgnlk[j] = 0.0;
	            sgnrk[j] = 1.0;
	        }
	        else
	        {
	            sgnlk[j] = 1.0;
	            sgnrk[j] = 0.0;
	        }
	    }
	}

#if defined(DEBUG_MUSCL)
	if (debugging("eigen")) print_eigen(start,end,vegn,1);
#endif /* defined(DEBUG_MUSCL) */

	return;
}        /*end g_compute_eigens*/


/*
*			    g_half_step():
*
*	This routine solves the linearized gas dynamics equations in each
*	mesh with the above piecewise linear data for half time step and
*	before interaction.  The resulting left/right states at the edge
*	of each mesh are:
*
*		   uL[][j] = u[][j-i] + duL[][j]
*		   uR[][j] = u[][j] + duR[][j]
*
*		   duL[][j] = (sum over lambda[k][j-1] >= 0)
*			     { 0.5*(1-lambda[k][j-1]*dt/dx)*dq[k][j-1]
*			     + 0.5*dt*sqrt(1+lambda[k][j-1]**2)}*r[][k][j-1]
*
*		   duR[][j] = (sum over lambda[k][j] < 0)
*			     { 0.5*(-1-lambda[k][j]*dt/dx)*dq[k][j]
*			     + 0.5*dt*dqrt(1+lambda[k][j]**2}*r[][k][j]
*
*	These data will be put in the Riemann solver.
*/
/* NOTE: The g_half_step is modified according to
 * "Multidimensioanl Upwind Methods for Hyperbolic Conservation Laws"
 *  By P. Colella, J.C.P. 87 171-200 (1990).
 *  to introduce a reference state.
 */

LOCAL void g_half_step(
	int		start,
	int		end,
	float		dt,
	float		dn,
	Vec_Muscl	*vmuscl)
{
	Vec_Eigen	*vegn;
	int		j, k;
	int		dim, kmax;
	float		**u, **source, **dq, **uL, **uR;
	float		**right_wvs, **left_wvs, **source_sum;
	float		***l, ***r;
	float		**lambda;
	float		**sgnl, **sgnr;
	float		*ssk;
	float		*rk0, *rk1, *rk2;
	float		*ws0, *ws1, *ws2, *wsk;
	float		*lambdak, *lambda1;
	float		*dqk;
	float		*sgnrk, *sgnr1, *sgnlk, *sgnl1;
	float		*uk, *uLk, *uRk;
	float		dtdni = dt/dn;
	static Locstate tmpst = NULL;
	float           lambdaM, lambdaL;
	int             num_comps = 0;

	dim =		vmuscl->dim;
	kmax =		dim+2;
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* Compute differences for mass fraction */
            Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
            if(params->n_comps != 1)
            {
                num_comps = params->n_comps;
                kmax += num_comps;
            }
        }
	u =		vmuscl->u;
	source =	vmuscl->source;
	dq =		vmuscl->dq;
	uL =		vmuscl->uL;
	uR =		vmuscl->uR;
	right_wvs =	vmuscl->right_wvs;
	left_wvs =	vmuscl->left_wvs;
	source_sum =	vmuscl->source_sum;
	vegn =		&vmuscl->Vegn;
	l =		vegn->l;
	r =		vegn->r;
	lambda =	vegn->lambda;
	sgnl =		vegn->sgnl;
	sgnr =		vegn->sgnr;

	/* add up source term contributions for half time-step */

	if ((source != NULL) && (vmuscl->src != NULL))
	{
	    /* New 051005, change argument "source,dim" to "source,dim+num_comps" */
	    left_multiply_state_by_matrix(source_sum,l,source,dim+num_comps,start-1,end);
	    for (k = 0; k < kmax; ++k)
	    {
	    	ssk = source_sum[k];
		dqk = dq[k];
	    	for (j = start - 1; j < end; ++j)
		{
		    if(fabs(dqk[j]) > MACH_EPS)
                        ssk[j] *= 0.5 * dt;
                    else
                        ssk[j] = 0.0;
		}
	    }
	}
	else
	{
	    for (k = 0; k < kmax; ++k)
	    {
	    	ssk = source_sum[k];
	    	for (j = start - 1; j < end; ++j)
	            ssk[j] = 0.0;
	    }
	}

	/* add up contributions from right moving waves */
	for (k = 0; k < 3; ++k)
	{
	    wsk = right_wvs[k];
	    ssk = source_sum[k];
	    lambdak = lambda[k];
	    dqk = dq[k];
	    sgnrk = sgnr[k];
	    for (j = start; j < end; ++j)
	    {
	    	wsk[j] = sgnrk[j-1]*(0.5*(1.0-dtdni*lambdak[j-1])*dqk[j-1] +
	    			ssk[j-1]);
	    }
	}

	for (k = 3; k < kmax; ++k)
	{
	    wsk = right_wvs[k];
	    ssk = source_sum[k];
	    lambda1 = lambda[1];
	    dqk = dq[k];
	    sgnr1 = sgnr[1];
	    for (j = start; j < end; ++j)
	    {
	    	wsk[j] = sgnr1[j-1]*(0.5*(1.0-dtdni*lambda1[j-1])*dqk[j-1] +
					ssk[j-1]);
	    }
	}

	/* compute uL */
	for (k = 0; k < 3; ++k)
	{
	    rk0 = r[k][0];    rk1 = r[k][1];    rk2 = r[k][2];
	    ws0 = right_wvs[0]; ws1 = right_wvs[1]; ws2 = right_wvs[2];
	    uk = u[k];
	    uLk = uL[k];
	    for (j = start; j < end; ++j)
	    {
	    	uLk[j] = uk[j-1] +
			 (ws0[j]*rk0[j-1] + ws1[j]*rk1[j-1] + ws2[j]*rk2[j-1]);
	    }
	}
	for (k = 3; k < kmax; ++k)
	{
	    wsk = right_wvs[k];
	    uk = u[k];
	    uLk = uL[k];
	    for (j = start; j < end; ++j)
	    {
	    	uLk[j] = uk[j-1] + wsk[j];
	    }
	}
	/*Check for negative densities */
	for (j = start; j < end; ++j)
	{
	    Locstate st0 = vmuscl->vst->state[j-1+vmuscl->offset];
	    if (vmuscl->VLst.rho[j] < 0.0)
	    {
		float    min_p = Min_pressure(st0);

		if (tmpst == NULL)
		     g_alloc_state(&tmpst,Params(st0)->sizest); 	

		state_on_adiabat_with_pr(st0,min_p,tmpst,EGAS_STATE);
	        vmuscl->VLst.rho[j] = Dens(tmpst);
	        vmuscl->VLst.e[j] = Energy(tmpst);
		vmuscl->VLst.v[0][j] = u[2][j-1]+riemann_wave_curve(st0,min_p);
	    }
            /* If the reconstruction gives bad state, down to 1st order scheme */
            if(vmuscl->VLst.rho[j]*vmuscl->VLst.e[j] < Min_energy(st0))
            {
                vmuscl->VLst.rho[j] = u[0][j-1];
                vmuscl->VLst.e[j] = u[1][j-1];
                vmuscl->VLst.v[0][j] = u[2][j-1];
            }
	}
	Vec_Gas_field_set(&vmuscl->VLst,rho) = YES;
	Vec_Gas_field_set(&vmuscl->VLst,e) = YES;
	Vec_Gas_field_set(&vmuscl->VLst,v) = YES;

        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* Scaling the mass fraction */
            float       sum_f;
            if(num_comps != 1)
            {
                for (j = start; j < end; ++j)
                {
                    sum_f = 0.0;
                    for(k = 0; k < num_comps; k++)
                        sum_f += vmuscl->VLst.rho0[k][j];
                    for(k = 0; k < num_comps; k++)
                        vmuscl->VLst.rho0[k][j] *= 1.0/sum_f;
                }
                Vec_Gas_field_set(&vmuscl->VLst,rho0) = YES;
            }
        }

	/* add up contributions from left moving waves */
	for (k = 0; k < 3; ++k)
	{
	    wsk = left_wvs[k];
	    ssk = source_sum[k];
	    lambdak = lambda[k];
	    dqk = dq[k];
	    sgnlk = sgnl[k];
	    for (j = start; j < end; ++j)
	    {
	    	wsk[j] = sgnlk[j]*(0.5*(-1.0 - dtdni*lambdak[j])*dqk[j] +
					ssk[j]);
	    }
	}
	for (k = 3; k < kmax; ++k)
	{
	    wsk = left_wvs[k];
	    ssk = source_sum[k];
	    lambda1 = lambda[1];
	    dqk = dq[k];
	    sgnl1 = sgnl[1];
	    for (j = start; j < end; ++j)
	    {
	    	wsk[j] = sgnl1[j]*(0.5*(-1.0 - dtdni*lambda1[j])*dqk[j] +
					ssk[j]);
	    }
	}

	/* compute uR */
	for (k = 0; k < 3; ++k)
	{
	    rk0 = r[k][0];    rk1 = r[k][1];    rk2 = r[k][2];
	    ws0 = left_wvs[0]; ws1 = left_wvs[1]; ws2 = left_wvs[2];
	    uk = u[k];
	    uRk = uR[k];
	    for (j = start; j < end; ++j)
	    {
	    	uRk[j] = uk[j] + (ws0[j]*rk0[j]+ws1[j]*rk1[j]+ws2[j]*rk2[j]);
	    }
	}
	for (k = 3; k < kmax; ++k)
	{
	    wsk = left_wvs[k];
	    uk = u[k];
	    uRk = uR[k];
	    for (j = start; j < end; ++j)
	    {
	    	uRk[j] = uk[j] + wsk[j];
	    }
	}
	/*Check for negative densities */
	for (j = start; j < end; ++j)
	{
	    Locstate st0 = vmuscl->vst->state[j+vmuscl->offset];
	    if (vmuscl->VRst.rho[j] < 0.0)
	    {
		float    min_p = Min_pressure(st0);

		if (tmpst == NULL)
		    g_alloc_state(&tmpst,Params(st0)->sizest); 	

		state_on_adiabat_with_pr(st0,min_p,tmpst,EGAS_STATE);
	        vmuscl->VRst.rho[j] = Dens(tmpst);
	        vmuscl->VRst.e[j] = Energy(tmpst);
		vmuscl->VRst.v[0][j] = u[2][j]-riemann_wave_curve(st0,min_p);
	    }
            /* If the reconstruction gives bad state, down to 1st order scheme */
            if(vmuscl->VRst.rho[j]*vmuscl->VRst.e[j] < Min_energy(st0))
            {
                vmuscl->VRst.rho[j] = u[0][j];
                vmuscl->VRst.e[j] = u[1][j];
                vmuscl->VRst.v[0][j] = u[2][j];
            }
	}
	Vec_Gas_field_set(&vmuscl->VRst,rho) = YES;
	Vec_Gas_field_set(&vmuscl->VRst,e) = YES;
	Vec_Gas_field_set(&vmuscl->VRst,v) = YES;

        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            /* Scaling the mass fraction */
            float       sum_f;
            if(num_comps != 1)
            {
                for (j = start; j < end; ++j)
                {
                    sum_f = 0.0;
                    for(k = 0; k < num_comps; k++)
                        sum_f += vmuscl->VRst.rho0[k][j];
                    for(k = 0; k < num_comps; k++)
                        vmuscl->VRst.rho0[k][j] *= 1.0/sum_f;
                }
                Vec_Gas_field_set(&vmuscl->VRst,rho0) = YES;
            }
        }
#if defined(DEBUG_MUSCL)
	if (debugging("half_step"))
	{
	    int i, j;

	    g_printout_vec_data("Left state on the edge",
				uL[vmuscl->index.density],
				uL[vmuscl->index.energy],
				uL+vmuscl->index.v[0],dim,start,end,"muncons");
	    g_printout_vec_data("Right state on the edge",
				uR[vmuscl->index.density],
				uR[vmuscl->index.energy],
				uR+vmuscl->index.v[0],dim,start,end,"muncons");

	    (void) output();
	    (void) printf("%-4s %-14s %-14s","n","density","int_en");
	    for (i = 0; i < dim; ++i)
		(void) printf(" %-4s%1d%-9s","vel[",i,"]");
	    (void) printf(" Solution before half step\n");
	    for (j = start; j < end; ++j)
	    {
		(void) printf("%-4d",j);
		(void) printf(" %-14g",u[vmuscl->index.density][j]);
		(void) printf(" %-14g",u[vmuscl->index.energy][j]);
		for (k = 0; k < dim; ++k)
		    (void) printf(" %-14g",u[vmuscl->index.v[k]][j]);
		(void) printf("\n");
	    }
	    (void) printf("\n");

	    (void) printf(" Solution after half step\n");
	    (void) output();
	    (void) printf("%-4s %-14s %-14s","n","l_density","l_int_en");
	    for (i = 0; i < dim; ++i)
		    (void) printf(" %-4s%1d%-9s","l_vel[",i,"]");
	    (void) printf(" %-14s %-14s","r_density","r_int_en");
	    for (i = 0; i < dim; ++i)
		    (void) printf(" %-4s%1d%-9s","r_vel[",i,"]");
	    (void) printf("\n");
	    for (j = start; j < end; ++j)
	    {
		(void) printf("%-4d",j);
		(void) printf(" %-14g",uL[vmuscl->index.density][j]);
		(void) printf(" %-14g",uL[vmuscl->index.energy][j]);
		for (k = 0; k < dim; ++k)
			(void) printf(" %-14g",uL[vmuscl->index.v[k]][j]);
		(void) printf(" %-14g",uR[vmuscl->index.density][j]);
		(void) printf(" %-14g",uR[vmuscl->index.energy][j]);
		for (k = 0; k < dim; ++k)
			(void) printf(" %-14g",uR[vmuscl->index.v[k]][j]);
		(void) printf("\n");
	    }
	    (void) printf("\n");
	}
#endif /* defined(DEBUG_MUSCL) */
}		/*end g_half_step*/


LOCAL void g_strong_wave_half_step(
	int		start,
	int		end,
	float		dt,
	float		dn,
	Vec_Muscl	*vmuscl)
{
	bool  redo = NO;
	float *rhol, *rhor;
	float *el, *er;
	float **dq, **uL, **uR;
	int   j, k, kmax;

	g_half_step(start,end,dt,dn,vmuscl);

	dq =		vmuscl->dq;
	uL =		vmuscl->uL;
	uR =		vmuscl->uR;
	kmax =		vmuscl->dim+2;

	rhol = uL[vmuscl->index.density];
	el = uL[vmuscl->index.energy];
	rhor = uR[vmuscl->index.density];
	er = uR[vmuscl->index.energy];
	for (j = start; j < end; ++j)
	{
	    if ((rhol[j] < 0.0) || (el[j] < 0.0))
	    {
	    	for (k = 0; k < kmax; ++k)
	    	    dq[k][j-1] = 0.0;
	    	redo = YES;
	    }
	    if ((rhor[j] < 0.0) || (er[j] < 0.0))
	    {
	    	for (k = 0; k < kmax; ++k)
	    	    dq[k][j] = 0.0;
	    	redo = YES;
	    }
	}
	if (redo == YES)
	    g_half_step(start,end,dt,dn,vmuscl);
}		/*end g_strong_wave_half_step*/


/*ARGSUSED*/
LOCAL void g_first_order_godunov_half_step(
	int		start,
	int		end,
	float		dt,
	float		dn,
	Vec_Muscl	*vmuscl)
{
	int		j, k;
	int		dim, kmax;
	float		**u, **uL, **uR;
	float		*uk, *uLk, *uRk;

	dim =		vmuscl->dim;
	kmax =		dim+2;
	u =		vmuscl->u;
	uL =		vmuscl->uL;
	uR =		vmuscl->uR;

	/* set uL and uR */
	for (k = 0; k < kmax; ++k)
	{
	    uk = u[k]; uLk = uL[k]; uRk = uR[k];
	    for (j = start; j < end; ++j)
	    {
	    	uLk[j] = uk[j - 1];
	    	uRk[j] = uk[j];
	    }
	}

#if defined(DEBUG_MUSCL)
	if (debugging("half_step"))
	{
	    int i, j;

	    g_printout_vec_data("Left state on the edge",
				uL[vmuscl->index.density],
				uL[vmuscl->index.energy],
				uL+vmuscl->index.v[0],dim,start,end,"muncons");
	    g_printout_vec_data("Right state on the edge",
				uR[vmuscl->index.density],
				uR[vmuscl->index.energy],
				uR+vmuscl->index.v[0],dim,start,end,"muncons");

	    (void) output();
	    (void) printf("%-4s %-14s %-14s","n","density","int_en");
	    for (i = 0; i < dim; ++i)
	        (void) printf(" %-4s%1d%-9s","vel[",i,"]");
	    (void) printf(" Solution before half step\n");
	    for (j = start; j < end; ++j)
	    {
		(void) printf("%-4d",j);
		(void) printf(" %-14g",u[vmuscl->index.density][j]);
		(void) printf(" %-14g",u[vmuscl->index.energy][j]);
		for (k = 0; k < dim; ++k)
		    (void) printf(" %-14g",u[vmuscl->index.v[k]][j]);
		(void) printf("\n");
	    }
	    (void) printf("\n");

	    (void) printf(" Solution after half step\n");
	    (void) output();
	    (void) printf("%-4s %-14s %-14s","n","l_density","l_int_en");
	    for (i = 0; i < dim; ++i)
		(void) printf(" %-4s%1d%-9s","l_vel[",i,"]");
	    (void) printf(" %-14s %-14s","r_density","r_int_en");
	    for (i = 0; i < dim; ++i)
		(void) printf(" %-4s%1d%-9s","r_vel[",i,"]");
	    (void) printf("\n");
	    for (j = start; j < end; ++j)
	    {
		(void) printf("%-4d",j);
		(void) printf(" %-14g",uL[vmuscl->index.density][j]);
		(void) printf(" %-14g",uL[vmuscl->index.energy][j]);
		for (k = 0; k < dim; ++k)
		    (void) printf(" %-14g",uL[vmuscl->index.v[k]][j]);
		(void) printf(" %-14g",uR[vmuscl->index.density][j]);
		(void) printf(" %-14g",uR[vmuscl->index.energy][j]);
		for (k = 0; k < dim; ++k)
		    (void) printf(" %-14g",uR[vmuscl->index.v[k]][j]);
		(void) printf("\n");
	    }
	    (void) printf("\n");
	}
#endif /* defined(DEBUG_MUSCL) */
}		/*end g_first_order_godunov_half_step*/

LOCAL void print_eigen(
	int		n1,
	int		n2,
	Vec_Eigen	*vegn,
	int		dim)
{
	int		i, j, k, kmax = dim+2;
	float		**lambda = vegn->lambda;
	float		***l = vegn->l, ***r = vegn->r;

	(void) printf("PRINTOUT OF EIGENVALUES AND EIGENVECTORS\n");
	for (j = n1; j < n2; ++j)
	{
	    (void) printf("Eigenvalues and eigenuni_arrays at index %d\n",j);
	    (void) printf("\tEigenvalues - ");
	    for (k = 0; k < kmax; ++k)
	    	(void) printf("%-14g",lambda[k][j]);
	    (void) printf("\n\tLeft eigenuni_array bi_array, L\n");
	    for (i = 0; i < kmax; ++i)
	    {
	    	(void) printf("\t\t[ ");
	    	for (k = 0; k < kmax; ++k)
	    	    (void) printf("%-14g ",l[i][k][j]);
	    	(void) printf("]\n");
	    }
	    (void) printf("\tRight eigenuni_array bi_array, R\n");
	    for (i = 0; i < kmax; ++i)
	    {
	    	(void) printf("\t\t[ ");
	    	for (k = 0; k < kmax; ++k)
	    	    (void) printf("%-14g ",r[i][k][j]);
	    	(void) printf("]\n");
	    }
	    (void) printf("\tL*R\n");
	    for (i = 0; i < kmax; ++i)
	    {
	    	(void) printf("\t\t[ ");
	    	for (k = 0; k < kmax; ++k)
	    	{
	    	    int	n;
	    	    float	tmp = 0.0;

	    	    for (n = 0; n < kmax; ++n)
	    	    	tmp += l[i][n][j]*r[n][k][j];
	    	    (void) printf("%-14g ",tmp);
	    	}
	    	(void) printf("]\n");
	    }
	    (void) printf("\tR*L\n");
	    for (i = 0; i < kmax; ++i)
	    {
	    	(void) printf("\t\t[ ");
	    	for (k = 0; k < kmax; ++k)
	    	{
	    	    int	n;
	    	    float	tmp = 0.0;

	    	    for (n = 0; n < kmax; ++n)
	    	    	tmp += r[i][n][j]*l[n][k][j];
	    	    (void) printf("%-14g ",tmp);
	    	}
	    	(void) printf("]\n");
	    }
	    (void) printf("End Eigenvalues and eigenuni_arrays ");
	    (void) printf("at index %d\n\n",j);
	}
	(void) printf("END PRINTOUT OF EIGENVALUES AND EIGENVECTORS\n\n\n");
}		/*end print_eigen*/

LOCAL	void print_linear_reconstruction(
	int	  start,
	int	  end,
	Vec_Muscl *vmuscl)
{
	int	dim = vmuscl->dim, i, j, k, kmax;
	int	negn = vmuscl->Vegn.negn;
	float	**dq = vmuscl->dq, *dqk;
	float	*qk;
	float	*qLk, *qRk;
	float	**u = vmuscl->u;
	float	***l = vmuscl->Vegn.l, ***r = vmuscl->Vegn.r;
	static	float ***rm1 = NULL;
	static  float **q = NULL, **qL = NULL, **qR = NULL,
	              **uL = NULL, **uR = NULL;
	static  int   vsize = 0;

	kmax = dim+2;
	if (rm1 == NULL)
	    bi_array(&rm1,negn,negn,sizeof(float*));

	if (vsize < vmuscl->max_vsize)
	{
	    if (q != NULL)
	        free_these(5,q,qL,qR,uL,uR);
	    vsize = vmuscl->max_vsize;
	    bi_array(&q,kmax,vsize,FLOAT);
	    bi_array(&qL,kmax,vsize,FLOAT);
	    bi_array(&qR,kmax,vsize,FLOAT);
	    bi_array(&uL,kmax,vsize,FLOAT);
	    bi_array(&uR,kmax,vsize,FLOAT);
	}

	g_printout_vec_data("Monotonic MUSCL slopes",
			    dq[0],dq[1],dq+2,dim,start,end,"");

	left_multiply_state_by_matrix(q,l,u,dim,start,end);

	for (k = 0; k < kmax; ++k)
	{
	    qk = q[k]; qLk = qL[k]; qRk = qR[k];
	    dqk = dq[k];
	    for (j = start+1; j < end; ++j)
	    {
	    	qLk[j] = qk[j-1] + 0.5 *dqk[j-1];
	    	qRk[j] = qk[j]	 - 0.5 *dqk[j];
	    }
	}

	(void) output();
	(void) printf("%-4s %-14s %-14s","n","density","int_en");
	for (i = 0; i < dim; ++i)
	    (void) printf(" %-4s%1d%-9s","vel[",i,"]");
	(void) printf(" Input solution in physical state coordinates\n");
	for (j = start+1; j < end; ++j)
	{
	    (void) printf("%-4d",j);
	    (void) printf(" %-14g",u[vmuscl->index.density][j-1]);
	    (void) printf(" %-14g",u[vmuscl->index.energy][j-1]);
	    for (k = 0; k < dim; ++k)
	    	(void) printf(" %-14g",u[vmuscl->index.v[k]][j-1]);
	    (void) printf("\n");
	}
	(void) printf("\n");

	(void) output();
	(void) printf("%-4s","n");
	for (k = 0; k < kmax; ++k)
	    (void) printf(" u%-13d",k);
	(void) printf(" Input solution in eigen coordinates\n");
	for (j = start+1; j < end; ++j)
	{
	    (void) printf("%-4d",j);
	    for (k = 0; k < kmax; ++k)
	    	(void) printf(" %-14g",q[k][j-1]);
	    (void) printf("\n");
	    (void) printf("%-4d",j);
	    for (k = 0; k < kmax; ++k)
	    	(void) printf(" %-14g",q[k][j]);
	    (void) printf("\n");
	}
	(void) printf("\n");

	(void) output();
	(void) printf("%-4s","n");
	for (k = 0; k < kmax; ++k)
	    (void) printf(" u%-13d",k);
	(void) printf(" Reconstructed solution in eigen coordinates\n");
	for (j = start+1; j < end; ++j)
	{
	    (void) printf("%-4d",j);
	    for (k = 0; k < kmax; ++k)
	    	(void) printf(" %-14g",qL[k][j]);
	    (void) printf("\n");
	    (void) printf("%-4d",j);
	    for (k = 0; k < kmax; ++k)
	    	(void) printf(" %-14g",qR[k][j]);
	    (void) printf("\n");
	}
	(void) printf("\n");

	for (k = 0; k < negn; ++k)
	    for (j = 0; j < negn; ++j)
	    	rm1[k][j] = r[k][j]-1;
	left_multiply_state_by_matrix(uL,rm1,qL,dim,start+1,end);
	left_multiply_state_by_matrix(uR,r,qR,dim,start+1,end);

	(void) printf(" Reconstructed solution in "
	              "physical state coordinates\n");
	(void) output();
	(void) printf("%-4s %-14s %-14s","n","l_density","l_int_en");
	for (i = 0; i < dim; ++i)
	    (void) printf(" %-4s%1d%-9s","l_vel[",i,"]");
	(void) printf(" %-14s %-14s","r_density","r_int_en");
	for (i = 0; i < dim; ++i)
	    (void) printf(" %-4s%1d%-9s","r_vel[",i,"]");
	(void) printf("\n");
	for (j = start+1; j < end; ++j)
	{
	    (void) printf("%-4d",j);
	    (void) printf(" %-14g",uL[vmuscl->index.density][j]);
	    (void) printf(" %-14g",uL[vmuscl->index.energy][j]);
	    for (k = 0; k < dim; ++k)
	    	(void) printf(" %-14g",uL[vmuscl->index.v[k]][j]);
	    (void) printf(" %-14g",uR[vmuscl->index.density][j]);
	    (void) printf(" %-14g",uR[vmuscl->index.energy][j]);
	    for (k = 0; k < dim; ++k)
	    	(void) printf(" %-14g",uR[vmuscl->index.v[k]][j]);
	    (void) printf("\n");
	}
	(void) printf("\n");
}		/*end print_linear_reconstruction*/

LOCAL	void	limit_strong_wave_in_cell(
	Vec_Muscl	*vmuscl,
	int		start,
	int		end)
{
	int   dim = vmuscl->dim, i, j, k;
	float **du = vmuscl->du;
	float **dq = vmuscl->dq;
	float ***r = vmuscl->Vegn.r;
	float **lambda = vmuscl->Vegn.lambda;
	static float *dq_p,*dq_n;
	static float *du_p,*du_n;

	if (dq_p == NULL)
        {
	    uni_array(&dq_p,3,FLOAT);
	    uni_array(&dq_n,3,FLOAT);
	    uni_array(&du_p,3,FLOAT);
	    uni_array(&du_n,3,FLOAT);
	}
	for (j = start; j < end; ++j)
	{
	    if ((lambda[0][j] > 0.0 && 
	         lambda[1][j] > 0.0 &&
	         lambda[2][j] > 0.0) || 
		(lambda[0][j] <= 0.0 &&
		 lambda[1][j] <= 0.0 && 
		 lambda[2][j] <= 0.0))
		continue;
	    for (k = 0; k < 3; ++k)
	    {
	    	if (lambda[k][j] > 0.0)
		{
		    dq_p[k] = dq[k][j];
		    dq_n[k] = 0.0;
		}
		else
		{
		    dq_p[k] = 0.0;
		    dq_n[k] = dq[k][j];
		}
	    }
	    for (k = 0; k < 3; ++k)
	    {
	    	du_p[k] = dq_p[0]*r[k][0][j] + dq_p[1]*r[k][1][j]
			+ dq_p[2]*r[k][2][j];
	    	du_n[k] = dq_n[0]*r[k][0][j] + dq_n[1]*r[k][1][j]
			+ dq_n[2]*r[k][2][j];
	    }
	    for (k = 0; k < 3; ++k)
	    {
	    	if (fabs(du[k][j]) < fabs(du_p[k]) ||
		    fabs(du[k][j]) < fabs(du_n[k]))
		{
	    	    for (i = 0; i < 3; ++i)
		    	dq[i][j] = 0.0;
		    break;
	    	}
	    }
	}
}	/* end limit_strong_wave_in_cell */

LOCAL	void	left_multiply_state_by_matrix(
	float		**q,
	float		***l,
	float		**u,
	int		dim,
	int		start,
	int		end)
{
	float		*u0, *u1, *u2, *uk;
	float		*lk0, *lk1, *lk2;
	float		*qk;
	int		j, k;
	int		kmax = 2+dim;

	u0 = u[0]; u1 = u[1]; u2 = u[2];
	for (k = 0; k < 3; ++k)
	{
	    lk0 = l[k][0]; lk1 = l[k][1]; lk2 = l[k][2];
	    qk = q[k];
	    for (j = start; j < end; ++j)
	    	qk[j] = lk0[j]*u0[j] + lk1[j]*u1[j] + lk2[j]*u2[j];
	}
	for (k = 3; k < kmax; ++k)
	{
	    qk = q[k]; uk = u[k];
	    for (j = start; j < end; ++j)
	    	qk[j] = uk[j];
	}
}		/*end left_multiply_state_by_matrix*/




/*
*			g_compute_art_visc_coefs();
*/

LOCAL void g_compute_art_visc_coefs(
	int		start,
	int		end,
	Vec_Muscl	*vmuscl)
{
	Vec_Avisc	*avisc = vmuscl->avisc;
	float		*g0;
	Locstate	*state;
	Gas_param	*prms;
	float		**lambda;
	float		dlambda, max_dlambda;
	int		i, j;
	int		negn;

	if (avisc == NULL)
	    return;
	if (avisc->use_lapidus)
	    g_compute_lapidus_art_visc_coefs(start,end,vmuscl);
	else
	{
	    g0 = avisc->g[0];
	    for (j = start; j < end; ++j)
	    	g0[j] = 0.0;
	}
	state = vmuscl->vst->state + vmuscl->offset;
	g0 = avisc->g[0];
	if (avisc->use_linear)
	{
	    float *c = vmuscl->vst->c;
	    float *vn = vmuscl->vst->v[0];
	    prms = NULL;
	    for (j = start; j < end; ++j)
	    {
	        if (Params(state[j]) != prms)
	    	    prms = Params(state[j]);
	        if (prms != NULL)
		{
		    float sp = (fabs(vn[j]) + c[j])*prms->avisc.sp_coef;;
	            g0[j] += 2.0*sp*prms->avisc.linear_visc_coef;
		}
	    }
	}
	if (avisc->use_upwind)
	{
	    lambda = vmuscl->Vegn.lambda;
	    negn = vmuscl->Vegn.negn;
	    prms = NULL;
	    for (j = start; j < end; ++j)
	    {
	        if (Params(state[j]) != prms)
	            prms = Params(state[j]);
	        if (prms == NULL)
	            continue;
	        max_dlambda = 0.0;
	        for (i = 0; i < negn; ++i)
	        {
	            dlambda = lambda[i][j] - lambda[i][j+1];
	            if (dlambda > max_dlambda)
	                max_dlambda = dlambda;
	        }
	        g0[j] += 2.0*prms->avisc.upwind_visc_coef*max_dlambda;
	    }
	}
}		/*end g_compute_art_visc_coefs*/

LOCAL void g_compute_lapidus_art_visc_coefs(
	int		start,
	int		end,
	Vec_Muscl	*vmuscl)
{
	Vec_Avisc	*avisc = vmuscl->avisc;
	Vec_Eigen	*vegn;
	Vec_Gas		*vst;
	int		offset;
	float		**g;
	float		*cs_ave;
	float		*c_ave;
	float		*vn_ave;
	float		**b;
	float		*visc;
	float		*vn;
	float		*c2;
	float		*bi, *b0, *b1, *b2;
	float		*g0, *g1, *g2;
	Locstate	*state;
	Gas_param	*prms;
	float		**lambda, *lambdai;
	int		i, j, negn;

	        /* Read arrays from Vec_Muscl structure */
	offset = vmuscl->offset;
	vst = vmuscl->vst;
	vegn = &vmuscl->Vegn;
	g = avisc->g;
	cs_ave = avisc->cs_ave;
	c_ave = avisc->c_ave;
	vn_ave = avisc->vn_ave;
	b = avisc->b;
	visc = avisc->visc;
	vn = vst->v[0] + offset;
	c2 = vst->c2 + offset;
	state = vst->state + offset;
	lambda = vegn->lambda;
	negn = vegn->negn;

	prms = NULL;
	for (j = start; j < end; ++j)
	{
	    if (Params(state[j]) != prms)
	        prms = Params(state[j]);
	    visc[j] = (prms == NULL) ? 0.0 : prms->avisc.lapidus_visc_coef;
	}

	for (j = start; j < end; ++j)
	{
	    cs_ave[j] = 0.5*(c2[j+1] + c2[j]);
	    vn_ave[j] = 0.5*(vn[j+1] + vn[j]);
	}
	for (j = start; j < end; ++j)
	    c_ave[j] = sqrt(cs_ave[j]);
	
	for (i = 0; i < negn; ++i)
	{
	    bi = b[i];
	    lambdai = lambda[i];
	    for (j = start; j < end; ++j)
	        bi[j] = visc[j]*fabs(lambdai[j+1] - lambdai[j]);
	}

	b0 = b[0]; b1 = b[1]; b2 = b[2];
	g0 = g[0]; g1 = g[1]; g2 = g[2];
	for (j = start; j < end; ++j)
	{
	    g2[j] = 0.5*(b0[j] - 2.0*b1[j] + b2[j])/cs_ave[j];
	    g1[j] = -0.5/cs_ave[j]*(
	        		(2.0*vn_ave[j]+c_ave[j])*b0[j] -
	        		           4.0*vn_ave[j]*b1[j] +
	        		(2.0*vn_ave[j]-c_ave[j])*b2[j]);
	    g0[j] = 0.5/cs_ave[j]*(
	        		vn_ave[j]*(vn_ave[j]+c_ave[j])*b0[j] -
	        		2.0*(sqr(vn_ave[j])-cs_ave[j])*b1[j] +
	        		vn_ave[j]*(vn_ave[j]-c_ave[j])*b2[j]);
	}
#if defined(DEBUG_MUSCL)
	if (debugging("art_visc"))
	{
	    (void) printf("Artificial viscosity coefficients\n");
	    (void) printf("%-4s %-14s %-14s %-14s\n","n","g0","g1","g2");
	    for (j = start; j < end; ++j)
	        (void) printf("%-4d %-14g %-14g %-14g\n",j,g0[j],g1[j],g2[j]);
	    (void) printf("\nEnd Artificial viscosity coefficients\n");
	    (void) printf("Wave speed gradients\n");
	    (void) printf("%-4s %-14s %-14s %-14s\n","n","b0","b1","b2");
	    for (j = start; j < end; ++j)
	        (void) printf("%-4d %-14g %-14g %-14g\n",j,b0[j],b1[j],b2[j]);
	    (void) printf("End Wave speed gradients\n");
	    (void) printf("Mid wave speeds\n");
	    (void) printf("%-4s %-14s %-14s %-14s\n","n","v-c","v","v+c");
	    for (j = start; j < end; ++j)
	    {
	        (void) printf("%-4d %-14g %-14g %-14g\n",j,vn_ave[j]-c_ave[j],
	        	      vn_ave[j],vn_ave[j]+c_ave[j]);
	    }
	    (void) printf("End Mid wave speeds\n");
	    (void) printf("Interpolation test\n");
	    (void) printf("%-4s %-14s %-14s %-14s\n",
	        	  "n","G(v-c) - b0","G(v) - b1","G(v+c) - b2");
	    for (j = start; j < end; ++j)
	    {
	        float	vmc, v, vpc;
	        float	G0, G1, G2;

	        vmc = vn_ave[j] - c_ave[j];
	        v = vn_ave[j];
	        vpc = vn_ave[j] + c_ave[j];
	        G0 = g0[j] + g1[j]*vmc + g2[j]*vmc*vmc;
	        G1 = g0[j] + g1[j]*v + g2[j]*v*v;
	        G2 = g0[j] + g1[j]*vpc + g2[j]*vpc*vpc;
	        (void) printf("%-4d %-14g %-14g %-14g\n",j,G0-b0[j],
	    		  G1-b1[j],G2-b2[j]);
	    }
	    (void) printf("End Interpolation test\n");
	}
#endif /* defined(DEBUG_MUSCL) */
}		/*end g_compute_lapidus_art_visc_coefs*/


/*ARGSUSED*/
LOCAL	void g_muscl_lin_approx_rsolver(
	int	   start,
	int	   end,
	float      **uL,
	Vec_Gas    *vlst,
	float      **uR,
	Vec_Gas    *vrst,
	float      **uM,
	Vec_Gas    *vmst,
	MUSCL_FLUX *Flux,
	Vec_Muscl  *vmuscl)
{
	g_lin_rsoln(start,end,uL,uR,uM,vmuscl);

        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            Gas_param *params = Params(vmuscl->vst->state[vmuscl->offset]);
            int    num_comps,i,j;
            if((num_comps = params->n_comps) != 1)
            {
                for (j = start; j < end; ++j)
                {
                    if(vmst->v[0][j] > 0.0)
                    {
                        for(i = 0; i < num_comps; i++)
                            vmst->rho0[i][j] = uL[vmuscl->index.prho[i]][j];
                    }
                    else
                    {
                        for(i = 0; i < num_comps; i++)
                            vmst->rho0[i][j] = uR[vmuscl->index.prho[i]][j];
                    }
                }
                Vec_Gas_field_set(vmst,rho0) = YES;
            }
        }

	Vec_Gas_field_set(vmst,rho) = YES;
	Vec_Gas_field_set(vmst,e) = YES;
	Vec_Gas_field_set(vmst,v) = YES;
	g_load_muscl_flux(start,end,uM,vmst,Flux,vmuscl);
}		/*end g_muscl_lin_approx_rsolver*/


/*
*		      g_lin_rsoln():
*
*	This routine is the approximate Riemann solver. It solves the
*	linearized system. It does not handle sonic rarefactions properly.
*/

LOCAL void g_lin_rsoln(
	int	  start,
	int	  end,
	float     **uL,
	float     **uR,
	float     **uM,
	Vec_Muscl *vmuscl)
{
	float		**awv;
	Vec_Eigen	*vegn;
	int		j, k;
	int		kmax;
	float		***l, ***r;
	float		**lambda;
	float		*lambdak;
	float		**sgnm;
	float		*awv0, *awv1, *awv2, *awvk;
	float		*sgnmk, *sgnm1;
	float		*uR0, *uR1, *uR2, *uRk;
	float		*uL0, *uL1, *uL2, *uLk;
	float		*uMk;
	float		*lk0, *lk1, *lk2;
	float		*rk0, *rk1, *rk2;

	kmax =		vmuscl->dim+2;
	vegn =		&vmuscl->Vegn;
	l =		vegn->l;
	r =		vegn->r;
	sgnm =		vegn->sgnm;
	lambda =	vegn->lambda;
	awv =		vmuscl->awv;

	/* Do left contributions */
	for (k = 0; k < 3; ++k)
	{
	    sgnmk = sgnm[k];
	    lambdak = lambda[k];
	    for (j = start; j < end; ++j)
	    {
	    	sgnmk[j] = (lambdak[j-1] < 0.0) ? 1.0 : 0.0;
	    }
	}

	for (k = 0; k < 3; ++k)
	{
	    awvk = awv[k];
	    sgnmk = sgnm[k];
	    uR0 = uR[0]; uR1 = uR[1]; uR2 = uR[2];
	    uL0 = uL[0]; uL1 = uL[1]; uL2 = uL[2];
	    lk0 = l[k][0]; lk1 = l[k][1]; lk2 = l[k][2];
	    for (j = start; j < end; ++j)
	    {
	    	awvk[j] = sgnmk[j]*( (uR0[j] - uL0[j])*lk0[j-1] +
	    		             (uR1[j] - uL1[j])*lk1[j-1] +
	    		             (uR2[j] - uL2[j])*lk2[j-1]);
	    }
	}
	for (k = 3; k < kmax; ++k)
	{
	    awvk = awv[k];
	    sgnm1 = sgnm[1];
	    uRk = uR[k];
	    uLk = uL[k];
	    for (j = start; j < end; ++j)
	    {
	    	awvk[j] = sgnm1[j]*(uRk[j] - uLk[j]);
	    }
	}

	for (k = 0; k < 3; ++k)
	{
	    uMk = uM[k];
	    uLk = uL[k];
	    uRk = uR[k];
	    rk0 = r[k][0]; rk1 = r[k][1]; rk2 = r[k][2];
	    awv0 = awv[0]; awv1 = awv[1]; awv2 = awv[2];
	    for (j = start; j < end; ++j)
	    {
	    	uMk[j] = 0.5*(uLk[j] + ( awv0[j]*rk0[j-1] +
	        			 awv1[j]*rk1[j-1] +
	        			 awv2[j]*rk2[j-1]));
	    }
	}
	for (k = 3; k < kmax; ++k)
	{
	    uLk = uL[k];
	    uRk = uR[k];
	    uMk = uM[k];
	    awvk = awv[k];
	    for (j = start; j < end; ++j)
	    {
	    	uMk[j] = 0.5*(uLk[j] + awvk[j]);
	    }
	}

	/* Do right contributions */
	for (k = 0; k < 3; ++k)
	{
	    sgnmk = sgnm[k]; lambdak = lambda[k];
	    for (j = start; j < end; ++j)
	    {
	    	sgnmk[j] = (lambdak[j] > 0.0) ? -1.0 : 0.0;
	    }
	}

	for (k = 0; k < 3; ++k)
	{
	    awvk = awv[k];
	    sgnmk = sgnm[k];
	    uR0 = uR[0]; uR1 = uR[1]; uR2 = uR[2];
	    uL0 = uL[0]; uL1 = uL[1]; uL2 = uL[2];
	    lk0 = l[k][0]; lk1 = l[k][1]; lk2 = l[k][2];
	    for (j = start; j < end; ++j)
	    {
	    	awvk[j] = sgnmk[j]*( (uR0[j] - uL0[j])*lk0[j] +
	        		     (uR1[j] - uL1[j])*lk1[j] +
	        		     (uR2[j] - uL2[j])*lk2[j]);
	    }
	}
	for (k = 3; k < kmax; ++k)
	{
	    awvk = awv[k];
	    sgnm1 = sgnm[1];
	    uRk = uR[k];
	    uLk = uL[k];
	    for (j = start; j < end; ++j)
	    {
	    	awvk[j] = sgnm1[j]*(uRk[j] - uLk[j]);
	    }
	}

	for (k = 0; k < 3; ++k)
	{
	    uMk = uM[k];
	    uLk = uL[k];
	    uRk = uR[k];
	    rk0 = r[k][0]; rk1 = r[k][1]; rk2 = r[k][2];
	    awv0 = awv[0]; awv1 = awv[1]; awv2 = awv[2];
	    for (j = start; j < end; ++j)
	    {
	    	uMk[j] += 0.5*(uRk[j] + ( awv0[j]*rk0[j] +
	    			          awv1[j]*rk1[j] +
	    			          awv2[j]*rk2[j]));
	    }
	}
	for (k = 3; k < kmax; ++k)
	{
	    uLk = uL[k];
	    uRk = uR[k];
	    uMk = uM[k];
	    awvk = awv[k];
	    for (j = start; j < end; ++j)
	    {
	    	uMk[j] += 0.5*(uRk[j] + awvk[j]);
	    }
	}
}		/*end g_lin_rsoln*/

/*
*			g_add_art_visc1();
*
*	This function adds the first half of artificial viscosity.  It is
*	added to the mid states which will be used to compute conservative
*	fluxes and then in the conservative difference scheme.
*/

LOCAL void g_add_art_visc1(
	int		start,
	int		end,
	float           **uM,
	Vec_Muscl	*vmuscl)
{
	Vec_Avisc	*avisc;
	int		dim;
	int		j, k;
	float		*g2, *g1;
	float		**u, **ucon;
	float		*p, **F, **g, **uconM;
	float		*uconMk, *uconk, *Fk;
	float		du;
	int		id = vmuscl->index.density;
	int		ie = vmuscl->index.energy;

	if ((avisc = vmuscl->avisc) == NULL)
	    return;
	g = avisc->g;
	g2 = g[2], g1 = g[1];
	if (g2 == NULL || g1 == NULL)
	    return;

	u =		vmuscl->u;
	dim =		vmuscl->dim;
	ucon =		vmuscl->ucon;
	F =		vmuscl->Flux.F;
	p =		vmuscl->vst->p + vmuscl->offset;
	uconM =		avisc->uconM;

	flux_vectors(start,end,u,p,&vmuscl->Flux,vmuscl);
	g_uncons_to_cons(start,end,dim,uM,uconM);

	for (k = 0; k < dim+2; ++k)
	{
	    uconk = ucon[k];
	    uconMk = uconM[k];
	    Fk = F[k];
	    for (j = start; j < end; ++j)
	    {
		du = 0.5*(g2[j-1]*(Fk[j] - Fk[j-1]) +
		                  g1[j-1]*(uconk[j] - uconk[j-1]));
	    	if ((k == id || k == ie) && du >= uconMk[j])
                {
                    const char *name = (k == id) ? "density" : "energy";
                    screen("ERROR in g_add_art_visc1():\n");
                    screen("nonlinead artificial viscosity causes middle ");
                    screen("state negative %s,\n",name);
                    screen("recommend to use default (no viscosity) or ");
                    screen("linear artificial viscosity\n");
                    clean_up(ERROR);
                }
		uconMk[j] -= du;
	    }
	}
	g_cons_to_uncons(start,end,dim,uconM,uM);
}		/*end g_add_art_visc1*/

/*
*			g_add_art_visc2();
*
*	This function adds the second half of the artificial viscosity.	 All
*	we need do here is to modify appropriately the flux uni_arrays already
*	computed.  Thus the conservative differencing algorithm does not 
*	need to be modified.
*/

LOCAL void g_add_art_visc2(
	int		start,
	int		end,
	Vec_Muscl	*vmuscl)
{
	Vec_Avisc	*avisc;
	int		dim, j, k;
	float		*g0;
	float		*uconk, *Fk;
	float		**ucon, **F;

	if ((avisc = vmuscl->avisc) == NULL)
	    return;
	dim = vmuscl->dim;
	ucon = vmuscl->ucon;
	F = vmuscl->Flux.F;
	g0 = avisc->g[0];

	for (k = 0; k < dim+2; ++k)
	{
	    uconk = ucon[k];
	    Fk = F[k];
	    for (j = start; j < end; ++j)
	    	Fk[j] -= 0.5*g0[j-1]*(uconk[j] - uconk[j-1]);
	}
}		/*end g_add_art_visc2*/


LOCAL void g_uncons_to_cons(
	int		start,
	int		end,
	int		dim,
	float		**uncons,
	float		**cons)
{
	int		j;
	float		*rho = uncons[0],	*cons0 = cons[0];
	float		*e = uncons[1],		*en_den = cons[1];
	float		*v0 = uncons[2],	*m0 = cons[2];
	float		*v1, *v2, *m1, *m2, v_norm;

	switch (dim)
	{
	case 1:
	    for (j = start; j < end; ++j)
	    {
	    	m0[j] = v0[j]*rho[j];
	    	v_norm = sqr(v0[j]);
	    	en_den[j] = (e[j] + 0.5*v_norm)*rho[j];
	    	cons0[j] = rho[j];
	    }
	    break;
	case 2:
	    v1 = uncons[3]; m1 = cons[3];
	    for (j = start; j < end; ++j)
	    {
	    	m1[j] = v1[j]*rho[j];
	    	m0[j] = v0[j]*rho[j];
	    	v_norm = sqr(v1[j]) + sqr(v0[j]);
	    	en_den[j] = (e[j] + 0.5*v_norm)*rho[j];
	    	cons0[j] = rho[j];
	    }
	    break;
	case 3:
	    v1 = uncons[3]; m1 = cons[3];
	    v2 = uncons[4]; m2 = cons[4];
	    for (j = start; j < end; ++j)
	    {
	    	m2[j] = v2[j]*rho[j];
	    	m1[j] = v1[j]*rho[j];
	    	m0[j] = v0[j]*rho[j];
	    	v_norm = sqr(v2[j]) + sqr(v1[j]) + sqr(v0[j]);
	    	en_den[j] = (e[j] + 0.5*v_norm)*rho[j];
	        cons0[j] = rho[j];
	    }
	    break;
	}
}		/*end g_uncons_to_cons*/

LOCAL void g_cons_to_uncons(
	int		start,
	int		end,
	int		dim,
	float		**cons,
	float		**uncons)
{
	int		j;
	float		*rho = cons[0],		*uncons0 = uncons[0];
	float		*en_den = cons[1],	*e = uncons[1];
	float		*m0 = cons[2],		*v0 = uncons[2];
	float		*m1, *m2, *v1, *v2, v_norm;

	switch (dim)
	{
	case 1:
	    for (j = start; j < end; ++j)
	    {
	    	v0[j] = m0[j]/rho[j];
	    	v_norm = sqr(v0[j]);
	    	e[j] = en_den[j]/rho[j] - 0.5*v_norm;
	    	uncons0[j] = rho[j];
	    }
	    break;
	case 2:
	    m1 = cons[3]; v1 = uncons[3];
	    for (j = start; j < end; ++j)
	    {
	    	v1[j] = m1[j]/rho[j];
	    	v0[j] = m0[j]/rho[j];
	    	v_norm = sqr(v0[j]) + sqr(v1[j]);
	    	e[j] = en_den[j]/rho[j] - 0.5*v_norm;
	    	uncons0[j] = rho[j];
	    }
	    break;
	case 3:
	    m1 = cons[3]; v1 = uncons[3];
	    m2 = cons[4]; v2 = uncons[4];
	    for (j = start; j < end; ++j)
	    {
	    	v2[j] = m2[j]/rho[j];
	    	v1[j] = m1[j]/rho[j];
	    	v0[j] = m0[j]/rho[j];
	    	v_norm = sqr(v0[j]) + sqr(v1[j]) + sqr(v2[j]);
	    	e[j] = en_den[j]/rho[j] - 0.5*v_norm;
	    	uncons0[j] = rho[j];
	    }
	    break;
	}
}		/*end g_cons_to_uncons*/


#if defined(DEBUG_MUSCL)

/*
*		EXPLANATION OF DEBUGGING SWITCHES
*
*	The printout of the interior state debugging is so extensive
*	that it is impractical in most runs to printout unconditional
*	debug lines.  Therefore the following control strings for
*	debugging have been added to allow the restriction of debugging
*	to only certain areas of the program.  In particular it is possible
*	to restrict the printout as a function of the grid cell index
*	being computed.  The use of these debugging control strings is
*	enabled by turning on the debugging string MUSCL_SOLVER.
*
*		CONTROLLING DEBUGGING BY SWEEP DIRECTION
*
*	Key debug strings:
*		xonly	yonly	zonly
*
*	If any of these debug strings is turned on,  then the debugging
*	will be restricted to that corresponding sweep direction.  These
*	values are obviously mutually exclusive.
*
*		CONTROLLING DEBUGGING BY SWEEP TYPE
*
*	Key debug strings:
*		mtsten_all	mtsten_some	mtsten_none
*		mirreg_all	mirreg_some	mirreg_none
*		mreg_all	mreg_some	mreg_none
*
*	The debug strings are formed in two parts by combining the
*	prefixes mtsten,  mirreg,  and mreg,  with the suffixes
*	_all, _some, or _none.  The suffix arguments are mutually
*	exclusive with the least restrictive flags overriding the
*	more restrictive flags.  (For example requesting all and some 
*	yields all.)
*	
*	DEBUG STRING PREFIX	EXPLANATION
*
*	mtsten			Controls debugging for the tangential sweep.
*
*	mirreg			Controls debugging for the irregular stencil
*				sweeps.
*
*	mreg			Controls debugging for the regular stencil
*				sweep.
*
*	DEBUG STRING SUFFIX	EXPLANATION
*
*	_all			Turns on all debugging strings for MUSCL code.
*				The printout of these strings will be 
*				controled by the cell index and time step
*				control desribed below.
*
*	_some			Allows controled output of MUSCL code
*				debugging.  If this string is active you
*				can turn on specific debug strings
*				in the MUSCL code (eg cg_rsolve,  etc)
*				and the corresponding debug lines will
*				be controled by the cell index and
*				time step control described below.
*
*	_none			Turns off all debugging for MUSCL code
*				tangential sweep. In this case no debugging
*				output will be printed whether the individual
*				debug strings are turned on or not.
*
*
*		CONTROLLING DEBUGGING BY CELL INDEX
*
*	Key debug strings:
*		mrgrid
*		LX (followed by an integer not separated by a space, eg LX10)
*		LY (followed by an integer not separated by a space, eg LY10)
*		LZ (followed by an integer not separated by a space, eg LZ10)
*		UX (followed by an integer not separated by a space, eg UX10)
*		UY (followed by an integer not separated by a space, eg UY10)
*		UZ (followed by an integer not separated by a space, eg UZ10)
*
*	The grid cell index restriction is communicated by multiple debug 
*	strings.  If the debug string mrgrid is input this should be followed
*	by debug strings of the form LXnn, UXNN,  LYnn,  UYNN,  LZnn, UZNN,
*	where the nn and NN (which need not be all equal) are the indicies of
*	the cells for which debugging is desired to be turned on.  If any of
*	these are omited,  they default to the ends of the computational grid.
*
*			CONTROLLING DEBUGGING BY TIME STEP
*
*	Key debug strings:
*		mrtstep
*		mt (followed by an integer not separated by a space, eg mt10)
*		MT (followed by an integer not separated by a space, eg MT10)
*
*	The time step restriction is also communicated by multiple debug
*	strings.  If the debug string mrtstep is input this should be followed
*	by the debug strings mtnn and MTNN,  where nn is the first time to
*	print debuggin and NN is the last.
*
*	Example:	The following collection of debug strings
*	will print only output during the irregular sweep for cells with
*	indicies ix = 28 and iy = 50 for times 0 and 1.
*
*	: MUSCL_SOLVER
*	: mreg_none
*	: mtsten_none
*	: mirreg_all
*	: mrgrid
*	: LY28
*	: UY28
*	: LX50
*	: UX50
*	: mrtstep
*	: mt0
*	: MT1
*/


LOCAL	const char **set_muscl_debugging(
	int		sn,
	int		*ip,
	int		*ic,
	Tan_stencil	*tsten,
	Front		*fr,
	int		vsize)
{
	const char	**debug_strings = NULL;
	int		dim;

	if (!debugging("MUSCL_SOLVER"))
	    return NULL;
	if (tsten != NULL)
	{
	    if (debugging("mtsten_all"))
	    {
	    	debug_strings = toggle_muscl_debugging("ALL",sn,ip,ic,fr,NO);
	    }
	    else if (debugging("mtsten_some"))
	    {
	    	debug_strings = toggle_muscl_debugging("SOME",sn,ip,ic,fr,NO);
	    }
	    else if (debugging("mtsten_none"))
	    {
	    	debug_strings = toggle_muscl_debugging("NONE",sn,ip,ic,fr,NO);
	    }
	    if (debugging("oned_MUSCL"))
	    {
	    	(void) printf("oned_MUSCL called from oblique solver.\n");
	    	dim = fr->rect_grid->dim;
	    	print_general_vector("At coords ",Coords(tsten->p[0]),dim,
	        		     "\n\n");
	    }
	}
	else if (vsize == 5)
	{
	    if (debugging("mirreg_all"))
	    {
	    	debug_strings = toggle_muscl_debugging("ALL",sn,ip,ic,fr,NO);
	    }
	    else if (debugging("mirreg_some"))
	    {
	    	debug_strings = toggle_muscl_debugging("SOME",sn,ip,ic,fr,NO);
	    }
	    else if (debugging("mirreg_none"))
	    {
	    	debug_strings = toggle_muscl_debugging("NONE",sn,ip,ic,fr,NO);
	    }
	    if (debugging("oned_MUSCL"))
	    {
	    	dim = fr->rect_grid->dim;
	    	(void) printf("oned_MUSCL called from irregular solver.\n");
	    	(void) printf("Sweep direction = %d, ",ip[sn]);
	    	print_int_vector("at ic ",ic,dim,"\n");
	    	(void) printf("\n");
	    }
	}
	else 
	{
	    if (debugging("mreg_all"))
	    {
	    	debug_strings = toggle_muscl_debugging("ALL",sn,ip,ic,fr,YES);
	    }
	    else if (debugging("mreg_some"))
	    {
	    	debug_strings = toggle_muscl_debugging("SOME",sn,ip,ic,fr,YES);
	    }
	    else if (debugging("mreg_none"))
	    {
	    	debug_strings = toggle_muscl_debugging("NONE",sn,ip,ic,fr,YES);
	    }
	    if (debugging("oned_MUSCL"))
	    {
	    	dim = fr->rect_grid->dim;
	    	(void) printf("oned_MUSCL called from vector solver.\n");
	    	(void) printf("Sweep direction = %d, ",ip[sn]);
	    	print_int_vector("at ic ",ic,dim,"\n");
	    	(void) printf("\n");
	    }
	}
	return debug_strings;
}		/*end set_muscl_debugging*/

LOCAL	const char **toggle_muscl_debugging(
	const char	*howmuch,
	int		swp_num,
	int		*iperm,
	int		*icoords,
	Front		*fr,
	int		reg_sweep)
{
	int		  i, j, idir;
	int		  dim = fr->rect_grid->dim;
	int		  ndb;
	long		  int mt, MT;
	long		  int mg, MG;
	char		  **db;
	const char        **c, **s;
	static const char *debug_strings[100];
	static const char        *swpdb[3] = {"xonly", "yonly", "zonly"};
	static const char        *mnames[3] = {"LX", "LY", "LZ"};
	static const char        *Mnames[3] = {"UX", "UY", "UZ"};
	static const char        *muscl_debug_strings[] = { "cg_rsolve",
	        					    "MUSCL",
	        					    "oned_MUSCL",
	        					    "MUSCLob",
	        					    "lriem",
	        					    "vflux",
	        					    "load_Vgas",
	        					    "eigen",
	        					    "lcnst",
	        					    "half_step",
	        					    "csrc",
	        					    "vgdnv",
	        					    "vgdnv_ie",
	        					    "art_visc",
	        					    NULL};

	if (iperm != NULL)
	{
	    idir = iperm[swp_num];
	    if (!debugging(swpdb[idir]))
	    {
	    	for (i = 0; i < dim; ++i)
	    	{
	    	    if (i == idir)
	        	continue;
	            if (debugging(swpdb[i]))
	        	howmuch = "NONE";
	        }
	    }
	    if (debugging("mrgrid"))
	    {
	    	db = debugging_names(&ndb);
	    	for (i = 0; i < dim; ++i)
	    	{
	    	    if (reg_sweep == YES && i == idir)
	        	continue;
	    	    mg = 0;
	    	    MG = INT_MAX;
	    	    for (j = 0; j < ndb; ++j)
	    	    {
	    	        if (strncmp(db[j],mnames[i],2) == 0)
	    	            mg = atoi(db[j]+2);
	    	        else if (strncmp(db[j],Mnames[i],2) == 0)
	    	    	    MG = atoi(db[j]+2);
	    	    }
	    	    if (icoords[i] < mg || icoords[i] > MG)
	    	    	howmuch = "NONE";
	    	}
	    }
	}

	if (debugging("mrtstep"))
	{
	    mt = 0;
	    MT = INT_MAX;
	    db = debugging_names(&ndb);
	    for (i = 0; i < ndb; ++i)
	    {
	    	if (strncmp(db[i],"mt",2) == 0)
	    	    mt = atoi(db[i]+2);
	    	else if (strncmp(db[i],"MT",2) == 0)
	    	    MT = atoi(db[i]+2);
	    }
	    if (fr->step < mt || fr->step > MT)
	    	howmuch = "NONE";
	}

	c = debug_strings;
	if (strcmp(howmuch,"ALL") == 0)
	{
	    for (s = muscl_debug_strings; *s != NULL; ++s)
	    {
	    	if (debugging(*s))
	    	    continue;
	    	add_to_debug(*s);
	    	*c++ = "TURNED_ON";
	    	*c++ = *s;
	    }
	}
	if (strcmp(howmuch,"NONE") == 0)
	{
	    for (s = muscl_debug_strings; *s != NULL; ++s)
	    {
	    	if (!debugging(*s))
	            continue;
	    	remove_from_debug(*s);
	    	*c++ = "TURNED_OFF";
	    	*c++ = *s;
	    }
	}
	*c = NULL;
	return debug_strings;
}		/*end toggle_muscl_debugging*/

LOCAL	void reset_muscl_debugging(
	const char **debug_strings)
{
	const char **s;

	if (debug_strings == NULL)
	    return;

	for (s = debug_strings; *s != NULL; ++s)
	{
	    if (strcmp(*s++,"TURNED_OFF") == 0)
	    	add_to_debug(*s);
	    else
	    	remove_from_debug(*s);
	}
}		/*end reset_muscl_debugging*/

LOCAL void oned_turbulence_boundarylayer_solver(
		int dim,
		float u, 
                float v,
		float w, 
		float y0,
		float yw,
		float *flux,
                Locstate statep)	
{      
		float Tmatch = temperature(statep);
		debug_print("WALL_MODEL","Entered oned_turbulence_boundarylayer_solver()\n");		
  		float delta = fabs(y0 - yw);
		float ymatch = delta;
  		float *Urans, *Trans, *Rrans, *MUrans, *MUTrans, *Krans, *aw, *ap, *ae, *qp;
		float *dyrans, *yfrans;
		float ycrans[200];
  		float wallunit= 0.001;
		
		float dyrans0;
		int NY = 1, j, i, k, kk ,iter;
		float rnm1 = 1;
		float lvmin = 1000;
		float Umatch;
		float Tw;
		float muw,kw;
		float tauwold, qwold, converged, damping;
		float intfac,muint,mutint,kint,factor,unit;
		float err_Vel,err_Temp;
		float tauw, rhow, utau, lv, qw;
		float *pmu, *pk;
		float uint;
		float rho0;
		float sigma13, sigma23, sigma33, usigma, q;;
	
		converged = HUGE_VAL;	

		//  Set up auxiliary quantities in all cells in 1D grid
		//  The first grid point is at about y^+ = 1. 
		wallunit = 0.0001;//FIXME
                dyrans0 = wallunit;
		ycrans[0] = 0.5;
  		while( ycrans[NY-1] < ymatch/dyrans0)
		{
			ycrans[NY] = ycrans[NY-1] + 0.5*rnm1*(1.0 + WALLMODEL_STRETCHING);
        		rnm1 *= WALLMODEL_STRETCHING;
       			NY++;
  		}
		dyrans0 = ymatch/ycrans[NY-1];

		uni_array(&dyrans,NY+4,sizeof(float));
		uni_array(&yfrans,NY+4,sizeof(float));
		
		ycrans[0] *= dyrans0;
		dyrans[0] = dyrans0;
		yfrans[0] = 0.0;
		for(j=1;j < NY; j++){
			ycrans[j] *= dyrans0;
			dyrans[j] = dyrans[j-1]*WALLMODEL_STRETCHING;
			yfrans[j] = yfrans[j-1] + dyrans[j-1];
		}


		//if(fabs(ycrans[NY-1] - ymatch) > 0.00001)
		//	printf( "Error in viscousWallY : ycrans[NY-1] does not match ymathc\n");

		uni_array(&Urans,NY+4,sizeof(float));
		uni_array(&Trans,NY+4,sizeof(float));
		uni_array(&Rrans,NY+4,sizeof(float));
		uni_array(&MUrans,NY+4,sizeof(float));
		uni_array(&MUTrans,NY+4,sizeof(float));
		uni_array(&Krans,NY+4,sizeof(float));
		uni_array(&aw,NY+4,sizeof(float));
		uni_array(&ap,NY+4,sizeof(float));
		uni_array(&ae,NY+4,sizeof(float));
		uni_array(&qp,NY+4,sizeof(float));
		
	        rho0 = Dens(statep);
		float C_P = C_P(statep);

	  	//The actual wall-model --- solve ODEs for each wall face
	  	Umatch = sqrt(u*u + v*v);
		if(Tmatch >= Wall_Temp)
	  	Tw = Wall_Temp;
		else Tw = Tmatch;
	  	//Initial Guess -- linear profiles	  	
	  	for (j=0; j < NY; j++)
		{
			Urans[j] = Umatch/ymatch*ycrans[j];
		  	Trans[j] = Tw + (Tmatch - Tw)*ycrans[j]/ymatch;
	 	}
		dynamic_viscosity_thermalconduct(statep,Tw,&muw,&kw);
	  	tauw = muw*Urans[0]/ycrans[0];
	  	rhow = rho0*Tmatch/Tw;
	  	utau = sqrt(tauw/rhow);
	  	lv = muw/(rhow*utau);
	  	qw = - kw*(Trans[0] - Tw)/ycrans[0];

		//Iterate until convergence
		iter = 0;
	 	do{
			
			//  Set up auxiliary quantities in all cells in 1D grid
			tauwold = tauw ;
	 		qwold = qw ;
		
	  		for (j=0 ; j < NY ; j++) 
			{
				dynamic_viscosity_thermalconduct(statep,Trans[j],&MUrans[j],&Krans[j]);
	    			damping = 1.0 - exp( - ycrans[j] / ( lv * WALLMODEL_APLUS)) ;
				damping = damping*damping;
	    			Rrans[j] =  rho0*Tmatch/Trans[j];
	    			MUTrans[j] = WALLMODEL_KAPPA * ycrans[j] * sqrt( Rrans[j] * tauw ) * damping ;
	    			ap[j] = 0.0 ;
	    			qp[j] = 0.0 ;   //  if LHS not zero, then should set qp to LHS * dyrans[j]...
	  		}

	  		//  Set up momentum equation using finite volume
			/*
			intfac = 1.0/(1.0+WALLMODEL_STRETCHING);
			muint = (1.0 - intfac)*muw + intfac*MUrans[0];
			kint = (1.0 - intfac)*kw + intfac*Krans[0];
			mutint = intfac*MUTrans[0];
			factor = (muint + mutint)/ycrans[0];
			ap[0] = - factor;
			*/	
			ap[0] = -muw/ycrans[0];
	 		for (j=1 ; j < NY ; j++) 
			{
            			intfac = ( yfrans[j] - ycrans[j-1] ) / ( ycrans[j] - ycrans[j-1] ) ;
             			muint  = ( 1.0 - intfac ) * MUrans[j-1]  + intfac * MUrans[j] ;
	     			kint = (1.0 - intfac ) * Krans[j-1] + intfac*Krans[j];
             			mutint = ( 1.0 - intfac ) * MUTrans[j-1] + intfac * MUTrans[j] ;
            			factor = ( muint + mutint ) / ( ycrans[j] - ycrans[j-1] ) ;
            			aw[j] = factor;
            			ap[j] -= factor ;
            			ae[j-1] = factor ;
            			ap[j-1] -= factor ;
          		}
          		Urans[NY-1] = Umatch ;

	  		//  Solve momentum equation by single TDMA sweep
	  		for (j=1 ; j < NY-1 ; j++) 
			{
	    			factor = aw[j] / ap[j-1] ;
	    			ap[j] -= factor * ae[j-1] ;
	    			qp[j] -= factor * qp[j-1] ;
	  		}
	  		for (j=NY-2 ; j >= 0 ; j--)
	    			Urans[j] = ( qp[j] - ae[j] * Urans[j+1] ) / ap[j] ;

	  		//  Auxiliary stuff
	  		tauw = muw * Urans[0] / ycrans[0] ;
	  		utau = sqrt( tauw / rhow ) ;
	  		lv = muw / ( rhow * utau ) ;
	  		for (j=0 ; j < NY ; j++) 
			{
	    			damping = 1.0 - exp( - ycrans[j] / ( lv * WALLMODEL_APLUS ) ) ;
				damping = damping*damping;
	   			MUTrans[j] = WALLMODEL_KAPPA * ycrans[j] * sqrt( Rrans[j] * tauw ) * damping ;
	    			ap[j] = 0.0 ;
	    			qp[j] = 0.0 ;   //  if LHS not zero, then should set qp to LHS * dyrans[j]...
	  		}
	
	  		//  Set up energy equation using finite volume
/*
			intfac = 1.0/(1.0+WALLMODEL_STRETCHING);
        		muint = (1.0 - intfac)*muw + intfac*MUrans[0];
			kint = (1.0 - intfac)*kw + intfac*Krans[0];
			uint = intfac*Urans[0];
        		mutint = intfac*MUTrans[0];
        		factor = (kint + C_P*mutint/WALLMODEL_PRT)/ycrans[0];
        		ap[0] = - factor;
			qp[0] = -(kint + C_P*mutint/WALLMODEL_PRT)*Tw/ycrans[0] + (muint + mutint) * uint * Urans[0]/ycrans[0];
*/
			ap[0] = - kw/ycrans[0];
			qp[0] = - kw/ycrans[0]*Tw;

	 		for (j=1 ; j < NY ; j++) 
			{
	     		    	intfac = ( yfrans[j] - ycrans[j-1] ) / ( ycrans[j] - ycrans[j-1] ) ;
	     			muint  = ( 1.0 - intfac ) * MUrans[j-1]  + intfac * MUrans[j] ;
	     			mutint = ( 1.0 - intfac ) * MUTrans[j-1] + intfac * MUTrans[j] ;
				kint = (1.0 - intfac ) * Krans[j-1] + intfac*Krans[j];
	     			factor = ( kint + C_P/WALLMODEL_PRT * mutint ) / ( ycrans[j] - ycrans[j-1] ) ;
		    		aw[j] = factor ;
	    			ap[j] -= factor ;
	    			ae[j-1] = factor ;
	    			ap[j-1] -= factor ;
	     			uint = ( 1.0 - intfac ) * Urans[j-1] + intfac * Urans[j] ;
	    			factor = ( muint + mutint ) * uint * ( Urans[j] - Urans[j-1] ) / ( ycrans[j] - ycrans[j-1] ) ;
	    			qp[j] += factor ;
	    			qp[j-1] -= factor ;
	  		}
	  		Trans[NY-1] = Tmatch ;

	  		//  Solve energy equation by single TDMA sweep
	  		for (j=1 ; j < NY-1 ; j++) 
			{
	    			factor = aw[j] / ap[j-1] ;
	    			ap[j] -= factor * ae[j-1] ;
	    			qp[j] -= factor * qp[j-1] ;
	  		}
	  		for (j=NY-2 ; j >= 0 ; j--)
	    			Trans[j] = ( qp[j] - ae[j] * Trans[j+1] ) / ap[j] ;

	 		if(Tmatch < Wall_Temp)
			    Tw = Trans[0];
	  		tauw = muw * Urans[0] / ycrans[0] ;
	  		rhow = rho0*Tmatch/ Tw ;
	  		utau = sqrt( tauw / rhow ) ;
	  		lv = muw / ( rhow * utau ) ;
	  		qw = kw* ( Trans[0] - Tw ) / ycrans[0] ;   // wall heat-flux

			converged = fabs((tauw - tauwold) / (tauwold + MACH_EPS));

			if( fabs((qw - qwold)/(qwold + MACH_EPS)) < converged )
				converged = fabs( (qw - qwold)/(qwold+MACH_EPS));
	 		iter ++; 
	          }while(converged > 1e-4);
      
                        flux[0] = - qw;
                        flux[1] = - tauw / Umatch * u ;
                        flux[2] = - tauw / Umatch * v ;
                        flux[3] = - muw * 4.0/3.0 * w/delta;

			free_these(10,Urans,Trans,Rrans,MUrans,MUTrans,Krans,aw,ap,ae,qp);
			free_these(2,dyrans,yfrans);
			debug_print("WALL_MODEL","Left oned_turbulence_boundarylayer_solver()\n");
}

EXPORT  void Set_up_merge_cell(
        Wave            *wave)
{
        RECT_GRID       *rgr = wave->rect_grid;
	int             *lbuf = rgr->lbuf;
	int             *ubuf = rgr->ubuf;
        int             dim = rgr->dim;
        int             *gmax = rgr->gmax;
        int             icoords[MAXD];
        int             imin[MAXD],imax[MAXD];
        int             i,j,k;
        Locstate        state;
        float           *coords;       
        float           dh[MAXD];
        float           vol;
	float		GL[3], GU[3];
	float dist_z;

        if (dim != 3) return;

        for (i = 0; i < dim; i++)
        {
            dh[i] = rgr->h[i];
            imin[i] = 0;
            imax[i] = gmax[i];
	    GL[3] = rgr->GL[3];
	    GU[3] = rgr->GU[3];
        }
	
	for (i = 0; i < dim; ++i)
	{
		if (lbuf[i] > 0) imin[i] -= lbuf[i];
		if (ubuf[i] > 0) imax[i] += ubuf[i];
	}

	for(i = imin[2]; i < imax[2]; i++)
	{
	    icoords[2] = i;
	    for( j = imin[1]; j < imax[1]; j++)
	    {
		icoords[1] = j;
		for( k = imin[0]; k < imax[0]; k++)
		{
		    icoords[0] = k;
		    coords = Rect_coords(icoords,wave);
                    state  = Rect_state(icoords,wave);		
		    if(coords[2] < 2.4 && coords[2] > 1.6)
		    {
			float xramp = (coords[2] - 1.6)/tan10 - 7.0;
			if(!is_obstacle_state(state) &&coords[0] < xramp && coords[0] + 0.104 > xramp)
			{
				Ramp_cellx(state) = YES;
			}
			else Ramp_cellx(state) = NO;				
		    }
		}
	    }
	}

        for (i = imin[0]; i < imax[0]; i++)
        {
            icoords[0] = i;

                for (j = imin[1]; j < imax[1]; j++)
                {
                    icoords[1] = j;
                    for (k = imin[2]; k < imax[2]; k++)
                    {
		        icoords[2] = k;
                        state  = Rect_state(icoords,wave);
                        coords = Rect_coords(icoords,wave);
			/*
			if(!is_obstacle_state(state) && coords[0] > GL[0] && coords[0] < GL[0] + dh[0])
				Inlet(state) = YES;
			else Inlet(state) = NO;

			 if(!is_obstacle_state(state) &&coords[0] < GU[0] && coords[0] > GU[0] - dh[0])
				 Outlet(state) = YES;
			 else  Outlet(state) = NO;
			*/

			  if(!is_obstacle_state(state) &&coords[2] < 3.9 && coords[2] > 3.9 - dh[2] && sqrt(coords[0]*coords[0] + coords[1]*coords[1]) > 0.1)
				   Reflect_wall(state) = YES;
			  else if(!is_obstacle_state(state) &&(coords[0] > ramp_x1 && coords[2] > 2.4  && coords[2] <  2.4 +  dh[2])||(coords[0] <ramp_x0 && coords[2] > 1.6 && coords[2] < 1.6 + dh[2]))  
			  {
				  Reflect_wall(state) = YES;
			   }
			/*
			if(!is_obstacle_state(state) && coords[2] > 3.9 && coords[2] < 3.9 + dh[2] && sqrt(coords[0]*coords[0] + coords[1]*coords[1]) < 0.1)
			{
				Nozzle(state) = YES;
				printf("coords = [%f  %f  %f]\n", coords[0], coords[1], coords[2]);
			}
			else
				Nozzle(state) = NO;*/

		    	dist_z = coords[2] - (ramp_z0 + tan10*(coords[0] - ramp_x0));
			if (dist_z > 0 && dist_z < 1.5*dh[2] && coords[0] > ramp_x0 && coords[0] < ramp_x1&& !is_obstacle_state(state))
				Merge_cell(state) = YES;
			else 	Merge_cell(state) = NO;	

			if ( coords[0] > ramp_x0 && coords[0] < ramp_x1 && dist_z > 0 && dist_z < dh[2] && !is_obstacle_state(state))
			{			
				Ramp_cellz(state) = YES;
			}	
			else 	Ramp_cellz(state) = NO;
		    }
		}
        }
}

EXPORT  void Update_merge_cell(
	float		dt,
        Wave            *wave,
	Front		*front)
{
        RECT_GRID       *rgr = wave->rect_grid;
	int		*lbuf = rgr->lbuf;
	int		*ubuf = rgr->ubuf;
        int             dim = rgr->dim;
        int             *gmax = rgr->gmax;
        int             icoords[3],icoords2[3];
        int             imin[3],imax[3];
        int             i, j, k, i0, i1, i2;
        Locstate        state,state2,state3;
	static Locstate 	newstate=NULL;
        float           *coords,*coords2;
        float           vol;
	float		du,dE,drho,pdens,p;
	float		ds0[6],ds1[6],dist[2];
	float		dtdh[3];
	float 		T,P,delta;
        if (dim != 3) return;
    	float GL[3],GU[3],dh[3], uvel[3],flux[4];
    	for(i = 0; i < dim; i++)
   	{
        	imin[i] = 0;
         	imax[i] = gmax[i];
         	dh[i] = rgr->h[i];
       	 	GL[i] = rgr->GL[i];
         	GU[i] = rgr->GU[i];
        }

	for (i = 0; i< 3; i++)  
		dtdh[i] = dt/dh[i];

        /*Boundary layer for the flat part. */
	float y0,yw;
	float   xm,ym;
	int     icoordsp[2][2][3];
        int     icoordsx[2],icoordsz[2];
        float   intpx, intpz;
	Locstate  statep[4];
	float   yd;       
	int      cover_by_nozzle;
	double   rate;
	/*Boundary layer treatment for side wall */	
        for (i0 = imin[0]; i0 < imax[0]; i0++)
        {
            icoords[0] = i0;
            for (i2 = imin[2]; i2 < imax[2]; i2++)
            {
                icoords[2] = i2;
                for (i1 = imin[1]; i1 < imax[1]; i1++)
                {
		    icoords[1] = i1;
		    state = Rect_state(icoords,wave);
		    coords = Rect_coords(icoords,wave);
		    //Front Wall
		    if(coords[1] < GU[1] && coords[1] >= GU[1] - dh[1] && !is_obstacle_state(state))
		    {
		    	icoords2[0] = icoords[0];
			icoords2[1] = icoords[1] - (tk - 2); //520m
			//icoords2[1] = icoords[1] - (tk - 1); //260m
			icoords2[2] = icoords[2];
			coords2 = Rect_coords(icoords2, wave);
			state2  = Rect_state(icoords2, wave);
			y0 = coords2[1];
			yw = GU[1];
               	    	for (i = 0; i < dim; ++i)
               	            uvel[i] = vel(i,state2);
		        for (i = 0; i < 4 ; i++)
			    flux[i] = 0.0;
		    	if ( fabs(uvel[0]) > MACH_EPS || fabs(uvel[1]) > MACH_EPS)
                    	    oned_turbulence_boundarylayer_solver(dim,uvel[2],uvel[0],uvel[1],y0,yw,flux,state2);
		
			Mom(state)[0] += flux[2]*dt/dh[1];
			Mom(state)[1] += flux[3]*dt/dh[1];
			Mom(state)[2] += flux[1]*dt/dh[1];
			Energy(state) += flux[0]*dt/dh[1];			
		     }
		    //Back Wall
		    else if( coords[1] > GL[1] && coords[1] <= GL[1] + dh[1] && !is_obstacle_state(state))
		    {
		    	icoords2[0] = icoords[0];
			icoords2[1] = icoords[1] + (tk - 2); //520m
			//icoords2[1] = icoords[1] + (tk - 1); //260m
			icoords2[2] = icoords[2];
			coords2 = Rect_coords(icoords2, wave);
			state2  = Rect_state(icoords2, wave);
			y0 = coords2[1];
			yw = GL[1];
               	    	for (i = 0; i < dim; ++i)
               	            uvel[i] = vel(i,state2);
		        for (i = 0; i < 4 ; i++)
			    flux[i] = 0.0;
		    	if ( fabs(uvel[0]) > MACH_EPS || fabs(uvel[1]) > MACH_EPS)
                    	    oned_turbulence_boundarylayer_solver(dim,uvel[2],uvel[0],uvel[1],y0,yw,flux,state2);

			Mom(state)[0] += flux[2]*dt/dh[1];
			Mom(state)[1] += flux[3]*dt/dh[1];
			Mom(state)[2] += flux[1]*dt/dh[1];
			Energy(state) += flux[0]*dt/dh[1];
		     }
		}
	    }
	}
	/*Boundary layer treatment for up wall, low wall and ramp*/
        for (i0 = imin[0]; i0 < imax[0]; i0++)
        {
            icoords[0] = i0;
            for (i1 = imin[1]; i1 < imax[1]; i1++)
            {
                icoords[1] = i1;
                for (i2 = imin[2]; i2 < imax[2]; i2++)
                {
                    icoords[2] = i2;
                    state = Rect_state(icoords,wave);
                    coords = Rect_coords(icoords,wave);
		    if( coords[2] < GU[2] - 0.1 && coords[2] >= GU[2] - 0.1 - dh[2] && sqrt(coords[0]*coords[0] + coords[1]*coords[1]) > 0.1)
		    {
			icoords2[0] = icoords[0];
			icoords2[1] = icoords[1];
			icoords2[2] = icoords[2] - (tk - 2); //520m
			//icoords2[2] = icoords[2] - (tk - 1);   //260m
			coords2 = Rect_coords(icoords2, wave);
			state2  = Rect_state(icoords2, wave);
			y0 = coords2[2];
			yw = GU[2] - 0.1;
               	    	for (i = 0; i < dim; ++i)
               	            uvel[i] = vel(i,state2);
		        for (i = 0; i < 4 ; i++)
			    flux[i] = 0.0;
		    	if ( fabs(uvel[0]) > MACH_EPS || fabs(uvel[1]) > MACH_EPS)
                    	    oned_turbulence_boundarylayer_solver(dim,uvel[0],uvel[1],uvel[2],y0,yw,flux,state2);
			compute_square_circle(0,0,0.099,coords[0],coords[1],0.5*dh[0],0.5*dh[0],&cover_by_nozzle, &rate);
			for( i = 0; i < 3; i++)
				Mom(state)[i] += rate*flux[i + 1]*dt/dh[2];
			Energy(state) += rate*flux[0]*dt/dh[2];

			    Boundary_flux(state)[0] = flux[0];
			    Boundary_flux(state)[1] = flux[1];
			    Boundary_flux(state)[2] = flux[2];
			    Boundary_flux(state)[3] = flux[3];
		    }
		    else if((coords[0] > ramp_x1 && coords[2]  > ramp_z1 && coords[2] <= ramp_z1 + dh[2])||(coords[0] <ramp_x0 && coords[2] > ramp_z0 && coords[2] <= ramp_z0 + dh[2]))
		    {
			icoords2[0] = icoords[0];
			icoords2[1] = icoords[1];
			icoords2[2] = icoords[2] + (tk - 2); //520
			//icoords2[2] = icoords[2] + (tk - 1); //260
			coords2 = Rect_coords(icoords2, wave);
			state2  = Rect_state(icoords2, wave);
			y0 = coords2[2];
			if( coords[0] > ramp_x1 && coords[2] - ramp_z1 > 0) yw = ramp_z1;
			else yw = ramp_z0;
               	    	for (i = 0; i < dim; ++i)
               	            uvel[i] = vel(i,state2);
		        for (i = 0; i < 4 ; i++)
			    flux[i] = 0.0;
		    	if ( fabs(uvel[0]) > MACH_EPS || fabs(uvel[1]) > MACH_EPS)
                    	    oned_turbulence_boundarylayer_solver(dim,uvel[0],uvel[1],uvel[2],y0,yw,flux,state2);
			for( i = 0; i < 3; i++)
				Mom(state)[i] += flux[i + 1]*dt/dh[2];
			Energy(state) += flux[0]*dt/dh[2];
			    Boundary_flux(state)[0] =  flux[0];
			    Boundary_flux(state)[1] =  flux[1];
			    Boundary_flux(state)[2] =  flux[2];
			    Boundary_flux(state)[3] =  flux[3];
		    }
		    else if(coords[0] > ramp_x0 && coords[0] < ramp_x1 && coords[2] - (ramp_z0 + tan10*(coords[0] - ramp_x0)) >0 && coords[2] - (ramp_z0 + tan10*(coords[0] - ramp_x0)) <= dh[2])
		    {
			delta = (tk  - 0.5)*dh[2];
			yd = ramp_z0 + tan10*(coords[0] - ramp_x0);
			xm = coords[0] - delta*sin10;
			ym = yd + delta*cos10;
			
			icoordsx[0] = icoords[0] - 1;
			icoordsx[1] = icoords[0];
			if(ym >= coords[2] + (tk-1)*dh[2] && ym < coords[2] + tk*dh[2]) 
			    icoordsz[0] = icoords[2] + tk -1;
			else if(ym > coords[2] + (tk -2)*dh[2] && ym < coords[2] + (tk-1)*dh[2])
			    icoordsz[0] = icoords[2] + tk - 2;
			else
			    printf("error in the icoordsy ym = %f  coords[2] = %f  yd = %f delta = %f cos10 = %f tk = %d \n",ym,coords[2],yd, delta,cos10,tk );
			icoordsz[1] = icoordsz[0] + 1;
			
			for( i = 0; i < 2; i++)
			    for( j = 0; j < 2; j++)
				icoordsp[i][j][1] = icoords[1];
			for( i = 0; i < 2; i++)
			{
			    icoordsp[i][0][2] = icoordsz[0];
			    icoordsp[i][1][2] = icoordsz[1];
			    icoordsp[0][i][0] = icoordsx[0];
			    icoordsp[1][i][0] = icoordsx[1];
			}
			
			for(i = 0; i < 2; i++)
			    for(j = 0; j < 2; j++)
		        	statep[2*i + j] = Rect_state(icoordsp[i][j], wave);
			
			intpx = (coords[0] - xm)/dh[0];
			intpz = (coords[2] + (icoordsz[1] - icoords[2])*dh[2] - ym)/dh[2];
			float *coords3;
			coords3 = Rect_coords(icoordsp[0][0],wave);			
			float pz = 1 - (ym - coords3[2])/dh[2];
			if(fabs(pz - intpz) > 10* MACH_EPS) printf("error in intpz computation!\n");
			if( intpx < 0 || intpx > 1 || intpz < 0 || intpz > 1) 
				printf("intpx is %f, intpy is %f, xm is %f  ym is %f coords[0] = %f, coords[2] = %f\n", intpx, intpz,xm,ym,coords[0],coords[2]);

			if (newstate == NULL)
				g_alloc_state(&newstate,Params(state)->sizest); 	

			//assign the momentum
        		for (i = 0; i < dim; ++i)
            		     Mom(newstate)[i] =  intpz*(intpx*Mom(statep[0])[i] + (1 - intpx)*Mom(statep[2])[i] ) + (1 - intpz)*(intpx*Mom(statep[1])[i]  + (1 - intpx)*Mom(statep[3])[i]);
			//assign the energy
			Energy(newstate) =  intpz*(intpx*Energy(statep[0]) + (1 - intpx)*Energy(statep[2])) + (1 - intpz)*(intpx*Energy(statep[1]) + (1 - intpx)*Energy(statep[3]));
			//assign the density
			Dens(newstate) =  intpz*(intpx*Dens(statep[0]) + (1 - intpx)*Dens(statep[2])) + (1 - intpz)*(intpx*Dens(statep[1]) + (1 - intpx)*Dens(statep[3]));
			//assign the partial density
			for( i = 0; i < Params(state)->n_comps ; i++)
			{
			    pdens(newstate)[i] =  intpz*(intpx*pdens(statep[0])[i] + (1 - intpx)*pdens(statep[2])[i]) + (1 - intpz)*(intpx*pdens(statep[1])[i] + (1 - intpx)*pdens(statep[3])[i]);
			}
			//assign the Eos parameters
        		Set_params(newstate,statep[0]);
        		set_type_of_state(newstate,GAS_STATE);

			
			//assin the velocity 
			for( i = 0; i < 3; i++)
			uvel[i] = intpz*(intpx*vel(i,statep[0]) + (1 - intpx)*vel(i,statep[2]))
				+(1 - intpz)*(intpx*vel(i,statep[1]) + (1 - intpx)*vel(i,statep[3]));
			float u[3];/*velocity which are parallel and normal to the ramp.*/
		    	u[0] = uvel[0]*cos10 + uvel[2]*sin10;
		    	u[1] = uvel[1];
		    	u[2] = uvel[2]*cos10 - uvel[0]*sin10;
		    	for (i = 0; i < 4 ; i++)
			flux[i] = 0.0;
 		        yw = yd;
			y0 = yd + delta;
		    	if ( fabs(u[0]) > MACH_EPS || fabs(u[1]) > MACH_EPS)
                    	    oned_turbulence_boundarylayer_solver(dim,u[0],u[1],u[2],y0,yw,flux,newstate);
			float f[3];
		        f[0] = (flux[1]*cos10 - flux[3]*sin10)/cos10;
		        f[1] = flux[2]/cos10;
		        f[2] = (flux[1]*sin10 + flux[3]*cos10)/cos10;

			Flux(state)[4][2] =  - Flux(state)[4][4]*tan10;
			for( i = 0; i < 3; i++)
				Flux(state)[4][2+i] += f[i];
			Flux(state)[4][1] += flux[0];
			    Boundary_flux(state)[0] =  flux[0];
			    Boundary_flux(state)[1] =  flux[1];
			    Boundary_flux(state)[2] =  flux[2];
			    Boundary_flux(state)[3] =  flux[3];
		    }
		}
	    }
	}      

        return;
	for (i0 = imin[0]; i0 < imax[0]; i0++)
        {
            icoords[0] = i0;
            for (i1 = imin[1]; i1 < imax[1]; i1++)
            {
                icoords[1] = i1;
                for (i2 = imin[2]; i2 < imax[2]; i2++)
                {
                    icoords[2] = i2;
                    state = Rect_state(icoords,wave);
                    coords = Rect_coords(icoords,wave);		   
		    dist[0]= coords[2] - (ramp_z0 + tan10*(coords[0] - ramp_x0));
		    if ( !is_obstacle_state(state) && Merge_cell(state) && flux_flag(state) == 1)
                    {
			dist[0] = coords[2] - (ramp_z0 + tan10*(coords[0] - ramp_x0));
			if(dist[0] < 0.5*dh[2])  /*The fisrt cell is a > 0.5 cell and then second cell is a whole cell.*/
			{
			    /*The area of six side of the first cell.*/
			    ds0[0] = (dist[0] + 0.5*dh[0]*tan10)/dh[2] + 0.5;
			    ds0[1] = (dist[0] - 0.5*dh[0]*tan10)/dh[2] + 0.5;
			    ds0[2] = dist[0]/dh[2] + 0.5;
			    ds0[3] = dist[0]/dh[2] + 0.5;
			    ds0[4] = 1.0;
			    ds0[5] = 1.0;
			    vol = dist[0]/dh[2] + 0.5;

			    /*The second cell.*/
			    icoords2[0] = i0;
			    icoords2[1] = i1;
			    icoords2[2] = i2 + 1;
			    state2 = Rect_state(icoords2,wave);
			    coords2 = Rect_coords(icoords2,wave);
			    dist[1] = dist[0] + dh[2];
			    for( i= 0; i< 6; i++) 
				ds1[i] = 1.0;

			    /*The third cell*/
			    int icoords3[3];
			    for( i = 0; i < 3; i++) 
				icoords3[i] = icoords[i];
			    icoords3[0]+= 1;
			    state3 = Rect_state(icoords3,wave);
			    if(is_obstacle_state(state3))
			    {
				for(j = 0; j < 5; j++)
					Flux(state)[1][j] = Flux(state2)[1][j];
				if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
				{
					Gas_param *params = Params(state);
					int    num_comps;
					if((params != NULL) && ((num_comps = params->n_comps) != 1))
					{
						for(j = 0; j < num_comps; j ++)
							Flux(state)[1][j + 5] = Flux(state2)[1][j+5];
					}
				}
			    }
                            //Update the merged cells by using Flux(state)[][] and then split.
			    /*average the increase in density,energy,and momentum for two cells.*/
			    /*update density*/
			    drho = 0;
			    for( i = 0; i < 3 ; i++)
			    {
				drho += (Flux(state )[2*i][0]*ds0[2*i] - Flux(state )[2*i+1][0]*ds0[2*i+1])*dt/dh[i];
				drho += (Flux(state2)[2*i][0]*ds1[2*i] - Flux(state2)[2*i+1][0]*ds1[2*i+1])*dt/dh[i];
			    }
			    drho = drho/(vol + 1);
			    Dens(state)   = 0.5*(Dens(state)+Dens(state2))+drho;
			    Dens(state2)  = Dens(state);						   
			    /*update partial density.*/
            		    if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            		    {
                		Gas_param *params = Params(state);
                		int    num_comps;
                		if((params != NULL) && ((num_comps = params->n_comps) != 1))
                		{
                		    for(j = 0; j < num_comps; j ++)
				    {
					drho = 0;
					for( i = 0; i < 3; i ++)
					{
					 drho +=(Flux(state )[2*i][5+j]*ds0[2*i]-Flux(state )[2*i+1][5+j]*ds0[2*i+1])*dt/dh[i];
					 drho +=(Flux(state2)[2*i][5+j]*ds1[2*i]-Flux(state2)[2*i+1][5+j]*ds1[2*i+1])*dt/dh[i];
					}
					drho = drho/(1 + vol);
					pdens(state)[j]  = 0.5*(pdens(state)[j]+pdens(state2)[j])+drho;
					pdens(state2)[j] = pdens(state)[j];
				    } 
                		}
            		    }
			    /*update the momentum.*/			    
			    for(j = 0;j < 3;j++ )
			    {
				du = 0;
				for( i = 0; i < 3 ; i++)
			        {
				    du += (Flux(state )[2*i][j+2]*ds0[2*i] - Flux(state )[2*i+1][j+2]*ds0[2*i+1])*dt/dh[i];
				    du += (Flux(state2)[2*i][j+2]*ds1[2*i] - Flux(state2)[2*i+1][j+2]*ds1[2*i+1])*dt/dh[i];
			        }
				du = du/(1 + vol);
				Mom(state )[j] = 0.5*(Mom(state)[j]+Mom(state2)[j]) + du;
				Mom(state2)[j] = Mom(state)[j];
			    }
			    /*update the energy.*/
		      	    dE = 0;
			    for( i = 0; i < 3 ; i++ )
			    {
				dE += (Flux(state )[2*i][1]*ds0[2*i]- Flux(state )[2*i+1][1]*ds0[2*i+1])*dt/dh[i];
				dE += (Flux(state2)[2*i][1]*ds1[2*i]- Flux(state2)[2*i+1][1]*ds1[2*i+1])*dt/dh[i];
			    }
			    dE = dE/(vol + 1);
			    Energy(state)   = 0.5*(Energy(state)+Energy(state2)) + dE;
			    Energy(state2)  = Energy(state); 
			}
		        else if(dist[0] > 0.5*dh[2] && dist[0] <= dh[2])/*Merge the cell with a < half cell behind it.*/
			{
			    vol = dist[0]/dh[2] + 0.5;
			    ds0[0] = (dist[0] + 0.5*dh[0]*tan10)/dh[2] + 0.5;
			    ds0[1] = (dist[0] - 0.5*dh[0]*tan10)/dh[2] + 0.5;
			    ds0[2] = vol;
			    ds0[3] = vol;
			    ds0[4] = 1.0;
			    ds0[5] = 1.0;
			
                            //Update the merged cells by using Flux(state)[][] and then split.
			    	/*update the density in each cell.*/
			    	drho = 0;
			    	for( i = 0; i < 3; i++)
			        	drho += (Flux(state)[2*i][0]*ds0[2*i] - Flux(state)[2*i+1][0]*ds0[2*i+1])*dt/dh[i]; 
			    	drho =  drho/vol;
			    	Dens(state) += drho ;
			    	/*update the partial density.*/
            		    	if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            		    	{
                			Gas_param *params = Params(state);
                			int    num_comps;
                			if((params != NULL) &&
                			   ((num_comps = params->n_comps) != 1))
                			{				    
                		    	    for(j = 0; j < num_comps; j ++)
				    	    {
						drho = 0;				
						for( i = 0; i < 3; i ++)
						{
					    	    drho += (Flux(state)[2*i][5 + j]*ds0[2*i]
						    	   - Flux(state )[2*i+1][5 + j]*ds0[2*i+1])*dt/dh[i];
						}
					 	drho = drho/vol;
					 	pdens(state)[j]  += drho;
				    	    } 
                			}
            		    	}
			    	/*update the energy in each cell.*/
			    	dE = 0;
			    	for( i = 0; i < 3; i++)
					dE += (Flux(state)[2*i][1]*ds0[2*i]- Flux(state)[2*i+1][1]*ds0[2*i+1])*dt/dh[i]; 
			    	dE = dE/vol;
			    	Energy(state) += dE;

			    	/*update the momentum in each cell.*/
			    	for(j =0;j<3;j++)
			    	{
				    du = 0;
				    for( i = 0; i < 3; i++)
				    	du += (Flux(state)[2*i][2+j]*ds0[2*i]- Flux(state)[2*i+1][2+j]*ds0[2*i+1])*dt/dh[i]; 
				    du = du/vol;
				    Mom(state)[j] += du;
			    	}
			    }

			 	

                    }/*end of the update of merged cell.*/
                }
	    }
	}/*end of 3 space loop*/

}

void	get_constant_state(Locstate, int, float*, float);

void turbulent_inflow(
	Wave *wave,
	Front *front)
{
        int 	d, i, j, k, l, m, n, dim = front->rect_grid->dim;
	float 	dh = front->rect_grid->h[0]; //mesh width fixed for x,y,z
        int 	Ix, Iy, Iz; //goes to input 
	int 	nx, ny, nz;     // Number of grid points 
	int 	NFx , NFy, NFz; //filter size 
	int 	Mx, My, Mz; 
	
	// Integral length scales
        Ix = 1;
        Iy = 1;
        Iz = 1;

	nx = (int)Ix/dh;
       	ny = (int)Iy/dh;
	nz = (int)Iz/dh;

	NFx = 2*nx;
	NFy = 2*ny;
	NFz = 2*nz;

	// My X Mz : dimensions of the computational grid in the inflow plane
	Mx = front->rect_grid->gmax[0];
	My = front->rect_grid->gmax[1];
	Mz = front->rect_grid->gmax[2];

	// Convolution coefficients
        float b_base_x[2*NFx+1], b_x[2*NFx+1], b_sum_x = 0;
        float b_base_y[2*NFy+1], b_y[2*NFy+1], b_sum_y = 0;
        float b_base_z[2*NFz+1], b_z[2*NFz+1], b_sum_z = 0;

        for(i=0; i<2*NFx+1; ++i)
        {
            b_base_x[i] = exp(-PI*i/nx);
            b_sum_x = b_sum_x + b_base_x[i]*b_base_x[i];
        }

        b_sum_x = sqrt(b_sum_x);
        for(i=0; i<2*NFx+1; ++i)
            b_x[i] = b_base_x[i]/b_sum_x;

        for(i=0; i<2*NFy+1; ++i)
        {
            b_base_y[i] = exp(-PI*i/ny);
            b_sum_y = b_sum_y + b_base_y[i]*b_base_y[i];
        }

        b_sum_y = sqrt(b_sum_y);
        for(i=0; i<2*NFy+1; ++i)
            b_y[i] = b_base_y[i]/b_sum_y;

	for(i=0; i<2*NFz+1; ++i)
        {
            b_base_z[i] = exp(-PI*i/nz);
            b_sum_z = b_sum_z + b_base_z[i]*b_base_z[i];
        }

        b_sum_z = sqrt(b_sum_z);
        for(i=0; i<2*NFz+1; ++i)
            b_z[i] = b_base_z[i]/b_sum_z;

	//Initialize and store random fields of dimension
	/*
	float**** random_numbers;
	quad_array(&random_numbers, 3, 2*NFx+1, My+2*NFy, Mz+2*NFz, sizeof(float) );
        
	for(d=0;d<dim;++d)
        {
            for(i=0; i<2*NFx+1; ++i)
	    {
            	for(j=0; j<My+2*NFy; ++j)
              	  for(k=0; k<Mz+2*NFz; ++k)
                	random_numbers[d][i][j][k] = random_gaussian(0,1,seed);
	    }
        }*/
        unsigned short int seed[3]={pp_mynode()+front->step,pp_mynode()+front->step,pp_mynode()+front->step};

        float	 S11, S22, S33, intensity=0.023; //goes to input 
	Locstate st;
	int   	 st_type, icoords[3],icoords1[3];
	float 	 *coords;
	float    GLx = front->rect_grid->GL[0];

        float 	u_field[3][My][Mz], u_perturb[3][My][Mz];	
        float 	dt = front->dt; // time step	
	float   u_x = 189.4;
        float 	tau_t = Ix/u_x; // Lagrangian time step
        float 	filter_exp1 = sqrt(1-exp(-PI*dt/tau_t));
        float 	filter_exp2 = exp(-PI*dt/(2*tau_t));
/* TK for nozzle + outflow */
	for(i=0; i<Mx+front->rect_grid->ubuf[0]; ++i)
        {	
	   icoords[0] = i;
           for(j=0; j<My; ++j)
	   {	
	     icoords[1] = j;
	     for(k=0; k<Mz; ++k)
	     {	
	      icoords[2] = k;

	      coords = Rect_coords(icoords,wave);
	      st = Rect_state(icoords, wave);
	      if (coords[2]> (front->rect_grid->GU[2]-0.1) && !is_obstacle_state(st))
		get_constant_state(st, 2, coords, front->time);
	      if(coords[0]> (front->rect_grid->GU[0]-dh) && !is_obstacle_state(st))
	      {
		st_type = state_type(st);
		set_state(st, TGAS_STATE, st);
		Press(st) = wall_pressure_avg;
		set_state(st, st_type, st);
	      }	

	    }	
	   }
	 }			
/* TK for nozzle */	
        if(dt == 0) //means  initial step
        {	
	   icoords[0] = 0;
           for(j=0; j<My; ++j)
	   {	
	     icoords[1] = j;
	     for(k=0; k<Mz; ++k)
	     {	
	      icoords[2] = k;

	      st = Rect_state(icoords, wave);
		
	      for(d=0; d<3; ++d)
	         Vel_Field(st)[d] = 0.0; 	
	     }
	   }	
	   return ;
        }	

        for(d=0; d<dim; ++d)
          for(j=0; j<My; ++j)
            for(k=0; k<Mz; ++k)
		    u_field[d][j][k] = 0;
/*
	icoords[0] = 0;
        for(d=0; d<dim; ++d){
        	for(j=0; j<My; ++j){
	    			icoords[1] = j;
            			for(k=0; k<Mz; ++k){
	    				icoords[2] = k;
	    				//coords = Rect_coords(icoords, wave);  DO WE NEED THIS??????
	    				st = Rect_state(icoords, wave);
            				for(l=0; l<2*NFx+1; ++l){
              					for(m=0; m<2*NFy+1; ++m){
                 					for(n=0; n<2*NFz+1; ++n){
           	    						u_field[d][j][k] += b_x[l]*b_y[m]*b_z[n]*random_numbers[d][l][j+m][k+n];
							}
						}
					}
	      	    			u_perturb[d][j][k] = Vel_Field(st)[d]*filter_exp2 + u_field[d][j][k]* filter_exp1;
	    			}	
	  	}
	}
*/

	float** random_numbers;
	bi_array(&random_numbers, My+2*NFy, Mz+2*NFz, sizeof(float) );
	
	icoords[0] = 0;
        for(d=0; d<dim; ++d){
            	for(l=0; l<2*NFx+1; ++l){
			load_random_nums( random_numbers, My+2*NFy, Mz+2*NFz, 0, 1, seed);
        		for(j=0; j<My; ++j){
	    			icoords[1] = j;
            			for(k=0; k<Mz; ++k){
	    				icoords[2] = k;
	    				st = Rect_state(icoords, wave);
              					for(m=0; m<2*NFy+1; ++m){
                 					for(n=0; n<2*NFz+1; ++n){
           	    						u_field[d][j][k] += b_x[l]*b_y[m]*b_z[n]*random_numbers[j+m][k+n];
							}
						}
					}
	      	    			u_perturb[d][j][k] = Vel_Field(st)[d]*filter_exp2 + u_field[d][j][k]* filter_exp1;
	    			}	
	  	}
	}

	//Prescribed Reynolds stress tensor (isotropic)
	S11= S22 = S33 = (intensity)*(intensity)*u_x*u_x/3; 

	icoords[0] = 0;
	for(j=0; j<My; ++j)
	{
	    icoords[1] = j;
	    for(k=0; k<Mz; ++k)
	    {
	    	icoords[2] = k;

	  	coords = Rect_coords(icoords, wave);
	    	st = Rect_state(icoords, wave);
	    if(coords[0] < GLx+(icoords[0]+1)*dh && coords[2]>(front->rect_grid->GL[2]+0.1+dh) && coords[2]<(front->rect_grid->GU[2]-0.1-dh))
	    {
	      if(!is_obstacle_state(st))
	      {
		st_type = state_type(st);
		set_state(st, TGAS_STATE, st);

	    	Vel(st)[0] = u_x + sqrt(S11)*u_perturb[0][j][k];
	    	Vel(st)[1] = sqrt(S22)*u_perturb[1][j][k];
	    	Vel(st)[2] = sqrt(S33)*u_perturb[2][j][k];

	     	Vel_Field(st)[0] = u_field[0][j][k]; 	
		Vel_Field(st)[1] = u_field[1][j][k]; 	
		Vel_Field(st)[2] = u_field[2][j][k]; 	

	        Press(st) = 0.4;
		Dens(st) = 12.41e-5;
                for(i = 0; i < Params(st)->n_comps; i++)
                    pdens(st)[i] = 0.0;
                pdens(st)[0] = Dens(st);

		set_state(st, st_type, st);
	      }
	      else
		printf("TMP_TK some obstacle state %d,%d,%d\n",icoords[0],j,k);	
	    }
	  }
	}
	free(random_numbers);

} 

LOCAL	void	load_random_nums ( 
			float** random_numbers,    	/* Stores random numbers */   
			int m, 				/* Number of rows in random_numbers*/
			int n, 				/* Number of columns in random_numbers*/
			float avg, 			/* Average of Gaussian distribution*/
			float stdv,			/* Standard Deviation of Gaussian Distribution */
			unsigned short int seed[3])	/* Seed to be used for producing random variable */
{
	int i, j;
	for ( i = 0 ; i < m; i++ ){
		for ( j = 0; j < n; j++ ){
			random_numbers[i][j] = random_gaussian(avg, stdv, seed);
		}
	}
}

LOCAL	void compute_square_circle(
		double xc,//center of the nozzle
		double yc,
		double radius,//radius of the nozzle
		double xs,//center of the cell
		double ys,
		double xl,//lenght of the cell
		double yl,
		int *cover_by_nozzle,// whether the cell is in the nozzle or not
		double *rate)//ratio of the cell in the nozzle
{
	double xf,yf,xn,yn,l0,l1;
	bool sdir[2][2];
	int i,j;
	for(i = 0; i < 2;i++)
		for(j = 0; j < 2; j++)
			sdir[i][j] = 1;
	
	if(xs >xc)
	{
	    xf = xs + 0.5*xl;
	    xn = xs - 0.5*xl;
	    for( i = 0; i < 2; i++)
		    sdir[0][i]= 0;
	}
	else
	{
	    xf = xs - 0.5*xl;
	    xn = xs + 0.5*xl;
	    for( i = 0; i < 2; i++)
		    sdir[1][i]=0;
	}
	if(ys > yc)
	{
	    yf = ys + 0.5*yl;
	    yn = ys - 0.5*yl; 
	    for(i = 0; i <2; i++)
		    sdir[i][0] =0;
	}
	else
	{
	    yf = ys - 0.5*yl;
	    yn = ys + 0.5*yl;
	    for( i = 0; i < 2; i++)
		    sdir[i][1] = 0;
	}
	
	if( (xf - xc)*(xf - xc) + (yf - yc)*(yf - yc) < radius*radius)//The square lies in the nozzle
	{
	    *cover_by_nozzle = 1;
	    *rate = 0.0;
	}
	else if( (xn - xc)*(xn - xc) + (yn - yc)*(yn - yc) > radius*radius)// The square lies out the nozzle
	{
	    *cover_by_nozzle = -1;
	    *rate = 1.0;
	}
	else //part of the square is in the nozzle
	{
	    *cover_by_nozzle = 0;
	    float sq[2][2][2];
	    int sq_cover = 0;
	    bool cover[2][2];
	    for(i = 0; i < 2; i++)
 	    {
		sq[0][i][0] = xs - 0.5*xl;
		sq[1][i][0] = xs + 0.5*xl;
		sq[i][0][1] = ys - 0.5*yl;
		sq[i][1][1] = ys + 0.5*yl;
	    }
	    for(i = 0; i < 2; i++)
	    {
		for( j = 0; j < 2; j++)
		{
		    if( (sq[i][j][0] - xc)*(sq[i][j][0] - xc) + (sq[i][j][1] - yc )*(sq[i][j][1] - yc) < radius*radius)
			    cover[i][j] = 1;//whether the point is in the nozzle
		    else
			    cover[i][j] = 0;
		}
	    }
	
	
	    int ri,rj;
	    int test=0;
	for(i = 0; i < 2; i++)
	{
	   for( j = 0; j < 2; j ++)
		   test += sdir[i][j];
	}
	if(test != 1) printf("error in sdir!\n");
	for( i = 0; i < 2; i++)
	{
	    for( j = 0; j < 2; j++)
	    {
		if(sdir[i][j] == 1)
		{
		   ri = 1 - i;
		   rj = 1 - j;
		   if(cover[i][j] == 1 || cover[ri][rj] == 0)
		   {
			printf("error in the rate compuation! cover = [%d  %d  %d  %d]\n", cover[0][0],cover[0][1],cover[1][0],cover[1][1]);
			printf("sq[i][j] = [%f  %f], sq[ri][rj] = [%f  %f] sf = [%f  %f] sn = [%f %f]\n", sq[i][j][0], sq[i][j][1], sq[ri][rj][0], sq[ri][rj][1], xf, yf, xn, yn);
		   }
		   if(cover[i][rj] == 1 && cover[ri][j] == 1)
		   {
			l0 = fabs(sq[i][j][0] - xc) - sqrt( radius*radius - (sq[i][j][1] - yc)*(sq[i][j][1] - yc));
			l1 = fabs( sq[i][j][1] - yc) - sqrt( radius*radius - (sq[i][j][0] - xc)*(sq[i][j][0] - xc));
			*rate = 1 - 0.5*l0*l1/(xl*yl);
		   }
		   else if( cover[i][rj] == 1 && cover[ri][j] == 0)
		   {
			l0 = sqrt(radius*radius - (sq[i][j][0] - xc)*(sq[i][j][0] - xc)) - fabs(sq[ri][rj][1] - yc);	
			l1 = sqrt(radius*radius - (sq[ri][rj][0] - xc)*(sq[ri][rj][0] - xc)) - fabs(sq[ri][rj][1] - yc);
			*rate = 0.5*(l0 + l1)*xl/(xl*yl);
		   }
		   else if(cover[i][rj] == 0 && cover[ri][j] == 1)
		   {
			l0 = sqrt( radius * radius - (sq[ri][j][1] - yc)*(sq[ri][j][1] - yc)) - fabs(sq[ri][j][0] - xc);
			l1 = sqrt( radius * radius - (sq[ri][rj][1] - yc)*(sq[ri][rj][1] - yc)) - fabs(sq[ri][rj][0] - xc);
			*rate = 0.5*(l0+l1)*yl/(xl*yl);
		   }
		   else
		   {
			l0 = sqrt(radius * radius - (sq[ri][rj][1] - yc)*(sq[ri][rj][1] - yc)) - fabs(sq[ri][rj][0] - xc);
			l1 = sqrt(radius * radius - (sq[ri][rj][0] - xc)*(sq[ri][rj][0] - xc)) - fabs(sq[ri][rj][1] - yc);
			*rate = 0.5*l0*l1/(xl*yl);
		   }
		}
	    }
	}
	}
	if(*rate <0 || *rate > 1.0) printf("error! area less than 0 or area largier than 1\n");
	return;

}

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)<(b))?(b):(a))
double average_p_if_mach(Wave *wave, Grid *grid)
{
    RECT_GRID *gr = grid->rect_grid;
    double xr[2]={(gr->GU[0]-gr->h[0]), gr->GU[0]};
    double yr[2]={gr->GL[1], gr->GU[1]};
    double zr[2]={gr->GL[2], gr->GU[2]};

    int ir[2];
    int jr[2];
    int kr[2];

    ir[0] = max(0, (int)((xr[0] - gr->L[0])/gr->h[0]));
    ir[1] = min(gr->gmax[0],ceil((xr[1]-gr->L[0])/gr->h[0]));

    jr[0] = max(0, (int)((yr[0] - gr->L[1])/gr->h[1]));
    jr[1] = min(gr->gmax[1],ceil((yr[1]-gr->L[1])/gr->h[1]));

    kr[0] = max(0, (int)((zr[0] - gr->L[2])/gr->h[2]));
    kr[1] = min(gr->gmax[2], ((zr[1]-gr->L[2])/gr->h[2]));


    //printf("ir = %d : %d\n",ir[0],ir[1]);
    //printf("jr = %d : %d\n",jr[0],jr[1]);
    //printf("kr = %d : %d\n",kr[0],kr[1]);
    //printf("gmax = %d \n",gr->gmax[0]);
    //printf("L = %f %f %f\n",gr->L[0],gr->L[1],gr->L[2]);
    //printf("U = %f %f %f\n",gr->U[0],gr->U[1],gr->U[2]);
    //exit(0);
    int num_cells=0;
    double p_sum = 0.0f;
    for(int i=ir[0]; i < ir[1]; ++i)
    for(int j=jr[0]; j < jr[1]; ++j)
    for(int k=kr[0]; k < kr[1]; ++k)
    {
	int icoords[3] = {i,j,k};
	int comp = Rect_comp(icoords,wave);
//	printf("-------\n");
//	printf("icoords= %d %d %d\n",icoords[0],icoords[1],icoords[2]);
//	printf("comp=%d\n",comp);
	if(comp ==3)
	{

	    Locstate st = Rect_state(icoords,wave);
	    double *vel = Vel(st);
	    double abs_vel[3];
	    for(int i=0; i < 3; ++i)
	    {
		abs_vel[i]=fabs(vel[i]);
	    }
	    double m = mach_number(st, abs_vel);
	    if(m>=1.0)
	    {
		double p = pressure(st);
		p_sum +=p;
		++num_cells;
	    }
	    else
	    {
//		printf("Mach too small: %f\n",m);
	    }
	}

    }
    pp_global_sum(&p_sum, 1);
    pp_global_isum(&num_cells, 1);
    double p_avg = p_sum / num_cells;
    if(p_avg > 2.0)
	p_avg = 2.0;
    wall_pressure_avg = p_avg; // Set the global variable
	printf("wall_pressure_avg = %f\n",p_avg);
    return p_avg;
}

#endif /* defined(DEBUG_MUSCL) */
