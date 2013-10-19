/************************************************************************************
FronTier is a set of libraries that implements differnt types of Front Traking algorithms.
Front Tracking is a numerical method for the solution of partial differential equations 


ne
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
*				gmuscl_chem.c
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*       New and improved version. Contains the drivers for the MUSCL scheme.
*/


#include <ghyp/ghyp.h>
#include <gdecs/vecdecs.h>
#include "gchem.c"
#include <gdecs/geos.h>
#include <gdecs/gstate.h>
#include <geos/mpoly.h>
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
LOCAL   void    oned_turbulence_boundarylayer_solver(int,double,double,double,double,double,double*,Locstate,int);

#if defined(DEBUG_MUSCL)
LOCAL	const char **set_muscl_debugging(int,int*,int*,Tan_stencil*,Front*,int);
LOCAL	const char **toggle_muscl_debugging(const char*,int,int*,int*,
					    Front*,int);
LOCAL	void	reset_muscl_debugging(const char**);
#endif /* defined(DEBUG_MUSCL) */


/*
*			oned_MUSCL():
*
*	One dimensional MUSCL code.   The input data are the conservative
*	state variables vst and the source terms src.   The MUSCL is a five
*	point scheme, therefore it updates the states from n1+2 to n2-3.
*	(n2 is the first index that does NOT get processed, thus at a
*	boundary n2-1 is the boundary state).
*/

#if defined(__cplusplus)
extern "C" {
#endif /* defined(__cplusplus) */
    LOCAL   LSODE_FUNC  ChemHeatSources;
//    int tag;
#if defined(__cplusplus)
}
#endif /* defined(__cplusplus) */


/*ARGSUSED*/
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
//	printf("In ONED_muscl:\n");

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
		Locstate state = Rect_state(icoords, wave);
		for (i = 0; i < kmax; ++i)
		{
		    if (vmuscl->du[i][j] != 0)
			Flag[i] = YES;
		    Lflag(state)[i] = Flag[i];
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
	    for (i = 0; i < dim; ++i) idirs[i] = iperm[(i+swp_num)%dim];

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
        if(swp_num == 2)
        {
	    float *coords;
	    int ix, iy, iz;
	    ix = icoords[0];
	    iy = icoords[1];
	    iz = icoords[2];
	    for(j = start; j < end; ++j)
	    {
		float T = temperature(vmuscl->vst->state[j]);
		if (T>=200)
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
		            printf("iCoords: (%d, %d, %d)\n", icoords[0], icoords[1], icoords[2]);/*wrong,xiaoxue*/
		    	    printf("Coords: (%f, %f, %f)\n", coords[0], coords[1], coords[2]);/*?*/
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
     
        *neq = NUM_OF_SPECIES + 1;   //num_comp+1
        dens = Y[NUM_OF_SPECIES + 2];
        temp = Y[NUM_OF_SPECIES + 3];
	pres = Y[NUM_OF_SPECIES + 4];
        for (i = 0; i < NUM_OF_SPECIES; ++i)
            YY[i] = Y[i+1]/dens;       
	//TMP_1109
/*//10	if (flag == 1)
	{
	    printf("1. pdens[H] = %e\n", Y[2]);	
	    printf("Y[H] = %e\n", YY[1]);
	}    *///10
             
        ReadThermoData();

	//TMP_XY
        //printf("In ChemHeatSources: array Y: %e %e %e %e %e %e %e %e %e %e %e %e %e\n",Y[0],Y[1],Y[2],Y[3],Y[4],Y[5],Y[6],Y[7],Y[8],Y[9],Y[10],Y[11],Y[12]);
        //fprintf(stdout, "1-before ODE SOLVER %e %e %e %e %e %e %e %e %e %e %e %e %e\n",Y[0],Y[1],Y[2],Y[3],Y[4],Y[5],Y[6],Y[7],Y[8],Y[9],Y[10],Y[11],Y[12]);
	//printf("Density: %e Temperature: %e\n", dens, temp);
          
        //Chemical source term
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

	double *coords;
    	int    *icoords;
    	double GL[3],GU[3],dh[3], uvel[3];
    	for(i = 0; i< 3; i++)
   	{
       	 GL[i] = vmuscl->GL[i];
         GU[i] = vmuscl->GU[i];
         dh[i] = vmuscl->dh[i];
        }
        double 	wall_z[2];
	double 	lwall = GL[2] + 0.1;
        wall_z[0] = GL[2] + 0.8 + 0.1;
        wall_z[1] = GU[2] - 0.1;
        double 	ramp_x[2],ramp_z[2];
        float 	tan10 =  tan(3.1415926*10/180);
        ramp_x[0] = -7.0;
        ramp_x[1] = -7 + 0.8/tan10;
        ramp_z[0] = 1.6;
        ramp_z[1] = 1.6 + 0.8;
        double 	dist_z[2];
        double 	T, P;
        double 	flux[4];
        int    	idir = vmuscl->idir;
        double 	delta;
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
		        if(vmuscl->turb_boundary_layer == YES && idirs[0] == 2 && coords[0] > ramp_x[0] && coords[0] < ramp_x[1] && dist_z[0] > 0 && dist_z[0] <= dh[2] && swp_num != -1)
			{
				float cos10 = cos(3.1415926*10/180);
		    		float sin10 = sin(3.1415926*10/180);
		    		/*Set the physical flux across a wall to be zero here.*/
		    		F0[j] = 0;
	            		F1[j] = 0;
		    		F2[j] = p[j];
		    		F3[j] = 0;
		    		F4[j] = 0;
				/*FIXME*/
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
		    alloc_state(vmuscl->front->interf,&tmpst,
		                vmuscl->front->sizest);
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
		    alloc_state(vmuscl->front->interf,&tmpst,
		                vmuscl->front->sizest);
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
		double u, 
                double v,
		double w, 
		double y0,
		double yw,
		double *flux,
                Locstate statep)	
{      
		double Tmatch = temperature(statep);
		if( Tmatch < Wall_Temp) 
		{
			flux[0] = 0;
			flux[1] = 0;
			flux[2] = 0;
			flux[3] = 0;
			return;
		}

		debug_print("WALL_MODEL","Entered oned_turbulence_boundarylayer_solver()\n");		
  		float delta = fabs(y0 - yw);
		float ymatch = delta;
		//double *Urans, *Trans, *Rrans, *MUrans, *MUTrans, *Krans, *aw, *ap, *ae, *qp, *Uransold, *Transold;//12 variables
		//TK
  		double *Urans, *Trans, *Rrans, *MUrans, *MUTrans, *Krans, *aw, *ap, *ae, *qp;//10 variables
		double *dyrans, *yfrans;
		double ycrans[200];
		double yplusmax = 200;
  		double wallunit= 0.001;
		
		double dyrans0;
		int NY = 1, j, i, k, kk ,iter;
		double rnm1 = 1;
		double lvmin = 1000;
		double Umatch;
		double Tw;
		double muw,kw;
		double tauwold, qwold, converged, damping;
		double intfac,muint,mutint,kint,factor,unit;
		double err_Vel,err_Temp;
		double tauw, rhow, utau, lv, qw;
		double *pmu, *pk;
		double uint;
		double rho0;
		double sigma13, sigma23, sigma33, usigma, q;
		int    normal; //normal direction of the wall
		if( yw > y0) normal = 1;
		if( yw < y0) normal = -1;
	
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

		uni_array(&dyrans,NY+4,sizeof(double));
		uni_array(&yfrans,NY+4,sizeof(double));
		
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

		uni_array(&Urans,NY+4,sizeof(double));
		uni_array(&Trans,NY+4,sizeof(double));
		uni_array(&Rrans,NY+4,sizeof(double));
		uni_array(&MUrans,NY+4,sizeof(double));
		uni_array(&MUTrans,NY+4,sizeof(double));
		uni_array(&Krans,NY+4,sizeof(double));
		uni_array(&aw,NY+4,sizeof(double));
		uni_array(&ap,NY+4,sizeof(double));
		uni_array(&ae,NY+4,sizeof(double));
		uni_array(&qp,NY+4,sizeof(double));
		/*
		uni_array(&Uransold,NY+4,sizeof(double));
		uni_array(&Transold,NY+4,sizeof(double));
		*/
      		Locstate state;
      		state = statep
;
	        rho0 = Dens(state);
		double C_P = C_P(state);

	  	//The actual wall-model --- solve ODEs for each wall face
	  	Umatch = sqrt(u*u + v*v);
	  	Tw = Wall_Temp;

	  	//Initial Guess -- linear profiles	  	
	  	for (j=0; j < NY; j++)
		{
			Urans[j] = Umatch/ymatch*ycrans[j];
		  	Trans[j] = Tw + (Tmatch - Tw)*ycrans[j]/ymatch;
	 	}
		//muw = dynamic_viscosity(state,Tw);
                //kw = dynamic_thermal_conductivity(state,Tw);
		//TK
		dynamic_viscosity_thermalconduct(state,Tw,&muw,&kw);
	  	tauw = muw*Urans[0]/ycrans[0];
	  	rhow = rho0*Tmatch/Tw;
	  	utau = sqrt(tauw/rhow);
	  	lv = muw/(rhow*utau);
	  	qw = - kw*(Trans[0] - Tw)/ycrans[0];

		//Iterate until convergence
	 	//for(i =1;i < MAXSLOOPSTEP;i++)
	 	//{
		iter = 0;
	 	do{
			
			//  Set up auxiliary quantities in all cells in 1D grid
			tauwold = tauw ;
	 		qwold = qw ;
			/*	
          		for (j=0; j < NY; j++)
			{
				Uransold[j] = Urans[j];
		  		Transold[j] = Trans[j];
         		}
			*/
	  		
	  		for (j=0 ; j < NY ; j++) 
			{
				// MUrans[j] = dynamic_viscosity(state,Trans[j]);
                                // Krans[j] = dynamic_thermal_conductivity(state,Trans[j]);
				//TK
				dynamic_viscosity_thermalconduct(state,Trans[j],&MUrans[j],&Krans[j]);
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

	 		
	  		tauw = muw * Urans[0] / ycrans[0] ;
	  		rhow = rho0*Tmatch/ Tw ;
	  		utau = sqrt( tauw / rhow ) ;
	  		lv = muw / ( rhow * utau ) ;
	  		qw = - kw* ( Trans[0] - Tw ) / ycrans[0] ;   // wall heat-flux

			/*	  
			//  Judge whether converged
         		err_Vel = 0;
	  		err_Temp = 0; 

        		for(k = 0; k < NY-1;k++)

			{
      				err_Vel += (Urans[k] - Uransold[k])*(Urans[k] - Uransold[k]);
				err_Temp += (Trans[k] - Transold[k])*(Trans[k] - Transold[k]);
			}

	 		if( (err_Vel < 0.01)&& (err_Temp < 0.01))
				break;
			*/
			// TK
			converged = fabs((tauw - tauwold) / (tauwold + MACH_EPS));

			if( fabs((qw - qwold)/(qwold + MACH_EPS)) > converged )
				converged = fabs( (qw - qwold)/(qwold+MACH_EPS));
	 		iter ++; 
	          }while(converged > 1e-4);

                  	sigma13 = tauw / Umatch * u  * normal;//tauw/Umatch * u is positive
                        sigma23 = tauw / Umatch * v  * normal;
                        sigma33 = muw * 4.0/3.0 * w/delta*normal;//muw*4/3*w/delta will be positive
                        q = -qw * normal  ; //- qw is positive
                    //    vmuscl->vst->wall_unit[index] = lv;
                    //    Vec_Gas_field_set(vmuscl->vst,wall_unit) = YES;
        

                        flux[0] = q;
                        flux[1] = sigma13;
                        flux[2] = sigma23 ;
                        flux[3] = sigma33 ;

			//free_these(12,Urans,Trans,Rrans,MUrans,MUTrans,Krans,aw,ap,ae,qp,Uransold,Transold);
			//TK
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

	float ramp_x[2],ramp_z[2];
	float tan10 = tan(3.1415926*10/180);
	ramp_x[0] = -7.0;
	ramp_x[1] = -7.0 + 0.8/tan10;
	ramp_z[0] = 1.6;
	ramp_z[1] = 1.6 + 0.8;
	float dist_z;

        if (dim != 3) return;

        for (i = 0; i < dim; i++)
        {
            dh[i] = rgr->h[i];
            imin[i] = 0;
            imax[i] = gmax[i];
        }
	
	for (i = 1; i < dim; ++i)
	{
		if (lbuf[i] > 0) imin[i] -= lbuf[i];
		if (ubuf[i] > 0) imax[i] += ubuf[i];
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
		    	dist_z = coords[2] - (ramp_z[0] + tan10*(coords[0] - ramp_x[0]));
			if (dist_z > 0 && dist_z < 1.5*dh[2] && coords[0] > ramp_x[0] && coords[0] < ramp_x[1]&& !is_obstacle_state(state))
				Merge_cell(state) = YES;
			else 	Merge_cell(state) = NO;	
		    }
		}
        }
}

EXPORT  void Update_merge_cell(
	float		dt,
        Wave            *wave)
{
        RECT_GRID       *rgr = wave->rect_grid;
	int		*lbuf = rgr->lbuf;
	int		*ubuf = rgr->ubuf;
        int             dim = rgr->dim;
        int             *gmax = rgr->gmax;
        int             icoords[3],icoords2[3];
        int             imin[3],imax[3];
        int             i, j, k, i0, i1, i2;
        Locstate        state,state2;
        float           *coords,*coords2;
        float           dh[MAXD];
        float           vol;
	float		du,dE,drho,pdens,p;
	float		ramp_x[2],ramp_z[2];
	float		ds0[6],ds1[6],dist[2];
	float 		tan10 = tan(3.1415926*10/180);
	float		dtdh[3];
	float 		T,P,delta;
	ramp_x[0] = -7.0;
        ramp_x[1] = -7 + 0.8/tan10;
	ramp_z[0] = 1.6;
	ramp_z[1] = 1.6 + 0.8;	
        
        //if(vmuscl->turb_boundary_layer =! YES) return;/*FIXME*/
        if (dim != 3)
            return;

        for (i = 0; i < dim; i++)
        {
            dh[i] = rgr->h[i];
            imin[i] = 0;
            imax[i] = gmax[i];
        }
	for( i = 0; i < dim; i++)
	{
	 	if (lbuf[i] > 0)imin[i] -= lbuf[i];
		if (ubuf[i] > 0)imax[i] += ubuf[i];
	}
	for (i = 0; i< 3; i++)  {dtdh[i] = dt/dh[i];}


        /*Boundary layer for the flat part. */
	double y0,yw;
        double  uvel[3],flux[4];
	float   xm,ym;
	int     icoordsp[2][2][3];
        int     icoordsx[2],icoordsz[2];
        float   intpx, intpz;
	Locstate  statep[4];
	float   cos10 =  cos(3.1415926*10/180);
	float   sin10 =  sin(3.1415926*10/180); 
	float   yd;       
	/*upwall*/
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
		    if( coords[2] < 3.9 && coords[2] > 3.9 - dh[2] && sqrt(coords[0]*coords[0] + coords[1]*coords[1]) > 0.1)
		    {
			icoords2[0] = icoords[0];
			icoords2[1] = icoords[1];
			icoords2[2] = icoords[2] - (tk - 1);
			coords2 = Rect_coords(icoords2, wave);
			state2  = Rect_state(icoords2, wave);
			y0 = coords2[2];
			yw = 3.9;
               	    	for (i = 0; i < dim; ++i)
               	            uvel[i] = vel(i,state2);
		        for (i = 0; i < 4 ; i++)
			    flux[i] = 0.0;
		    	if ( fabs(uvel[0]) > MACH_EPS || fabs(uvel[1]) > MACH_EPS)
                    	    oned_turbulence_boundarylayer_solver(dim,uvel[0],uvel[1],uvel[2],y0,yw,flux,state2);
			for( i = 0; i < 3; i++)
				Mom(state)[i] -= flux[i + 1]*dt/dh[2];
			Energy(state) -= flux[0]*dt/dh[2];

			    Boundary_flux(state)[0] = -flux[0];
			    Boundary_flux(state)[1] = -flux[1];
			    Boundary_flux(state)[2] = -flux[2];
			    Boundary_flux(state)[3] = -flux[3];
			    printf("up wall: flux[0] = %f, flux[1] = %f, flux[2] = %f,flux[3] = %f\n", -flux[0],-flux[1],-flux[2],-flux[3]);
		    }
		    else if((coords[0] > ramp_x[1] && coords[2] - 2.4 > 0 && coords[2] - 2.4 < dh[2])||(coords[0] <ramp_x[0] && coords[2] > 1.6 && coords[2] < 1.6 + dh[2]))
		    {
			icoords2[0] = icoords[0];
			icoords2[1] = icoords[1];
			icoords2[2] = icoords[2] + (tk - 1);
			coords2 = Rect_coords(icoords2, wave);
			state2  = Rect_state(icoords2, wave);
			y0 = coords2[2];
			if( coords[0] > ramp_x[1] && coords[2] - 2.4 > 0) yw = 2.4 ;
			else yw = 1.6;
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
			    printf("low wall flux[0] = %f, flux[1] = %f, flux[2] = %f,flux[3] = %f\n", flux[0],flux[1],flux[2],flux[3]);
		    }
		    else if(coords[0] > ramp_x[0] && coords[0] < ramp_x[1] && dist[0] > 0 && coords[2] - (ramp_z[0] + tan10*(coords[0] - ramp_x[0])) >0 && coords[2] - (ramp_z[0] + tan10*(coords[0] - ramp_x[0])) < dh[2])
		    {
			delta = (tk - 0.5)*dh[2];
			yd = ramp_z[0] + tan10*(coords[0] - ramp_x[0]);
			xm = coords[0] - delta*sin10;
			ym = yd + delta*cos10;

			icoordsx[0] = icoords[0] - 1;
			icoordsx[1] = icoords[0];
			if(ym >= coords[2] + (tk-1)*dh[2] && ym < coords[2] + tk*dh[2]) 
			    icoordsz[0] = icoords[2] + tk - 1;
			else if(ym > coords[2] + (tk-2)*dh[2] && ym < coords[2] + (tk-1)*dh[2])
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
			if( intpx < 0 || intpx > 1 || intpz < 0 || intpz > 1) 
				printf("intpx is %f, intpy is %f, xm is %f  ym is %f coords[0] = %f, coords[2] = %f\n", intpx, intpz,xm,ym,coords[0],coords[2]);

			if(intpx < 0.5) i = 0; else i = 1;
			if(intpz < 0.5) j = 0; else j = 1;
			state2  = statep[2*i + j];//FIMEME
			for( i = 0; i < 3; i++)
			uvel[i] = intpz*(intpx*vel(i,statep[0]) + (1 - intpx)*vel(i,statep[2]))
				+(1 - intpz)*(intpx*vel(i,statep[1]) + (1 - intpx)*vel(i,statep[3]));
			
			float u[3];/*velocity which are parallel and normal to the ramp.*/
		    	u[0] = uvel[0]*cos10 + uvel[2]*sin10;
		    	u[1] = uvel[1];
		    	u[2] = uvel[2]*cos10 - uvel[0]*sin10;
			printf("At the ramp, u[0] =%f  ,u[1] = %f ,u[2] = %f\n", u[0],u[1],u[2]);
		    	for (i = 0; i < 4 ; i++)
			flux[i] = 0.0;
			printf("yd = %f ,yd + delta = %f\n",yd,yd + delta);
		    	if ( fabs(u[0]) > MACH_EPS || fabs(u[1]) > MACH_EPS)
                    	    oned_turbulence_boundarylayer_solver(dim,u[0],u[1],u[2],yd + delta,yd,flux,state2);
			printf("after the oned_solver,ramp flux[0] = %f, flux[1] = %f, flux[2] = %f,flux[3] = %f\n", flux[0],flux[1],flux[2],flux[3]);
			float f[3];
			for( i = 0; i < 3; i++)
		        f[0] = (flux[1]*cos10 - flux[3]*sin10)/cos10;
		        f[1] = flux[2]/cos10;
		        f[2] = (flux[1]*sin10 + flux[3]*cos10)/cos10;
			for( i = 0; i < 3; i++)
			Flux(state)[5][1+i] += f[i];
			Flux(state)[5][1] += flux[0];

			    Boundary_flux(state)[0] =  flux[0];
			    Boundary_flux(state)[1] =  flux[1];
			    Boundary_flux(state)[2] =  flux[2];
			    Boundary_flux(state)[3] =  flux[3];
			 
		    }
		}
	    }
	}      

/*store boundary flux in the structure Gas for printing*/
        icoords[0] = imin[0];
        icoords[1] = imin[1];
        for( i2 = 0; i2 < gmax[2]; i2++)
        {
                icoords[2] = i2;
                coords = Rect_coords(icoords,wave);
                if(coords[2] > ramp_z[1] && coords[2] < ramp_z[1] + dh[2]) break;
        }
        int upwall_index = i2;
        icoords[2] = upwall_index;
        coords = Rect_coords(icoords,wave);
        for (i0 = imin[0]; i0 < imax[0]; i0++)
        {
            icoords[0] = i0;
            for (i1 = imin[1]; i1 < imax[1]; i1++)
            {
                icoords[1] = i1;
                for (i2 = imin[2]; i2 < imax[2]; i2++)
                {
                    icoords[2] = i2;
                    coords = Rect_coords(icoords,wave);
                    dist[0]= coords[2] - (ramp_z[0] + tan10*(coords[0] - ramp_x[0]));
                    if(dist[0] > 0 && dist[0] <= dh[2] && coords[0] > ramp_x[0] && coords[0] < ramp_x[1])
                    {
                        state = Rect_state(icoords,wave);
                        icoords2[0] = i0;
                        icoords2[1] = i1;
                        icoords2[2] = upwall_index;
                        state2 = Rect_state(icoords2,wave);
                        for( i = 0; i < 4; i++)
                        Boundary_flux(state2)[i] = Boundary_flux(state)[i];
                    }
                }
            }
        }



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
		    dist[0]= coords[2] - (ramp_z[0] + tan10*(coords[0] - ramp_x[0]));
		    if ( !is_obstacle_state(state) && Merge_cell(state) && flux_flag(state) == 1)
                    {
			dist[0] = coords[2] - (ramp_z[0] + tan10*(coords[0] - ramp_x[0]));
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
			    for( i= 0; i< 6; i++) ds1[i] = 1.0;
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
			    Dens(state)   += drho;
			    Dens(state2)  += drho;						   

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
					pdens(state)[j]  += drho;
					pdens(state2)[j] += drho;
				    } 
                		}
            		    }
			    /*update the momentum.*/
			    p = pressure(state);			    
			    for(j = 0;j < 3;j++ )
			    {
				du = 0;
				for( i = 0; i < 3 ; i++)
			        {
				    du += (Flux(state )[2*i][j+2]*ds0[2*i] - Flux(state )[2*i+1][j+2]*ds0[2*i+1])*dt/dh[i];
				    du += (Flux(state2)[2*i][j+2]*ds1[2*i] - Flux(state2)[2*i+1][j+2]*ds1[2*i+1])*dt/dh[i];
			        }
				if( j == 0)
				du += - p*(ds0[0] - ds0[1])*dt/dh[0];
				du = du/(1 + vol);
				Mom(state )[j] += du;
				Mom(state2)[j] += du;
			    }
			    /*update the energy.*/
		      	    dE = 0;
			    for( i = 0; i < 3 ; i++ )
			    {
				dE += (Flux(state )[2*i][1]*ds0[2*i]- Flux(state )[2*i+1][1]*ds0[2*i+1])*dt/dh[i];
				dE += (Flux(state2)[2*i][1]*ds1[2*i]- Flux(state2)[2*i+1][1]*ds1[2*i+1])*dt/dh[i];
			    }
			    dE = dE/(vol + 1);
			    Energy(state)   += dE  ;
			    Energy(state2)  += dE  ; 
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
			    /*The second cell.*/
			    bool l = 0;
			    if( i0 -1 > imin[0] && i2 - 1 > imin[2])
			    {
			    icoords2[0] = i0 - 1;
			    icoords2[1] = i1;
			    icoords2[2] = i2 - 1;
			    state2 = Rect_state(icoords2,wave);
			    if (!is_obstacle_state(state2) && flux_flag(state2)==1) l = 1;
			    }
			
                            //Update the merged cells by using Flux(state)[][] and then split.
		            if ( flux_flag(state) == 1 && l == 0)
			    {
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
				{
					dE += (Flux(state)[2*i][1]*ds0[2*i]- Flux(state)[2*i+1][1]*ds0[2*i+1])*dt/dh[i]; 
				}
			    	dE = dE/vol;
			    	Energy(state) += dE;

			    	/*update the momentum in each cell.*/
				p = pressure(state);
			    	for(j =0;j<3;j++)
			    	{
				    du = 0;
				    if(j == 0) du += -p*(ds0[0] - ds0[1])*dt/dh[0];
				    for( i = 0; i < 3; i++)
				    du += (Flux(state)[2*i][2+j]*ds0[2*i]- Flux(state)[2*i+1][2+j]*ds0[2*i+1])*dt/dh[i]; 
				    du = du/vol;
				    Mom(state)[j] += du;
			    	}
			    }
			    else if(flux_flag(state) == 1 && l == 1 )
			    {
			    	/*update the density in each cell.*/
			    	drho = 0;
			        drho += (Flux(state)[0][0] - Flux(state)[1][0]*ds0[1])*dt/dh[0]; 
				drho += (ds0[0] - 1)*Flux(state2)[1][0]*dt/dh[0];
			    	for( i = 1; i < 3; i++)
				{
			        	drho += (Flux(state)[2*i][0]*ds0[2*i]- Flux(state)[2*i+1][0]*ds0[2*i+1])*dt/dh[i]; 
				}
			    	drho = drho/vol;
			    	Dens(state) += drho ;				
			    	/*update the partial density.*/
			    	if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            		    	{
                		    Gas_param *params = Params(state);
                		    int    num_comps;
                		    if((params != NULL) && ((num_comps = params->n_comps) != 1))
                		    {
                		   	for(j = 0; j < num_comps; j ++)
				    	{
					    drho = 0;
				            drho += (Flux(state)[0][5+j]-Flux(state)[1][5+j]*ds0[1])*dt/dh[0]; 
					    drho +=  Flux(state2)[1][5+j]*(ds0[0] -1)*dt/dh[0]; 	
					    for( i = 1; i < 3; i ++)
					    {
						drho += (Flux(state)[2*i][5+j]*ds0[2*i]
							- Flux(state)[2*i+1][5+j]*ds0[2*i+1])*dt/dh[i]; 
					    }
					    drho = drho/vol;
					    pdens(state)[j]  += drho;
				   	} 
				    }
                		}
        
			    	/*update the energy in each cell.*/
			    	dE = 0;
			        dE += (Flux(state)[0][1] - Flux(state)[1][1]*ds0[1])*dt/dh[0]; 
				dE += (ds0[0] - 1)*Flux(state2)[1][1]*dt/dh[0];	
			    	for( i = 1; i < 3; i++)
				{
					dE += (Flux(state)[2*i][1]*ds0[2*i]  - Flux(state)[2*i+1][1]*ds0[2*i+1])*dt/dh[i]; 
				}
			    	dE = dE/vol;
			    	Energy(state) += dE  ;
			    	/*update the momentum in each cell.*/
				p = pressure(state);
			    	for(j = 0;j < 3;j++)
			    	{
					du = 0;
					du += (Flux(state )[0][2 + j] - Flux(state)[1][2 + j]*ds0[1])*dt/dh[0];
					du +=  Flux(state2)[1][2 + j]*(ds0[0]-1)*dt/dh[0];
					if( j == 0) du += - p*(ds0[0] - ds0[1])*dt/dh[0];
					for( i = 1; i < 3; i++)
					{
				   	du += (Flux(state)[2*i][2+j]*ds0[2*i]- Flux(state)[2*i+1][2+j]*ds0[2*i+1])*dt/dh[i];
					}
					du = du/vol;
					Mom(state)[j] += du;
			   	}
			   }
			}

                    }/*end of the update of merged cell.*/
                }
	    }
	}/*end of 3 space loop*/

}

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
	float random_numbers[3][My+2*NFy][Mz+2*NFz];

        for(d=0;d<dim;++d)
        {
            short unsigned int seed[2];
            for(j=0; j<My+2*NFy; ++j)
              for(k=0; k<Mz+2*NFz; ++k)
                random_numbers[d][j][k] = random_gaussian(0,1,seed);
        }

        float	 S11, S22, S33, intensity=0.016; //goes to input 
	float  	 u_x=180.0, u_y=0.0, u_z=0.0;
	Locstate st;
	int   	 st_type, icoords[3];
	float 	 *coords;
	float    GLx = front->rect_grid->GL[0];

        float 	u_field_old[3][My][Mz], u_field[3][My][Mz], u_perturb[3][My][Mz];	
        float 	dt = front->dt; // time step	
        float 	tau_t = Ix/u_x; // Lagrangian time step
        float 	filter_exp1 = sqrt(1-exp(-PI*dt/tau_t));
        float 	filter_exp2 = exp(-PI*dt/(2*tau_t));
	static  int countt = 0; 
        if(countt == 0)
        {	
         for(d=0; d<dim; ++d)
           for(j=0; j<My; ++j)
	     for(k=0; k<Mz; ++k)
	     {	
	      icoords[0] = 0;
	      icoords[1] = j;
	      icoords[2] = k;

	      st = Rect_state(icoords, wave);
	      Vel_Field(st)[d] = 0.0; 	
	     }
	  countt ++;
        }	

        for(d=0; d<dim; ++d)
          for(j=0; j<My; ++j)
            for(k=0; k<Mz; ++k)
	    {
	    icoords[0] = 0;
	    icoords[1] = j;
	    icoords[2] = k;

	    coords = Rect_coords(icoords, wave);
	    st = Rect_state(icoords, wave);

              for(m=0; m<2*NFy+1; ++m)
                 for(n=0; n<2*NFz+1; ++n)
           	    u_field[d][j][k] = b_y[m]*b_z[n]*random_numbers[d][j+m][k+n];

//	      u_perturb[d][j][k] = u_field_old[d][j][k]*filter_exp2 + u_field[d][j][k]* filter_exp1;
	      u_perturb[d][j][k] = Vel_Field(st)[d]*filter_exp2 + u_field[d][j][k]* filter_exp1;
	    }	

	//Prescribed Reynolds stress tensor (isotropic)
	S11= S22 = S33 = (intensity)*(intensity)*u_x*u_x/3; 

	for(j=0; j<My; ++j)
	 for(k=0; k<Mz; ++k)
	  {
	    icoords[0] = 0;
	    icoords[1] = j;
	    icoords[2] = k;

	    coords = Rect_coords(icoords, wave);
	    st = Rect_state(icoords, wave);

	    if(coords[0] < GLx+dh)
	    {
	      if(!is_obstacle_state(st))
	      {
		st_type = state_type(st);
		set_state(st, TGAS_STATE, st);

	    	Vel(st)[0] = u_x + sqrt(S11)*u_perturb[0][j][k];
	    	Vel(st)[1] = u_y + sqrt(S22)*u_perturb[1][j][k];
	    	Vel(st)[2] = u_z + sqrt(S33)*u_perturb[2][j][k];

	     	Vel_Field(st)[0] = u_field[0][j][k]; 	
		Vel_Field(st)[1] = u_field[1][j][k]; 	
		Vel_Field(st)[2] = u_field[2][j][k]; 	

		set_state(st, st_type, st);
	      }
	    }
	  }

} 

  

#endif /* defined(DEBUG_MUSCL) */
