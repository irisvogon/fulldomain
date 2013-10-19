/*
*
*				gcurve.c:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Propagation algorithm for curve hypersurface boundaries in three
*	space dimensions.
*/

#if defined(THREED)

#include <gdecs/gdecs.h>



struct _O_BOND
{
  	BOND 		*old_bond;
 	BOND_TRI 	**_btris;
  	float		*angle;
  	int		*orientation; /*TODO change to type ORIENTATION */
  	int		num_btris;
};
typedef struct _O_BOND O_BOND;

LOCAL  	O_BOND   *make_o_bond(BOND*);
LOCAL  	int      attached_b_curve_propagate(Front*,POINTER,CURVE*,CURVE*,float);
LOCAL   SURFACE  *find_physical_surface_at_curve(CURVE*, ORIENTATION*);
LOCAL   SURFACE  *adjacent_surface(SURFACE*, ORIENTATION, CURVE*, 
                                   ANGLE_DIRECTION, ORIENTATION*);
LOCAL   bool     matchable_comps(COMPONENT , COMPONENT,INTERFACE*);
LOCAL   void     reset_fixed_btri_states(BOND*, O_BOND*, Front*);
LOCAL   int      fixed_curve_propagate(Front*,POINTER,CURVE*,CURVE*,float);
LOCAL   int      subdomain_curve_propagate(Front*,POINTER,CURVE*,CURVE*,float);
LOCAL	void	 difference3d(float*,float*,float*);
LOCAL	void	 vector_scale3d(float*,float);
LOCAL   int      fixed_curve_debug = NO;
LOCAL   void     set_curve_vel_to_zero(CURVE*);

//#bjet2
#define surf_ang_oriented_l_to_r(ang_dir,orient)                             \
        ((((ang_dir) == CLOCKWISE && (orient) == POSITIVE_ORIENTATION) ||     \
        ((ang_dir) == COUNTER_CLOCK && (orient) == NEGATIVE_ORIENTATION)) ?  \
        YES : NO)



// #bjet
EXPORT  void	g_curve_propagate_3d(
	Front		*front,
	POINTER 	wave,
	CURVE		*oldc,
	CURVE		*newc,
	float		dt)
{
	int     status;
		
	DEBUG_ENTER(g_curve_propagate_3d)
        debug_print("curve_propagate","Entered g_curve_propagate_3d()\n");
	
	/* Propagate curve according to its type */

        switch (curve_type(newc))
	{
	    case PASSIVE_CURVE:   // to do
		 status = GOOD_CURVE;
		 break;
    	    case DIRICHLET_CURVE:
            case FIXED_CURVE:
                status = fixed_curve_propagate(front,wave,oldc,newc,dt);
                break; 
 	    case ATTACHED_B_CURVE:
	        status = attached_b_curve_propagate(front,wave,oldc,newc,dt);
	        break;
	    case NEUMANN_CURVE_W:
	    case NEUMANN_CURVE_P:
	        status = attached_b_curve_propagate(front,wave,oldc,newc,dt);
	        break;
	    case SUBDOMAIN_HSBDRY:
	        status = subdomain_curve_propagate(front,wave,oldc,newc,dt);
	        break;
	    default:
	        (void) printf("ERROR: g_curve_propagate_3d() oldc = %p ",oldc); 
	        print_hsbdry_type("Unable to process HSBDRY type ",
			      curve_type(oldc),"\n",oldc->interface);
		clean_up(ERROR);
	        break;
	}

        debug_print("curve_propagate","Left g_curve_propagate(), \n");
	DEBUG_LEAVE(g_curve_propagate_3d)
	return;
} 		/*end g_curve_propagate_3d*/

LOCAL   int     attached_b_curve_propagate(
        Front           *fr,
        POINTER         wave,
        CURVE           *oldc,
        CURVE           *newc,
        float           dt)
{
        O_SURFACE       Si, NewSi;
	O_SURFACE       Sa, NewSa;
	O_SURFACE       Sb, NewSb;
	BOND            *oldb, *newb;
        POINT           *oldp, *newp;
        Locstate        ahead, behind;
        Locstate        new_ahead, new_behind;
        float           V[MAXD];
        int             dim = fr->interf->dim;
        int             status;

        debug_print("curve_propagate","Entered attached_b_curve_propagate(), \n");
        
        oldb = oldc->first;
        newb = newc->first;

	Si.surface = find_physical_surface_at_curve(oldc,&Si.orient);
	NewSi.surface = find_physical_surface_at_curve(newc,&NewSi.orient);
	
        if (Si.surface != NULL &&  wave_type(Si.surface) >= FIRST_VECTOR_PHYSICS_WAVE_TYPE)
        {    
             while(oldb)
             {
                 if ((oldb != oldc->last) && (!n_pt_propagated(newb->end)))
                 {
                     point_propagate(fr,wave,oldb->end,newb->end,oldb->next,
                                     oldc,dt,V);
                     n_pt_propagated(newb->end) = YES;
                     if(wave_type(oldc) == FORWARD_SHOCK_WAVE) //????????
                     {
                         ahead = right_state(newb->end);
                         behind = left_state(newb->end);
                     }
                     else
                     {
                         ahead = left_state(newb->end);
                         behind = right_state(newb->end);
                     }
                 }
                 //if (fr->bond_propagate != NULL)
                 //   (*fr->bond_propagate)(fr,wave,oldb,newb,oldc,dt);
                 else
                     set_bond_length(newb,dim); /* Update new bond length */
                 if (oldb == oldc->last)
                     break;
                 oldb = oldb->next;
                 newb = newb->next;
             }
             if (Si.orient == POSITIVE_ORIENTATION)
             {
                 Sa.surface = adjacent_surface(Si.surface,Si.orient,oldc,
                                               CLOCKWISE,&Sa.orient);
                 Sb.surface = adjacent_surface(Si.surface,Si.orient,oldc,
                                               COUNTER_CLOCK,&Sb.orient);
             }  
             else
             {
                 Sa.surface = adjacent_surface(Si.surface,Si.orient,oldc,
                                           COUNTER_CLOCK,&Sa.orient);
                 Sb.surface = adjacent_surface(Si.surface,Si.orient,oldc,
                                           CLOCKWISE,&Sb.orient);
             }
             
             if (NewSi.orient == POSITIVE_ORIENTATION)
             {
                 NewSa.surface = adjacent_surface(NewSi.surface,NewSi.orient,
                                             newc, CLOCKWISE,&NewSa.orient);
                 NewSb.surface = adjacent_surface(NewSi.surface,NewSi.orient,
                                             newc, COUNTER_CLOCK,&NewSb.orient);
             }
             else
             {
                 NewSa.surface = adjacent_surface(NewSi.surface,NewSi.orient,
                                             newc, COUNTER_CLOCK,&NewSa.orient);
                 NewSb.surface = adjacent_surface(NewSi.surface,NewSi.orient,
                                             newc, CLOCKWISE,&NewSb.orient);
             }
        }
        else
        {
            status = fixed_curve_propagate(fr,wave,oldc,newc,dt);
            return status;
        }
        debug_print("curve_propagate","Left attached_b_curve_propagate(), \n");

}

//#bjet2  find the btri according to comp, only right for jet3d
LOCAL  BOND_TRI * find_correspond_btri(BOND_TRI *btri, BOND *bond)
{
BOND_TRI    **newbtri;
SURFACE     *newsurf, *surf;

	surf = btri->surface;
	for(newbtri = Btris(bond); newbtri && *newbtri; newbtri++)
	{
	    newsurf = (*newbtri)->surface;
	    if(positive_component(newsurf) == positive_component(surf) &&
	       negative_component(newsurf) == negative_component(surf))
	       return *newbtri;
	}
	return NULL;
}

//set the intersection states between two btris at the end of a bond 
LOCAL  void reset_fixed_btri_states(
	BOND 	*bond, 
	O_BOND  *ordered_bond, 
	Front   *fr)
{
int       	i, j, k;
SURFACE   	*surf1, *surf2;
BOND_TRI  	*bt1, *bt2;
Locstate  	state1, state2;
static Locstate tmpst = NULL;
int	  	sizest = fr->sizest;
POINT		*pos;

	if(tmpst == NULL)
	    alloc_state(fr->interf, &tmpst, sizest);

	for (i = 0; i < ordered_bond->num_btris; i++)
	    for(k=0; k<1; k++)
	    {
	        bt1 = find_correspond_btri(ordered_bond->_btris[i], bond);
            	surf1 = bt1->surface;
                
	    	j = (i == ordered_bond->num_btris - 1) ? 0 : i+1;
            	bt2 = find_correspond_btri(ordered_bond->_btris[j], bond);
            	surf2 = bt2->surface;

                if (wave_type(surf1) == SUBDOMAIN_BOUNDARY ||
                    wave_type(surf2) == SUBDOMAIN_BOUNDARY)
                    continue;
	
		if (wave_type(surf1) == PASSIVE_BOUNDARY ||
                    wave_type(surf2) == PASSIVE_BOUNDARY)
                    continue;
		
		//only consider the end states and ASSUME one pos_comp and one neg_comp
		//uniquly determine one surface.
		if(k == 0)
		{
		    if(positive_component(surf1) == positive_component(surf2))
		    {
		        state1 = right_end_btri_state(bt1);
		        state2 = right_end_btri_state(bt2);
		    }
		    else  if(positive_component(surf1) == negative_component(surf2))
		    {
		        state1 = right_end_btri_state(bt1);
		        state2 = left_end_btri_state(bt2);
		    }
		    else  if(negative_component(surf1) == negative_component(surf2))
		    {
		        state1 = left_end_btri_state(bt1);
		        state2 = left_end_btri_state(bt2);
		    }
		    else  if(negative_component(surf1) == positive_component(surf2))
		    {
		        state1 = left_end_btri_state(bt1);
		        state2 = right_end_btri_state(bt2);
		    }
		    else
		    {
		        printf("ERROR reset_fixed_btri_state, inconsistent two surfaces.\n");
			clean_up(ERROR);
		    }

		    pos = bond->end;
		}
               
		if(wave_type(surf1) == NEUMANN_BOUNDARY &&
                   wave_type(surf2) == DIRICHLET_BOUNDARY)
                {
                    ft_assign(state1,state2,sizest);
                }
                else if (wave_type(surf2) == NEUMANN_BOUNDARY &&
                     wave_type(surf1) == DIRICHLET_BOUNDARY)
                {
                    ft_assign(state2,state1,sizest);
                }
                else if ((wave_type(surf1)  >= FIRST_PHYSICS_WAVE_TYPE) &&
                    (wave_type(surf2) <  FIRST_PHYSICS_WAVE_TYPE))
                {
                    ft_assign(state2,state1,sizest);
                }
                else if ((wave_type(surf2) >= FIRST_PHYSICS_WAVE_TYPE) &&
                         (wave_type(surf1)  <  FIRST_PHYSICS_WAVE_TYPE))
                {
                    ft_assign(state1,state2,sizest);
                }
                else
                {
		    //printf("#wtype %d  %d\n", wave_type(surf1), wave_type(surf2));
		    interpolate_states(fr,0.5,0.5,Coords(pos),state1,
                                   Coords(pos),state2,tmpst);
		    ft_assign(state1,tmpst,sizest);
		    ft_assign(state2,tmpst,sizest);
		}
            }
}

float	noz2d[2800][20][4];

void	get_perturb()
{
	int	i, j, comp;
	float	x, step, dens, pres, u, v;
	char	fname[200];
	FILE	*fp;
	
	//sprintf(fname, "%snozzle2d", get_directory());
	sprintf(fname, "nozzle2d");
	fp = fopen(fname, "r");

	for(i=0; i<2800; i++)
	    for(j=0; j<20; j++)
	    {
		fscanf(fp, "%lf%lf%d%lf%lf%lf%lf", 
			&x, &step, &comp, &dens, &pres, &u, &v);

		noz2d[i][j][0] = dens;
		noz2d[i][j][1] = pres;
		noz2d[i][j][2] = u;
		noz2d[i][j][3] = v;
	    }

	fclose(fp);
}

LOCAL  	int	Nx, Ny, Nz, My, Mz, NNX, XST;
LOCAL	float	****Ra, *bx, *by, *bz, ***pvel, RM, TM;

void	init_vel_mesh()
{
	My = 40;	//mesh in y direction
	//Mz = 60;	//mesh in z direction
	Mz = 1;		//mesh in z direction
	//Ny = My;	//mesh of turb length scale in y
	Ny = 20;	//mesh of turb length scale in y
	//Nx = 13*Ny;	//mesh of turb length scale in x or t direction
	//Nx = 200;	//mesh of turb length scale in x or t direction
	Nx = 50;	//mesh of turb length scale in x or t direction
	Nz = Mz;	//mesh of turb length scale in z
	RM = 0.01 + 0.01*0.07;   //Radius of pertb
	TM = 90.0;		 //angle rangle of pertb
}

void	cylic_shift_Ra(
	int	M1,
	int	M2)
{
	int		MN = 3*M1*M2, i, j, k, n;
	static  float	*ras = NULL, *rar;
	int		nproc, mynode, dst, src;
	
	DEBUG_ENTER(cylic_shift_Ra)

	nproc = pp_numnodes();
	mynode = pp_mynode();

	if(ras == NULL)
	{
	    uni_array(&ras, MN, FLOAT);
	    uni_array(&rar, MN, FLOAT);
	}
	
	for(i=0; i<M1; i++)
	    for(j=0; j<M2; j++)
		for(k=0; k<3; k++)
		    ras[i*M2*3+j*3+k] = Ra[k][0][i][j];

	//send the first layer of Ra to its previous proc
	dst = (mynode + nproc - 1)%nproc;
	pp_send(mynode, ras, MN*FLOAT, dst);
	
	//recv the last layer of Ra from next proc
	src = (mynode + nproc + 1)%nproc;
	pp_recv(src, src, rar, MN*FLOAT);
	
	pp_gsync();
	
	//shift Ra to left
	for(n=0; n<NNX-1; n++)
	    for(i=0; i<M1; i++)
		for(j=0; j<M2; j++)
		    for(k=0; k<3; k++)
			Ra[k][n][i][j] = Ra[k][n+1][i][j];

	//fill the last layer of Ra
	for(i=0; i<M1; i++)
	    for(j=0; j<M2; j++)
		for(k=0; k<3; k++)
		    Ra[k][NNX-1][i][j] = rar[i*M2*3+j*3+k];

	DEBUG_LEAVE(cylic_shift_Ra)
}

void	init_Ra_layer(
	int	M1,
	int	M2)
{
	int	n, i, j;
	int	nproc = pp_numnodes();
	
	//fill the last layer with random numbers.
	if(pp_mynode() == nproc-1)
	{
	    for(n=0; n<3; n++)
		for(i=0; i<M1; i++)
		    for(j=0; j<M2; j++)
		    {
			Ra[n][NNX-1][i][j] = rand()/((float)RAND_MAX + 1)-0.5;
		    }
	}
}

void	init_pert()
{
	int	i, j, k, n;
	int	mynode, nproc, nrx, nnx, M1, M2;
	float	s;

	nproc = pp_numnodes();
	mynode = pp_mynode();
	nnx = 2*Nx+1;
	NNX = nnx/nproc;
	nrx = nnx%nproc;
	
	if(mynode < nrx)
	{
	    NNX++;
	    XST = NNX*mynode;
	}
	else
	{
	    XST = NNX*mynode + nrx;
	}

	M1 = My + 2*Ny;
	M2 = Mz + 2*Nz;

	quad_array(&Ra, 3, NNX, M1, M2, FLOAT);
	tri_array(&pvel, 3, My, Mz, FLOAT);
	uni_array(&bx, 2*Nx+1, FLOAT);
	uni_array(&by, 2*Ny+1, FLOAT);
	uni_array(&bz, 2*Nz+1, FLOAT);

	printf("#init_pert init Ra %d %d\n", NNX, XST);
	
	for(k=0; k<nnx; k++)
	{
	    cylic_shift_Ra(M1, M2);
	    init_Ra_layer(M1, M2);
	}

	printf("#init_pert init bxyz\n");

	s = 0.0;
	for(i=-Nx; i<=Nx; i++)
	{
	    bx[i+Nx] = exp(-PI*i*i/2.0/Nx/Nx);
	    s += bx[i+Nx]*bx[i+Nx];
	}
	for(i=-Nx; i<=Nx; i++)
	    bx[i+Nx] /= sqrt(s);
	
	s = 0.0;
	for(i=-Ny; i<=Ny; i++)
	{
	    by[i+Ny] = exp(-PI*i*i/2.0/Ny/Ny);
	    s += by[i+Ny]*by[i+Ny];
	}
	for(i=-Ny; i<=Ny; i++)
	    by[i+Ny] /= sqrt(s);
	
	s = 0.0;
	for(i=-Nz; i<=Nz; i++)
	{
	    bz[i+Nz] = exp(-PI*i*i/2.0/Nz/Nz);
	    s += bz[i+Nz]*bz[i+Nz];
	}
	for(i=-Nz; i<=Nz; i++)
	    bz[i+Nz] /= sqrt(s);

	printf("#init_pert finish\n");

}

float	get_one_pos_vel(
	int	n,
	int	jj,
	int	kk)
{
	int	i, j, k;
	float	s;

if(debugging("turb2d"))
{
	s = 0.0;
	for(i=0; i<NNX; i++)
	    for(j=-Ny; j<=Ny; j++)
	    {
		s += bx[XST+i]*by[j+Ny]*Ra[n][i][j+Ny+jj][0];
	    }
}
else
{
	s = 0.0;
	for(i=0; i<NNX; i++)
	    for(j=-Ny; j<=Ny; j++)
		for(k=-Nz; k<=Nz; k++)
		{
		    s += bx[XST+i]*by[j+Ny]*bz[k+Nz]*Ra[n][i][j+Ny+jj][k+Nz+kk];
		}
}

	return s;
}

void	get_pert_vel()
{
	int		n, i, j, k, ind;
	static	float	*pvela = NULL;
	int		LEN;

	LEN = 3*My*Mz;
	if(pvela == NULL)
	    uni_array(&pvela, LEN, FLOAT);
	
	printf("#get_pert_vel bf\n");

	for(n=0; n<3; n++)
	    for(i=0; i<My; i++)
		for(j=0; j<Mz; j++)
		{
		    ind = n*My*Mz + i*Mz + j;
		    pvela[ind] = get_one_pos_vel(n, i, j);
		}

	pp_global_sum(pvela, LEN);
	printf("#get_pert_vel af\n");
	
	for(n=0; n<3; n++)
	    for(i=0; i<My; i++)
		for(j=0; j<Mz; j++)
		{
		    ind = n*My*Mz + i*Mz + j;
		    pvel[n][i][j] = pvela[ind];
		}

}

void show_pert_vel();
void show_pert_vel()
{
	int	n, i, j;
	FILE	*fp;

	fp = fopen("pvel", "w");
	for(n=0; n<3; n++)
	{
	    for(i=0; i<My; i++)
	    {
		for(j=0; j<Mz; j++)
		    fprintf(fp, "%15.8e  ", pvel[n][i][j]);
	        fprintf(fp, "\n");
	    }
	    fprintf(fp, "\n");
	}
	fclose(fp);
}

void	init_turb_vel();

void test_pert_vel();
void test_pert_vel()
{
	int	n, i, j;
	int	M1, M2, LEN;
	FILE	*fp;
	char	fname[500];

	init_turb_vel();
	
	M1 = My + 2*Ny;
	M2 = Mz + 2*Nz;
	LEN = 2*Nx;

	sprintf(fname, "pvel_%s.plt", right_flush(pp_mynode(),PP_NODE_FIELD_WIDTH));
	
	fp = fopen(fname, "w");
	fprintf(fp,"TITLE = \"pert vel\"\n"
		   "VARIABLES = \"r\", \"z\", " 
			    "\"u\", \"v\", \"w\"\n");
	fprintf(fp, "ZONE i=%d, j=%d\n", My, LEN);

	for(i=0; i<LEN; i++)
	{
	    get_pert_vel();
	
	    printf("#test_pert %d\n", i);
	    fflush(NULL);

	    for(j=0; j<My; j++)
	    {
		//fprintf(fp, "%15.8e %15.8e   %15.8e %15.8e %15.8e", j*0.01/My, i*0.01/Nx, 
		//    pvel[0][j][0], pvel[1][j][0], pvel[2][j][0]);
		fprintf(fp, "%15.8e %15.8e   %15.8e %15.8e %15.8e", j*1.0, i*1.0, 
		    pvel[0][j][0], pvel[1][j][0], pvel[2][j][0]);
	        fprintf(fp, "\n");
	    }

	    cylic_shift_Ra(M1, M2);
	    
	    init_Ra_layer(M1, M2);
	}
	
	fclose(fp);
}

float	out_pert_vel(
	int	num,
	int	len,
	int	shift)
{
	int	n, i, j;
	int	M1, M2, LEN;
	FILE	*fp;
	char	fname[500];
	float	s;
	
	M1 = My + 2*Ny;
	M2 = Mz + 2*Nz;
	LEN = len*Nx;

	sprintf(fname, "pvel_%s.plt", right_flush(pp_mynode(),PP_NODE_FIELD_WIDTH));
	if(pp_mynode() == 0)
	    sprintf(fname, "pvelin_%d.plt", num);

	fp = fopen(fname, "w");
	if(!fp)
	{
	    screen("can't find %s.",fname);
	    clean_up(-1);
	}
	fprintf(fp,"TITLE = \"pert vel\"\n"
		   "VARIABLES = \"r\", \"z\", " 
			    "\"u\", \"v\", \"w\"\n");
	fprintf(fp, "ZONE i=%d, j=%d\n", My, LEN);

	s = 0.0;
	for(i=0; i<LEN; i++)
	{
	    get_pert_vel();
	
	    printf("#test_pert %d\n", i);
	    fflush(NULL);

	    for(j=0; j<My; j++)
	    {
		fprintf(fp, "%15.8e %15.8e   %15.8e %15.8e %15.8e", j*0.01/My, i*1.0 + shift, 
		    pvel[0][j][0], pvel[1][j][0], pvel[2][j][0]);

		s += sqr(pvel[0][j][0]) + sqr(pvel[2][j][0]);
		
		//fprintf(fp, "%15.8e %15.8e   %15.8e %15.8e %15.8e", j*0.01/My, i*0.01/Nx, 
		//    pvel[0][j][0], pvel[1][j][0], pvel[2][j][0]);
	        fprintf(fp, "\n");
	    }

	    cylic_shift_Ra(M1, M2);
	    
	    init_Ra_layer(M1, M2);
	}
	
	fclose(fp);

	return s;
}

void	generate_pert_vel();

void	generate_pert_vel()
{
	int	i, num, nper;
	float	s, s0;

	init_turb_vel();
	nper = 20;    //number of Nx in each file
	num = 5;    //number of files

	s = 0.0;
	for(i=0; i<num; i++)
	{
	    printf("#out_pert_vel %d\n", i);
	    fflush(NULL);

	    s0 = out_pert_vel(i, nper, nper*Nx*i);
	    
	    printf("#s0 = %15.8e\n", sqrt(s0/(nper*Nx*My*3.0)));
	    s += s0;
	}
	
	printf("#s = %15.8e\n", sqrt(s/(num*nper*Nx*My*3.0)));
}

void	get_pert_vel();
void	init_turb_vel();

void	init_turb_vel()
{
	init_vel_mesh();
	init_pert();
}

void	interp_turb_vel(
	float	*vel,
	float	*coords)
{
	float	hr, ht, r, theta, fr, ft, u, v;
	int	ra, ta, i;

	hr = RM/My;
	ht = TM/Mz;

	vel[0] = 0.0;
	vel[1] = 0.0;
	vel[2] = 0.0;
	r = sqrt(sqr(coords[0]) + sqr(coords[1]));
	theta = angle(fabs(coords[0]), fabs(coords[1]));

	ra = (int)(r/hr);
	fr = r/hr - ra;
	if(ra >= My-1)
	{
	    ra = My-2;
	    fr = 1.0;
	}
	ta = (int)(theta/ht);
	ft = theta/ht - ta;
	if(ta >= Mz-1)
	{
	    ta = Mz-2;
	    ft = 1.0;
	}

	for(i=0; i<3; i++)
	{
	    vel[i] = (pvel[i][ra][ta]*(1.0-fr) + pvel[i][ra+1][ta]*fr)*(1.0-ft)
		    +(pvel[i][ra][ta+1]*(1.0-fr) + pvel[i][ra+1][ta+1]*fr)*ft;

	}
	
	u = vel[0]*cos(theta) - vel[1]*sin(theta);
	v = vel[0]*sin(theta) + vel[1]*cos(theta);
	vel[0] = u;
	vel[1] = v;

}

int  get_char(FILE *fp, int n, char c)
{
int    cnt;
char   ch;

        cnt = 0;
        while(!feof(fp))
        {
            fscanf(fp, "%c", &ch);
            if(ch == c)
                cnt++;

            if(cnt == n)
                break;
        }

        return cnt;
}

void	read_pert_file_header(
	FILE		*fp,
	int		*M,
	int		*N)
{
        if(3 != get_char(fp, 3, '='))
        {
            printf("ERROR read_pert_file_header, file read error =\n");
	    clean_up(ERROR);
        }

	//radius dimension
	fscanf(fp, "%d", M);
	
	if(1 != get_char(fp, 1, '='))
        {
            printf("ERROR read_pert_file_header file read error 1 =\n");
	    clean_up(ERROR);
        }

	//time dimension
	fscanf(fp, "%d", N);
}

void	init_pert_from_file();
void	read_pert_from_file(int, int);
bool	read_pert_cond(int);

LOCAL	float	***pertin = NULL;
LOCAL	float	***pertint = NULL;

void	init_pert_from_file()
{
	FILE	*fp;

	fp = fopen("pvelin_0.plt", "r");
	if(fp == NULL)
	{
	    printf("ERROR init_pvel_from_file. Can not find pvelin_0.plt.\n");
	    clean_up(ERROR);
	}
	
	read_pert_file_header(fp, &My, &Nx);

	printf("#init_pert_from_file %d %d\n", My, Nx);

	tri_array(&pertin, Nx, My, 2, FLOAT);
	tri_array(&pertint, Nx, My, 2, FLOAT);

	fclose(fp);
}

void	cp_into_pertint();

void	cp_into_pertint()
{
	int	i, j;

	for(i=0; i<Nx; i++)
	    for(j=0; j<My; j++)
	    {
		pertint[i][j][0] = pertin[i][j][0];
		pertint[i][j][1] = pertin[i][j][1];
	    }
}

void	read_pert_from_file(
	int		step,
	int		flag)
{
	FILE	*fp;
	int	M, N, i, j;
	float	u, v, w, r, theta;
	char	fname[500];

	if(pertin == NULL)
	{
	    printf("ERROR read_pert_from_file, pertin is NULL.\n");
	    clean_up(ERROR);
	}

	if(flag == 1)
	    step += Nx;

	sprintf(fname, "pvelin_%d.plt", step/Nx);
	fp = fopen(fname, "r");
	if(fp == NULL)
	{
	    printf("ERROR read_pert_from_file. Can not find %s.\n", fname);
	    clean_up(ERROR);
	}
	
	read_pert_file_header(fp, &M, &N);
	
	if(M != My || N != Nx)
	{
	    printf("ERROR read_pert_from_file. Number inconsistent.\n");
	    printf("My %d %d, Nx %d %d\n", My, M, Nx, N);
	    clean_up(ERROR);
	}

	for(i=0; i<Nx; i++)
	    for(j=0; j<My; j++)
	    {
		fscanf(fp, "%lf%lf%lf%lf%lf", &r, &theta, &u, &v, &w);
		
		if(flag == 0)
		{
		    pertin[i][j][0] = u;
		    pertin[i][j][1] = w;
		}
		else
		{
		    pertint[i][j][0] = u;
		    pertint[i][j][1] = w;
		}
	    }

	fclose(fp);
}

bool	read_pert_cond(
	int	step)
{
	if(step%Nx == 0)
	    return YES;

	return NO;
}

void	init_period_pertb(Front*, Wave*);
void	init_period_pertb(
	Front	*front,
	Wave    *wave)
{
	int		ix, iy, iz, gz, i, n, comp;
	int		M1, M2, MZ;
	int		imax[3], imin[3], icoords[3];
	float		*coords, vel[3], z, v;
	RECT_GRID	*rgr = wave->rect_grid;
	Locstate	st;
	int		st_type;

	init_turb_vel();
	
	for(i=0; i<3; i++)
	{
	    imin[i] = 0;
	    imax[i] = rgr->gmax[i];
		
	    if(rect_boundary_type(front->interf,i,0) == OPEN_BOUNDARY)
		imin[i] = -rgr->lbuf[i];
	    if(rect_boundary_type(front->interf,i,1) == OPEN_BOUNDARY)
		imax[i] = rgr->gmax[i] + rgr->ubuf[i];
	}
	    
	M1 = My + 2*Ny;
	M2 = Mz + 2*Nz;
	MZ = 40;
	v = 20.0;

	iz = imin[2];
	for(gz=0; gz<MZ; ++gz)
	{
	    //shift Ra
	    for(n=0; n<Nx/MZ*2; n++)
		cylic_shift_Ra(M1, M2);

	    //get pvel
	    get_pert_vel();

	    z = rgr->GL[2] + (gz + 0.5)*rgr->h[2];
	    if(z < rgr->L[2] || z > rgr->U[2])
		continue;

	    printf("#gz %d %d\n", gz, iz);

	    icoords[2] = iz;
	    for(iy=imin[1]; iy<imax[1]; ++iy)
	    {
		icoords[1] = iy;
		for (ix = imin[0]; ix < imax[0]; ++ix)
		{
		    icoords[0] = ix;
	    	    
		    coords = Rect_coords(icoords,wave);
	    	    comp = Rect_comp(icoords,wave);
	    	    st = Rect_state(icoords,wave);
		    
		    st_type = state_type(st);
		    set_state(st, TGAS_STATE, st);
	    
		    Press(st) = 10.0;
		    if(comp == 2)
		    {
			Dens(st) = 0.64;
			interp_turb_vel(vel, coords);
			
			Vel(st)[0] = v*0.3*vel[0];
			Vel(st)[1] = v*0.3*vel[1];
			Vel(st)[2] = v + v*0.3*vel[2];
		    }
		    else
		    {
			Dens(st) = 0.016;
			Vel(st)[0] = 0.0;
			Vel(st)[1] = 0.0;
			Vel(st)[2] = 0.0;
		    }
		    
		    set_state(st, st_type, st);
		}
	    }

	    iz++;
	}
	
}

LOCAL	int	pert_step;

void	set_pertn(int);
void	set_pertn(int	step)
{
	pert_step = step;
}

void    tmp_print_gas_state(char*,Locstate);

void    tmp_print_gas_state(
        char            *sname,
        Locstate        st)
{       
        verbose_print_state(sname, st); 
}  

void	get_constant_state(Locstate, int, float*, float);

void    get_constant_state_a(
        Locstate        sl,
        int             comp,
        float           *coords,
        float           t)
{
        int              st_type, tf, xk, yk;
        float            amp, r, ang, epsu, epsv, epsw, v, dh;
        static bool      first = YES;

        if(!is_obstacle_state(sl))
        {
            st_type = state_type(sl);
            set_state(sl, TGAS_STATE, sl);

            if(t<0.05)
                v = -20.0*t/0.05;
            else
                v = -20.0;

            Vel(sl)[0] = 0.0;
            Vel(sl)[1] = 0.0;
//            Vel(sl)[2] = v;
	    Vel(sl)[2] = -120.5;

            if(comp == 3)
            {
                Press(sl) = 40.0;
                Dens(sl) = 0.064;
		
		//MP
		Press(sl) = 0.324;
                Dens(sl) = 0.0000846;

            }
            else
            {
                Press(sl) = 40.0;
                Dens(sl) = 0.64;
		
		//MP
		Press(sl) = 4.9;
                Dens(sl) = 0.00048;
            }

            if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
            {
                if(comp == 3)
                {
                    pdens(sl)[0] = 0.0;
                    pdens(sl)[1] = Dens(sl);
                }
                else
                {
                    pdens(sl)[0] = Dens(sl);
                    pdens(sl)[1] = 0.0;
                }
            }

            set_state(sl, st_type, sl);
        }
}

bool	check_change_inflow_st(Locstate,float*);

bool	check_change_inflow_st(
	Locstate	st,
	float		*coords)
{
	static float	pos[] = {-0.36, 0.0, 0.0};

	if(!debugging("chinst"))
	    return NO;
	
	if(coords[0] < -0.276 && !is_obstacle_state(st))
	{
	    get_constant_state(st, 2, pos, 0.0);
	}
	return YES;
}

//FIXME: Put me somewhere
double 	blasius_ff(
	double eta)
{
//	return 0.5;

	// There exists a preexisting table relating eta values to f values
	// for the blasius function. They will be stored in these arrays
	static double *table_eta = NULL;
	static double *table_f = NULL;
	static int table_size;
	
	// Read in tabular blasius function
	// format is:
	// num_values(int)
	// eta(double)	f(double)
	// ...			...
	if(table_eta == NULL)
	{
		FILE *f = fopen("blas.data","rt");
		int i, numi;
		fscanf(f,"%d\n",&numi);
		table_size = numi;

		uni_array(&table_eta, numi, sizeof(double));
		uni_array(&table_f, numi, sizeof(double));
		for(i=0; i < numi; ++i)
		{
//			fscanf(f, "%f    %f\n",&table_eta[i], &table_f[i]);
			char buf[32];
			fscanf(f,"%s",buf);
			table_eta[i] = atof(buf);
			fscanf(f,"%s",buf);
			table_f[i] = atof(buf);
			printf("%f    %f\n",table_eta[i], table_f[i]);
		}
		fclose(f);
	}
	
	// Do bisection method to find which tabular indices
	// the desired point is on.
	int a=0, b=table_size-1;
	int g = (a+b)/2;
	while(b-a>1)
	{
		if(table_eta[g] < eta)
		{
			a=g;
		}
		else if(table_eta[g] > eta)
		{
			b=g;
		}
		else
		{
			return table_eta[g];
		}
		g = (a+b)/2;
		// FIXME: prevent infinite loop in case of error
	}

	double f1 = table_f[a], f2 = table_f[b];
	double e1 = table_eta[a], e2 = table_eta[b];
	
	// Do linear interpolation from the table data to the known point.
	double coeff = (eta-e1)/(e2-e1);
//	printf("coeff=%f\n",coeff);
	double f= coeff*(f1) + (1.-coeff)*f2;
//	printf("blas_f=%f\n",f);
	return f;
}

double 	get_blasius_vel_coefficient_a(
	double *coords,
	double wall)
{
	// Determine if the point is in the boundary layer. If so, return
	// The solution to the blasius equation for f'(x) at the point giving
	// the velocity coefficient at that point. If the point is outside the BL
	// then return 1,0. This function will always return a double between
	// 0.0 and 1.0. 

	// This function determines the scaling and nu value based on
	// the fact that at x=plate_start, the boundary layer has zero thickness,
	// and some other known point (x,thickness) stored in fixed_bl_thickness.

	// It is assumed that the xflow is in X direction, and plate is boundary 
	// layer is orthogonal to Z direction.


	double fixed_bl_thickness[2] = {0.0, 0.075}; // Set bl thickness to y at x
	double U = 180.0; // magnitude of flow velocity 
	//double nu = 17.4/0.000085; // kinematic viscosity = dynamic viscosity/rho
	double n = 5.0; // number of deltas to make the bdry layer
	double plate_start=-7.1;

	// calculate nu based on the fixed thickness and plate start
	double nu = sqr(fixed_bl_thickness[1]/(double)n)*U/(fixed_bl_thickness[0]-plate_start);
//	printf("nu = %f\n",nu);

	// calculate delta based on nu, coords, etc.
	double delta = sqrt(nu*(coords[0]-plate_start)/U);
	// eta = y/delta
	double eta = (wall-coords[2])/delta;
	// Buffer zones exist which may cause eta to be less than 0
	if(eta < 0.0f)eta=0.0f;

//	printf("for coords %f %f %f: eta=%f: delta=%f: thickness=%f\n", coords[0], coords[1], coords[2], eta, delta,n*delta);

	// if the point is outside the boundary layer, coeff=1.0f
	if((int)eta > n)
		return 1.0f;
	
	return blasius_ff(eta);
}

void	get_constant_state(
	Locstate	sl,
	int		comp,
	float		*coords,
	float		t)
{
	int      	st_type, rn, tn, pertn, i;
	float		r, f, ptheta, pr, pz, v, ang, ang1, dh;
	float		r_nozz, tin;

	if(!is_obstacle_state(sl))
	{
	    st_type = state_type(sl);
	    set_state(sl, TGAS_STATE, sl);
	    
	    if(coords[2] > 3.9)
	    {
        /*orig upto Nov12 */
	/*
		Press(sl) = 9.25;
		Dens(sl) = 7.6e-4;
		Vel(sl)[0] = 0.0;
		Vel(sl)[1] = 0.0;
		Vel(sl)[2] = -130.0;
		//pdens(sl)[0] = Dens(sl);
		//pdens(sl)[1] = 0.0;
                for(i = 0; i < Params(sl)->n_comps; i++)
                    pdens(sl)[i] = 0.0;
                pdens(sl)[4] = Dens(sl);
	*/
		Press(sl) = 12.5;
		Dens(sl) = 12.3e-4;
		Vel(sl)[0] = 0.0;
		Vel(sl)[1] = 0.0;
		Vel(sl)[2] = -123.0;
                for(i = 0; i < Params(sl)->n_comps; i++)
                    pdens(sl)[i] = 0.0;
                pdens(sl)[4] = Dens(sl);
	    }
	    else 
	    {
		Press(sl) = 0.4;
//		Dens(sl) = 8.7e-5;
		Dens(sl) = 12.41e-5; //TMP_TK
		Vel(sl)[0] = 189.4;
		Vel(sl)[1] = 0.0;
		Vel(sl)[2] = 0.0;
                for(i = 0; i < Params(sl)->n_comps; i++)
                    pdens(sl)[i] = 0.0;
                pdens(sl)[0] = Dens(sl);
	    }
	    OHS(sl) = 0.0;

	    set_state(sl, st_type, sl);
	}

	return;
}

void	get_constant_state_init(Locstate, int, float*, float);

void	get_constant_state_init(
	Locstate	sl,
	int		comp,
	float		*coords,
	float		t)
{
	int      	st_type, rn, tn, pertn, i;
	float		r, f, ptheta, pr, pz, v, ang, ang1, dh;
	float		r_nozz, tin;

	if(!is_obstacle_state(sl))
	{
	    st_type = state_type(sl);
	    set_state(sl, TGAS_STATE, sl);
	    
	    if(coords[2] > 3.9)
	    {
	/* upto Nov19
		Press(sl) = 9.25;
		Dens(sl) = 7.6e-4;
		Vel(sl)[0] = 0.0;
		Vel(sl)[1] = 0.0;
		Vel(sl)[2] = -130.0;
		//pdens(sl)[0] = Dens(sl);
		//pdens(sl)[1] = 0.0;
                for(i = 0; i < Params(sl)->n_comps; i++)
                        pdens(sl)[i] = 0.0;
                pdens(sl)[4] = Dens(sl);
	*/
		Press(sl) = 12.5;
		Dens(sl) = 12.3e-4;
		Vel(sl)[0] = 0.0;
		Vel(sl)[1] = 0.0;
		Vel(sl)[2] = -123.0;
                for(i = 0; i < Params(sl)->n_comps; i++)
                        pdens(sl)[i] = 0.0;
                pdens(sl)[4] = Dens(sl);
	    }
	    else 
	    {
		Press(sl) = 0.4;
//		Dens(sl) = 8.7e-5;
		Dens(sl) = 12.41e-5; //TMP_TK
		Vel(sl)[0] = 189.4;
		Vel(sl)[1] = 0.0;
		Vel(sl)[2] = 0.0;
                for(i = 0; i < Params(sl)->n_comps; i++)
                        pdens(sl)[i] = 0.0;
                pdens(sl)[0] = Dens(sl);
//		verbose_print_state("#TK init",sl);
	    }
	    OHS(sl) = 0.0;

	    set_state(sl, st_type, sl);
	}

	return;
}

void	get_state_shift(float*, Locstate);

void	get_state_shift(
	float		*coords,
	Locstate	s)
{
	int	st_type;
	float	x, f;

	st_type = state_type(s);
	set_state(s, TGAS_STATE, s);

	//Vel(s)[0] = 500.0*(coords[2] + 1.0)/2.0;
	//Vel(s)[1] = 0.0;
	//Vel(s)[2] = 0.0;

	x = coords[1] - coords[0];
	if(x > 1.5)
	    x -= 2.0;
	if(x < -0.5)
	    x += 2.0;

	if(x < 0.5)
	    f = 100.0*x;
	else
	    f = 100.0*(1.0-x);

	//f = 50.0*sin(PI*(coords[1]-coords[0])); //*sin(PI*(coords[2]-coords[0]));

	Vel(s)[0] = f;
	Vel(s)[1] = f;
	Vel(s)[2] = 0;

	set_state(s, st_type, s);
}

EXPORT  void    tecplot_surface_geom(
	const char	*bname,
	FILE		*file,
	SURFACE		*s)
{
	TRI	*tri;
	POINT	*p;
	int	i,npts,ntri,num_tris;
	Locstate  sl, sr;
	float	dist, lenk, max_cos;

	if (bname != NULL)//direct call
	{
	    if ((file = fopen(bname,"w")) == NULL)
	    {
		screen("WARNING in tecplot_surface(), "
		       "can't open %s\n",bname);
		return;
	    }
	    (void) fprintf(file,"TITLE = \"tecplot surface\"\n"
			   	"VARIABLES = \"x\", \"y\", \"z\", \"curvature\"," 
				"\"n1\", \"n2\", \"n3\"\n");
	}
	
	//called from tecplot_interface
	if (file == NULL)
	{
	    screen("ERROR, in tecplot_surface, file is NULL\n");
	    clean_up(ERROR);
	}
	if (!(first_tri(s)))
	{
	    screen("WARNING, first bond of the curve is NULL\n");
	    return;
	}

	//count number of points(npts) and number of tris(ntri)
	for (tri=first_tri(s),ntri=0; !at_end_of_tri_list(tri,s); tri=tri->next,ntri++)
	{
	    for (i = 0; i < 3; i++)
	    {
		Index_of_point(Point_of_tri(tri)[i]) = -1;
	    }
	}
	for (tri=first_tri(s),npts=0; !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		p = Point_of_tri(tri)[i];
		if (Index_of_point(p) == -1)
		{
		    Index_of_point(p) = ++npts;
		}
	    }
	}
	//end counting
	
	fprint_wave_type(file, "ZONE T=\"", wave_type(s), "\"", s->interface);
    	fprintf(file, " N=%d E=%d\nF=FEPOINT, ET=TRIANGLE\n",npts,ntri);

	for (tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		Index_of_point(Point_of_tri(tri)[i]) = -1;
	    }
	}
	for (tri=first_tri(s),npts=0; !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		p = Point_of_tri(tri)[i];
		if (Index_of_point(p) == -1)
		{
		    Index_of_point(p) = ++npts;
		    
		    if(wave_type(s) != PASSIVE_BOUNDARY)
		        //fprintf(file,"%15.8e %15.8e %15.8e   %15.8e %15.8e\n",Coords(p)[0],
		    	//   Coords(p)[1],Coords(p)[2], pressure(sl), pressure(sr));
		        fprintf(file,"%15.8e %15.8e %15.8e   %15.8e   %15.8e  %15.8e  %15.8e\n",Coords(p)[0],
		    	   Coords(p)[1],Coords(p)[2], p->curvature, p->_nor[0], p->_nor[1], p->_nor[2]);

		    else
		        fprintf(file,"%15.8e %15.8e %15.8e   %15.8e %15.8e %15.8e %15.8e\n",Coords(p)[0],
		    	 Coords(p)[1],Coords(p)[2], 0.0, 0.0, 0.0, 0.0);
		}
	    }
	}
	for (tri=first_tri(s); !at_end_of_tri_list(tri,s); tri=tri->next)
	{
	    for (i = 0; i < 3; i++)
	    {
		fprintf(file,"%d ",Index_of_point(Point_of_tri(tri)[i]));
	    }
	    fprintf(file,"\n");
	}

	if (ntri != s->num_tri)
	{
	    printf("WARNING, num of tri in surface is wrong\n"); 
	}
	if (bname != NULL)
	    fclose(file);
}	/* end tecplot_surface_states */


EXPORT  void    tecplot_interface_geom(
	const char	*bname,
	INTERFACE	*intfc)
{
	SURFACE	**s;
	CURVE	**cc;
	char    bname1[100];
	FILE	*file;
	int	*i;

	sprintf(bname1, "%s-%d.plt", bname, pp_mynode());
	if ((file = fopen(bname1,"w")) == NULL)
	{
	    screen("WARNING in tecplot_interface(), "
	           "can't open %s\n",bname1);
	    return;
	}
	(void) fprintf(file,"TITLE = \"tecplot interface\"\n"
			    "VARIABLES = \"x\", \"y\", \"z\", \"curvature\"," 
			    "\"n1\", \"n2\", \"n3\"\n");

	for (s = intfc->surfaces; s && *s; ++s)
	{
	    //if(wave_type(*s) != CONTACT)
	    //    continue;
	    printf("#show contact here\n");
	    tecplot_surface_geom(NULL,file,*s);
	}
	fclose(file);
}	/* end tecplot_interface */


int	interp_st_in_direction(Locstate,int*,int,int,Stencil*);

int	interp_st_in_direction(
	Locstate	st,
	int		*icoords,
	int		dir,
	int		new_comp,
	Stencil		*sten)
{
	Wave		*wave = sten->wave;
	INTERFACE	*intfc = sten->fr->interf;
	RECT_GRID	*rgr = wave->rect_grid;
	TRI_GRID	*tg = wave_tri_soln(wave)->tri_grid;
	int		icrds[3], imin[3], imax[3], i, ni, pi, nc;
	int		comp;
	size_t		sizest = wave->sizest;		
	Locstate	stl[2], str[2];
	float		*ptl[2], *ptr[2];
	GRID_DIRECTION  nside[]={EAST,NORTH,UPPER}, pside[]={WEST,SOUTH,LOWER};
	static CRXING	**crxs = NULL;

	if(crxs == NULL)
	    uni_array(&crxs, 20, sizeof(CRXING*));

	ft_assign(icrds, icoords, 3*INT);
	//new_comp = sten->newcomp;

	for(i=0; i<3; i++)
	{
	    imin[i] = -rgr->lbuf[i];
	    imax[i] = rgr->gmax[i] + rgr->ubuf[i];
	}

	//printf("pside  %d %d %d dir %d\n", pside[0], pside[1], pside[2], dir);
	
	for(i=0; i<2; i++)
	{
	    icrds[dir] = icoords[dir] + i + 1;

	    //next point is out of the computational domain.
	    if(icrds[dir] >= imax[dir])
		break;
	    nc = crossings_in_direction(crxs,icrds,pside[dir],tg);
	    comp = Rect_comp(icrds, wave);

	    //have the same comp and
	    //previous intfc have no crx on the edge or 
	    //only one crx (the intfc must move after prop) .
	    if(nc <= 1 && comp == new_comp)
	    {
	        str[i] = Rect_state(icrds, wave);
	        ptr[i] = Rect_coords(icrds, wave);
	    }
	    else
		break;
	}
	ni = i;

	//printf("nside%d %d %d\n", nside[0], nside[1], nside[2]);
	
	for(i=0; i<2; i++)
	{
	    icrds[dir] = icoords[dir] - i - 1;
	    if(icrds[dir] < imin[dir])
		break;
	    nc = crossings_in_direction(crxs,icrds,nside[dir],tg);
	    comp = Rect_comp(icrds, wave);

	    if(nc <= 1 && comp == new_comp)
	    {
		stl[i] = Rect_state(icrds, wave);
		ptl[i] = Rect_coords(icrds, wave);
	    }
	    else
		break;
	}
	pi = i;

	if(debugging("interp_npt"))
	{
	    printf("#st dir=%d  nipi %d %d st = %d\n",dir, ni, pi, st);
	    for(i=0; i<ni; i++)
	        verbose_print_state("str", str[i]);
	    for(i=0; i<pi; i++)
	        verbose_print_state("stl", stl[i]);
	}

	if(ni >= 1 && pi >= 1)
	{
	    bi_interpolate_intfc_states(intfc, 0.5, 0.5, ptl[0], stl[0], 
			ptr[0], str[0], st);
	    return 0;
	}
	else  if(ni == 2)
	{
	    bi_interpolate_intfc_states(intfc, 2.0, -1.0, ptr[0], str[0], 
			ptr[1], str[1], st);
	    return 1;
	}
	else  if(pi == 2)
	{
	    bi_interpolate_intfc_states(intfc, 2.0, -1.0, ptl[0], stl[0], 
			ptl[1], stl[1], st);
	    return 1;
	}
	else  if(ni == 1)
	{
	    ft_assign(st, str[0], sizest);
	    return 2;
	}
	else  if(pi == 1)
	{
	    ft_assign(st, stl[0], sizest);
	    return 2;
	}

	return 99;
}

void	get_trail_state(Locstate,int);

void	get_trail_state(
	Locstate	s, 
	int		comp)
{
	int	st_type, i;

	st_type = state_type(s);
	set_type_of_state(s, TGAS_STATE);

	for(i=0; i<3; i++)
	    Vel(s)[i] = 0.0;

	if(comp == 2)
	{
	    Dens(s) = 0.66;
	    Press(s) = 1.0;
	}
	else
	{
	    Dens(s) = 0.0165;
	    Press(s) = 1.0;
	}

	set_state(s, st_type, s);
}

void	get_nearest_state(
	Locstate	s,
	int		*icoords,
	Stencil		*sten)
{
	Wave		*wave = sten->wave;
	RECT_GRID	*rgr = wave->rect_grid;
	int		imin[3], imax[3], i, icrds[3];
	int		new_comp, comp, ix, iy, iz;
	size_t		sizest = wave->sizest;		
	float		*coords, *crds, minlen, len;
	bool		found;

	for(i=0; i<3; i++)
	{
	    imin[i] = -rgr->lbuf[i];
	    imax[i] = rgr->gmax[i] + rgr->ubuf[i];
	}

	new_comp = sten->newcomp;
	coords = Rect_coords(icoords, wave);
	found = NO;
	minlen = HUGE_VAL;

	print_general_vector("nearest st crds = ", coords, 3, "\n");
	
	for(iz=imin[2]; iz < imax[2]; iz++)
	    for(iy=imin[1]; iy < imax[1]; iy++)
		for(ix=imin[0]; ix < imax[0]; ix++)
		{
		    icrds[0] = ix;
		    icrds[1] = iy;
		    icrds[2] = iz;

		    crds = Rect_coords(icrds, wave);
		    comp = Rect_comp(icrds, wave);
		    if(comp != new_comp)
			continue;

		    len = distance_between_positions(crds, coords, 3);
		    if(len < minlen)
		    {
			found = YES;
			minlen = len;
			ft_assign(s, Rect_state(icrds, wave), sizest);
		    }
		}

	if(!found)
	{
	    printf("WARNING get_const_st, no suitable states.\n");
	    get_trail_state(s, new_comp);
	    //clean_up(ERROR);
	}

}

void	get_st_from_exterp_ave(
	Locstate	s,
	int		*icoords,
	Stencil		*sten)
{
	Wave		*wave = sten->wave;
	INTERFACE	*intfc = sten->fr->interf;
	RECT_GRID	*rgr = wave->rect_grid;
	TRI_GRID	*tg = wave_tri_soln(wave)->tri_grid;
	int		icrds[3], imin[3], imax[3], i, j, k, nc, num, nt[3];
	int		new_comp, comp;
	size_t		sizest = wave->sizest;
	Locstate	st[3], st0;
	float		v0[3], v[3][3], coef, f[3], sum;
	float		*pt = Coords(sten->p[0]);
	static CRXING	**crxs = NULL;
	GRID_DIRECTION  nside[]={EAST,NORTH,UPPER}, pside[]={WEST,SOUTH,LOWER};
	bool		found;

	if(crxs == NULL)
	    uni_array(&crxs, 20, sizeof(CRXING*));
	
	for(i=0; i<3; i++)
	{
	    imin[i] = -rgr->lbuf[i];
	    imax[i] = rgr->gmax[i] + rgr->ubuf[i];
	}

	new_comp = sten->newcomp;

	num = 0;
	for(i=0; i<3; i++)
	{
	    for(j=0; j<2; j++)
	    {
		ft_assign(icrds, icoords, 3*INT);
		if(j == 0)
		    icrds[i] = icoords[i] - 1;
		else
		    icrds[i] = icoords[i] + 1;

		if(icrds[i] < imin[i] || icrds[i] >= imax[i])
		    continue;

		if(j == 0)
		    nc = crossings_in_direction(crxs,icrds,nside[i],tg);
		else
		    nc = crossings_in_direction(crxs,icrds,pside[i],tg);
		comp = Rect_comp(icrds, wave);

		if(nc > 1 || comp != new_comp)
		    continue;

		for(k=0; k<3; k++)
		    v[num][k] = 0.0;
		v[num][i] = j == 0 ? -1.0 : 1.0;
		st[num] = Rect_state(icrds, wave);

		num++;
	    }
	}

	if(num == 0)
	{
	    printf("ERROR get_st_from_exterp, "
		   "get_st_from_inter should be called before.\n");
	    clean_up(ERROR);
	}

	//use the velocity as the interpolation vector.
	st0 = Rect_state(icoords, wave);
	for(i=0; i<3; i++)
	    v0[i] = -vel(i, st0);

	k = 0;
	sum = 0.0;
	for(i=0; i<num; i++)
	{
	    coef = Dot3d(v[i], v0);
	    if(coef > 0)
	    {
		f[k] = coef;
		nt[k] = i;
		k++;
		sum += coef;
	    }
	}

	found = YES;
	if(k == 0)
	{
	    for(i=0; i<num; i++)
	    {
		f[i] = 1.0;
		nt[i] = i;
	    }
	    found = NO;
	    k = num;
	    sum = k;
	}

	//normalize the coef
	for(i=0; i<k; i++)
	    f[i] /= sum;

	if(debugging("interp_npt") || 
	   (icoords[0] == 4 && icoords[1] == 15 && icoords[2] == 3))
	{
	    printf("#exterp k %d  num %d found %d \n", k, num, found);
	    print_general_vector("f=", f, k, "\n");
	    print_int_vector("nt=", nt, k, "\n");
	    print_general_vector("v0=", v0, 3, "\n");
	    for(i=0; i<num; i++)
		print_general_vector("v= ", v[i], 3, "\n");

	    for(i=0; i<k; i++)
		verbose_print_state("st[i]", st[nt[i]]);
	}

	if(k == 1)
	    ft_assign(s, st[nt[0]], sizest);
	else if(k == 2)
	    bi_interpolate_intfc_states(intfc, f[0], f[1], pt, st[nt[0]], 
			pt, st[nt[1]], s);
	else if(k == 3)
	    tri_interpolate_intfc_states(intfc, f[0], f[1], f[2], 
		pt, st[nt[0]], pt, st[nt[1]], pt, st[nt[2]], s);
	
}

int	get_st_from_interp(Locstate, int*, Stencil*);

//s must have its own memory allocated.
int	get_st_from_interp(
	Locstate	s,
	int		*icoords,
	Stencil		*sten)
{
	Wave		*wave = sten->wave;
	INTERFACE	*intfc = sten->fr->interf;
	int		i, k, nt[3], nk[3], num;
	size_t		sizest = wave->sizest;		
	float		*pt = Coords(sten->p[0]);
	static Locstate	st[3];
	static bool	first = YES;

	if(first)
	{
	    for(i=0; i<3; i++)
	    {
		alloc_state(intfc, &(st[i]), sizest);
		printf("#stinit %d \n", st[i]);
	    }
	    first = NO;
	}

	//printf("#interp bf\n");

	add_to_debug("enp_int");
	k = 99;
	for(i=0; i<3; i++)
	{
	    nk[i] = interp_st_in_direction(st[i], icoords, i, 
	    	sten->newcomp, sten);

	    if(nk[i] < k)
		k = nk[i];
	}
	
	//printf("#interp af %d\n", k);

	if(k == 99)
	{
	    get_nearest_state(s, icoords, sten);
	    remove_from_debug("enp_int");
	    return -1;
	}

	//extraplation case.
	if(NO && k == 1)
	{
	    get_st_from_exterp_ave(s, icoords, sten);
	    if(NO && pressure(s) > 100.0)
	    {
		print_int_vector("high pressure ", icoords, 3, "\n");
		clean_up(ERROR);
	    }
	    remove_from_debug("enp_int");
	    return 1;
	}

	num = 0;
	for(i=0; i<3; i++)
	{
	    if(nk[i] == k)
	    {
		nt[num] = i;
		num++;
	    }
	}
	
	if(debugging("interp_npt"))
	{
	    printf("#interp k %d  num %d \n", k, num);
	    for(i=0; i<num; i++)
		verbose_print_state("st[i]", st[nt[i]]);
	}

	if(num == 1)
	    ft_assign(s, st[nt[0]], sizest);
	else if(num == 2)
	{
	    bi_interpolate_intfc_states(intfc, 0.5, 0.5, pt, st[nt[0]], 
			pt, st[nt[1]], s);
	}
	else if(num == 3)
	    tri_interpolate_intfc_states(intfc, 0.333, 0.333, 0.334, 
		pt, st[nt[0]], pt, st[nt[1]], pt, st[nt[2]], s);
	
	if(NO && pressure(s) > 5000.0)
	{
	    print_int_vector("high pressure 1 ", icoords, 3, "\n");
	    clean_up(ERROR);
	}

	remove_from_debug("enp_int");
	return 1;
}

void	get_mid_ghost_pressure(
	float		*pml,
	float		*pmr,
	Locstate	sl,
	Locstate	sr,
	float		pjump)
{
	RIEMANN_SOLVER_WAVE_TYPE l_wave,r_wave;
	size_t		sizest;
	float		uml, umr, ml, mr;
	int		i;
	static Locstate left = NULL, right = NULL;

	if (right == NULL)
	{
	    Gas_param *params=(!is_obstacle_state(sl)) ? Params(sl):Params(sr);
	    
	    sizest = params->sizest;
	    (*params->_alloc_state)(&left,max(sizeof(VGas),sizest));
	    (*params->_alloc_state)(&right,max(sizeof(VGas),sizest));
	}
	
	set_state_for_find_mid_state(right,sr);
	set_state_for_find_mid_state(left,sl);

	for(i=0; i<3; i++)
	    Vel(left)[i] =  Vel(right)[i] = 0.0;

	if(find_mid_state(left,right,pjump,pml,pmr,&uml,&umr,&ml,&mr,
	                       &l_wave,&r_wave) != FUNCTION_SUCCEEDED)
	{
	    (void) printf("WARNING in get_mid_ghost_pressure(), find_mid_state() "
	                  "did not converge for left/right RP\n");
	}

}

void	get_st_from_exterp(Locstate,Locstate,Locstate,
		Stencil*,int,HYPER_SURF*,float,float*);

void	get_st_from_exterp(
	Locstate	s,
	Locstate	s1,
	Locstate	s2,
	Stencil		*sten,
	int		comp,
	HYPER_SURF	*hs,
	float		kappa,
	float		*nor)
{
	int	i, st_type;
	size_t  sizest = sten->wave->sizest;
	float	*h = sten->wave->rect_grid->h, hmax;
	float	ncor, pjump, nor0[3], W[3];
	float	p0, p1;
	static Locstate s1g = NULL, s2g = NULL;
	
	if(s1g == NULL)
	{
	    Gas_param *params=(!is_obstacle_state(s2)) ? Params(s2):Params(s2);
	    
	    (*params->_alloc_state)(&s1g,max(sizeof(VGas),sizest));
	    (*params->_alloc_state)(&s2g,max(sizeof(VGas),sizest));
	}

	hmax = max3(h[0], h[1], h[2]);

	st_type = state_type(s1);
	ft_assign(s, s1, sizest);
	set_state(s, TGAS_STATE, s);

	//if use ncor, only the normal on the other side is used.
	ncor = 0.0;
	if(debugging("ncor"))
	    for(i=0; i<3; i++)
		ncor += (Gvel(s2)[1][i] - Gvel(s2)[0][i])*Gn(s2)[i];

	//pjump = 2.0*Kap(s2)*surface_tension(hs);
	pjump = 2.0*kappa*surface_tension(hs);

	if(debugging("revjump"))
	    pjump *= -1.0;

	if(debugging("mgpress"))
	{
	    //p1 is the pressure on s2 side (ghost side)
	    if(comp == 3)
		get_mid_ghost_pressure(&p0, &p1, s1, s2, pjump);
	    else
		get_mid_ghost_pressure(&p1, &p0, s2, s1, pjump);

	    //printf("#mgp %15.8e  %15.8e    %15.8e  %15.8e\n",
	    //    pressure(s1), p0, pressure(s2), p1);
	}

	if(debugging("igfm"))
	{
	    if(comp == 3)
		g_w_speed(NULL, s1, s2, s1g, s2g, W, 
		          pjump,nor,CONTACT,sten->fr);
	    else
	    {
		for(i=0; i<3; i++)
		    nor0[i] = -nor[i];
	        g_w_speed(NULL, s2, s1, s2g, s1g, W, 
		          pjump,nor0,CONTACT,sten->fr);
	    }
	}
	
	if(the_point(sten->p[0]))
	{
	    printf("#pjump = %15.8e comp=%d\n", pjump, comp);
	    verbose_print_state("s1", s1);
	    verbose_print_state("s2", s2);
	    print_general_vector("Gvel(s2)[0]", Gvel(s2)[0], 3, "\n");
	    print_general_vector("Gvel(s2)[1]", Gvel(s2)[1], 3, "\n");
	    print_general_vector("Gn(s2)", Gn(s2), 3, "\n");
	    printf("ncor = %15.8e\n", ncor);
	    if(debugging("igfm"))
	    {
		print_general_vector("nor", nor, 3, "\n");
		verbose_print_state("s1g", s1g);
		verbose_print_state("s2g", s2g);
	    }
	    printf("\n");
	}

if(debugging("flgasexp"))
{
if(debugging("igfm"))
{
	if(comp == 2)
	{
	    //Dens(s) = Gd(s2)[0];
	    //Dens(s) = Dens(s1);
	    Dens(s) = Dens(s1g);
	    Press(s) = pressure(s1g);
	    for(i=0; i<3; i++)
		Vel(s)[i] = vel(i,s1g);
	}
	else
	{
	    //Dens(s) = Gd(s2)[1];
	    //Dens(s) = Dens(s1);
	    Dens(s) = Dens(s1g);

	    //gas is sensitive to large curvature, for large curvature
	    //use the original gfm alg.
	    //if(fabs(kappa) > 1.0/hmax)
	    if(debugging("gasgfm"))
	    {
		Press(s) = pressure(s1);
		for(i=0; i<3; i++)
		    Vel(s)[i] = Gvel(s2)[1][i] - ncor*Gn(s2)[i];
	    }
	    else
	    {
		Press(s) = pressure(s1g);
		for(i=0; i<3; i++)
		    Vel(s)[i] = vel(i,s1g);
	    }
	}
}
else
{
	if(comp == 2)
	{
	    if(debugging("mgpress"))
		Press(s) = p1 - pjump;
	    else
	        Press(s) = Gp(s2)[1] - pjump;
	   
	    //printf("#pres  %15.8e  %15.8e  %15.8e\n", p1, Gp(s2)[1], Press(s));

	    //if(Press(s) < 0.0)
	    //Press(s) = Gp(s2)[1]/2.0;

	    Dens(s) = Gd(s2)[0];
	    if(debugging("ncor"))
	    {
		for(i=0; i<3; i++)
		    Vel(s)[i] = Gvel(s2)[0][i] + ncor*Gn(s2)[i];
	    }
	    else
	    {
		for(i=0; i<3; i++)
		    Vel(s)[i] = Gvel(s2)[0][i];
	    }
	}
	else
	{
	    if(debugging("mgpress"))
		Press(s) = p1 + pjump;
	    else
		Press(s) = Gp(s2)[1];

	    Dens(s) = Gd(s2)[1];
	    if(debugging("ncor"))
	    {
		for(i=0; i<3; i++)
		    Vel(s)[i] = Gvel(s2)[1][i] - ncor*Gn(s2)[i];
	    }
	    else
	    {
		for(i=0; i<3; i++)
		    Vel(s)[i] = Gvel(s2)[0][i];
	    }
	}
}
}
else
{
if(debugging("igfm"))
{
	if(comp == 2)
	{
	    //Dens(s) = Gd(s2)[0];
	    Dens(s) = Dens(s1);
	    Press(s) = pressure(s1g);
	    for(i=0; i<3; i++)
		Vel(s)[i] = vel(i,s1g);
	}
	else
	{
	    //Dens(s) = Gd(s2)[1];
	    Dens(s) = Dens(s1);

	    //gas is sensitive to large curvature, for large curvature
	    //use the original gfm alg.
	    //if(fabs(kappa) > 1.0/hmax)
	    if(debugging("gasgfm"))
	    {
		Press(s) = pressure(s1);
		for(i=0; i<3; i++)
		    Vel(s)[i] = Gvel(s2)[1][i] - ncor*Gn(s2)[i];
	    }
	    else
	    {
		Press(s) = pressure(s1g);
		for(i=0; i<3; i++)
		    Vel(s)[i] = vel(i,s1g);
	    }
	}
}
else
{
	if(comp == 2)
	{
	    if(debugging("mgpress"))
		Press(s) = p1 - pjump;
	    else
		Press(s) = Gp(s2)[1] - pjump;

	    if(Press(s) < 0.0)
		Press(s) = Gp(s2)[1]/2.0;

	    Dens(s) = Gd(s2)[0];
	    if(debugging("ncor"))
	    {
		for(i=0; i<3; i++)
		    Vel(s)[i] = Gvel(s2)[0][i] + ncor*Gn(s2)[i];
	    }
	    else
	    {
		for(i=0; i<3; i++)
		    Vel(s)[i] = Gvel(s2)[1][i];
	    }
	}
	else
	{
	    if(debugging("mgpress"))
		Press(s) = p1 + pjump;
	    else
		Press(s) = Gp(s2)[0] + pjump;

	    if(Press(s) < 0.0)
		Press(s) = Gp(s2)[0]/2.0;

	    Dens(s) = Gd(s2)[1];
	    if(debugging("ncor"))
	    {
		for(i=0; i<3; i++)
		    Vel(s)[i] = Gvel(s2)[1][i] - ncor*Gn(s2)[i];
	    }
	    else
	    {
		for(i=0; i<3; i++)
		    Vel(s)[i] = Gvel(s2)[0][i];
	    }
	}
}	//if igfm
}

	if(debugging("interp_npt"))
	{
	    printf("#Gd %15.8e\n", Gd(s2)[1]);
	    verbose_print_state("#resu st", s);
	}

	set_state(s, st_type, s);
}

void	get_constant_state(Locstate, int, float*, float);

void interp_dirichlet_states(
	int	 start,
	int	 side,
	int	 endpt,
	int	 new_comp,
	Stencil	 *sten,
	CRXING	 *cross)
{
	int		i, j, indx, in_one, newcomp;	
	size_t		sizest = sten->fr->sizest;
	Locstate	s1, s2;
	int		inone, st_type;

	indx = start*side;
	inone = (start-1)*side;
	
	for (j = 0; j < 3; ++j)
		Coords(sten->p[indx])[j] = Coords(cross->pt)[j];

	s1 = sten->st[inone];
	s2 = sten->st[indx] = sten->worksp[indx];
	ft_assign(s2, s1, sizest);
	
	//the inflow dirichlet bdry velocity
	if(new_comp == 2)
	{
	    get_constant_state(s2, 2, Coords(sten->p[indx]), sten->fr->time);
	    st_type = state_type(s2);
	    set_state(s2, TGAS_STATE, s2);
	    //no this for the run
	    //Dens(s2) = density(s1);
	    //Press(s2) = pressure(s1);
	    set_state(s2, st_type, s2);
	}
	
	//crossflow bdry
	if(new_comp == 3)
	{
	    st_type = state_type(s2);
	    set_state(s2, TGAS_STATE, s2);
	    Dens(s2) = density(s1);
	    Press(s2) = pressure(s1);
	    set_state(s2, st_type, s2);
	}

	if(NO && new_comp == 2)
	{
	    st_type = state_type(s2);
	    set_state(s2, TGAS_STATE, s2);
	    Dens(s2) = density(s1);
	    Press(s2) = 40.0;
	    set_state(s2, st_type, s2);
	}

	ft_assign(left_state(sten->p[indx]),sten->st[indx],sizest);
	ft_assign(right_state(sten->p[indx]),sten->st[indx],sizest);

	for (i = start+1; i <= endpt; ++i)
	{
	    indx = i*side;
	    in_one = indx - side;
	    sten->st[indx] = sten->st[in_one];
	    for (j = 0; j < 3; ++j)
	    	Coords(sten->p[indx])[j] = Coords(sten->p[in_one])[j];
	    ft_assign(left_state(sten->p[indx]),sten->st[indx],sizest);
	    ft_assign(right_state(sten->p[indx]),sten->st[indx],sizest);
	}
	
	return;
}



/*
void	get_st_from_exterp(
	Locstate	s,
	Locstate	s1,
	Locstate	s2,
	Stencil		*sten,
	int		comp)
{
	int	i, st_type;
	size_t  sizest = sten->wave->sizest;

	//verbose_print_state("s1", s1);
	//verbose_print_state("s2", s2);
	
	st_type = state_type(s1);
	ft_assign(s, s1, sizest);
	set_state(s, TGAS_STATE, s);

	if(comp == 2)
	    Press(s) = pressure(s2);     //s1 is fluid
	else
	    for(i=0; i<3; i++)
		Vel(s)[i] = vel(i, s2);  //s1 is gas

	set_state(s, st_type, s);
}
*/

void    print_states_gas(Locstate);
                
void    print_states_gas(Locstate st)
{       
        verbose_print_state("st", st);
}

void	RI_get_minmax_radius(Front *);

void	RI_get_minmax_radius(Front *fr)
{
	INTERFACE		*intfc = fr->interf;
	HYPER_SURF              *hs;
        HYPER_SURF_ELEMENT      *hse;
        POINT                   *p;
	float			dist, mindist, maxdist;

	DEBUG_ENTER(RI_get_minmax_radius)

	mindist = HUGE_VAL;
	maxdist = -HUGE_VAL;
	
	(void) next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
	{
	    if(Boundary_point(p))
		continue;

	    dist = sqrt(sqr(Coords(p)[0]) + sqr(Coords(p)[1]));
	    if(dist < mindist)
		mindist = dist;
	    if(dist > maxdist)
		maxdist = dist;
	}

	pp_global_min(&mindist, 1);
	pp_global_max(&maxdist, 1);

	printf("#RIminmax  %5d  %15.8e  %15.8e  %15.8e\n", 
		fr->step, fr->time, mindist, maxdist);

	DEBUG_LEAVE(RI_get_minmax_radius)
}


EXPORT  void  tecplot_interior_normal(
	char *bname,
	Wave *wave)
{
	const RECT_GRID 	*rgrid = wave->rect_grid;
	Locstate		state;
	int		        xmax = rgrid->gmax[0];
	int			ymax = rgrid->gmax[1];
	int 			zmax = rgrid->gmax[2];
	int 			ix, iy, iz;
	int			icoords[MAXD];
	float			*coords;
	int 			comp;
	FILE			*fp;
	char			s[200];

	sprintf(s,"%s-%d.plt", bname,pp_mynode());
	printf("tecplot_interior_states  file name %s \n",s);

	fp = fopen(s, "w");
	if(fp == NULL)
	{
	    printf("WARNING tecplot_interior_states, can not open file %s\n", s);
	    return; 
	}

	fprintf(fp, "TITLE = \"inner states\" ");
	fprintf(fp, "VARIABLES = \"x\", \"y\", \"z\", \"comp\",  ");
	fprintf(fp, "\"nx\", \"ny\", \"nz\" \n " );
	fprintf(fp, "ZONE i=%d, j=%d, k=%d\n", xmax, ymax, zmax);

	for(iz=0; iz < zmax; iz++)
	    for(iy=0; iy < ymax; iy++)
		  for(ix=0; ix < xmax; ix++)
		  {
		  	icoords[0] = ix;
		  	icoords[1] = iy;
		  	icoords[2] = iz;
			coords = Rect_coords(icoords, wave);
			comp = Rect_comp(icoords, wave);
			state = Rect_state(icoords, wave);

			fprintf(fp, "%f ", coords[0]);
			fprintf(fp, "%f ", coords[1]);
			fprintf(fp, "%f ", coords[2]);
			fprintf(fp, "%d ", comp);

			if(is_obstacle_state(state))
			{
			    fprintf(fp, "%15.8e %15.8e %15.8e ", 0.0, 0.0, 0.0);
			}
			else
			{
			    fprintf(fp, "%15.8e %15.8e %15.8e ", Vel(state)[0],
			            Vel(state)[1], Vel(state)[2]);
			}

			fprintf(fp, "\n");
		  }
	
	fclose(fp);
}

EXPORT  void  tecplot_interior_ghost(char *, Wave*);

EXPORT  void  tecplot_interior_ghost(
	char *bname,
	Wave *wave)
{
	const RECT_GRID 	*rgrid = wave->rect_grid;
	Locstate		state;
	int		        xmax = rgrid->gmax[0];
	int			ymax = rgrid->gmax[1];
	int 			zmax = rgrid->gmax[2];
	int 			ix, iy, iz, i;
	int			icoords[MAXD];
	float			*coords;
	int 			comp, imin[3], imax[3];
	FILE			*fp;
	char			s[200];

	sprintf(s,"%s-%d.plt", bname,pp_mynode());
	printf("tecplot_interior_states  file name %s \n",s);

	fp = fopen(s, "w");
	if(fp == NULL)
	{
	    printf("WARNING tecplot_interior_states, can not open file %s\n", s);
	    return; 
	}
	
	for(i=0; i<3; i++)
	{
	    //imin[i] = -rgrid->lbuf[i];
	    //imax[i] = rgrid->gmax[i] + rgrid->ubuf[i];
	    imin[i] = 0;
	    imax[i] = rgrid->gmax[i];
	}

	fprintf(fp, "TITLE = \"inner states\" ");
	fprintf(fp, "VARIABLES = \"x\", \"y\", \"z\", \"comp\",  ");
	fprintf(fp, "\"dens\", \"press\", \"u\", \"v\", \"w\"  " );
	fprintf(fp, "\"dens1\", \"press1\", \"u1\", \"v1\", \"w1\" " );
	fprintf(fp, "\"nx\", \"ny\", \"nz\", \"kappa\"  \n " );
	fprintf(fp, "ZONE i=%d, j=%d, k=%d\n", xmax, ymax, zmax);
	
	//for(iz=0; iz < zmax; iz++)
	//    for(iy=0; iy < ymax; iy++)
	//	  for(ix=0; ix < xmax; ix++)
	for(iz=imin[2]; iz < imax[2]; iz++)
	    for(iy=imin[1]; iy < imax[1]; iy++)
		  for(ix=imin[0]; ix < imax[0]; ix++)
		  {
		  	icoords[0] = ix;
		  	icoords[1] = iy;
		  	icoords[2] = iz;
			coords = Rect_coords(icoords, wave);
			comp = Rect_comp(icoords, wave);
			state = Rect_state(icoords, wave);

			fprintf(fp, "%f ", coords[0]);
			fprintf(fp, "%f ", coords[1]);
			fprintf(fp, "%f ", coords[2]);
			fprintf(fp, "%d ", comp);
		
			for(i=0; i<2; i++)
			    if(is_obstacle_state(state))
			    {
				fprintf(fp, "%15.8e %15.8e %15.8e %15.8e %15.8e  ", 
					0.0, 0.0, 0.0, 0.0, 0.0);
			    }
			    else
			    {
				fprintf(fp, "%15.8e %15.8e %15.8e %15.8e %15.8e  ", Gd(state)[i], Gp(state)[i],
			            Gvel(state)[i][0], Gvel(state)[i][1], Gvel(state)[i][2]);
			    }

			if(is_obstacle_state(state))
			    fprintf(fp, "%15.8e %15.8e %15.8e %15.8e ", 
					0.0, 0.0, 0.0, 0.0);
			else
			    fprintf(fp, "%15.8e %15.8e %15.8e %15.8e ", 
					Gn(state)[0], Gn(state)[1], Gn(state)[2], Kap(state));

			fprintf(fp, "\n");
		  }
	
	fclose(fp);
}



void	check_interior_val(Wave*);

void	get_normal_from_front(Front*, Wave*);
void	solve_exp_value(Front*, Wave*);
void	ghost_rk_forward(Front*,Wave*,int,int,int);

void	solve_exp_value(
	Front		*fr,
	Wave		*wave)
{
	int		i, j, k, len, num, imin[3], imax[3], icrds[3];
	RECT_GRID	*rgr = wave->rect_grid;
	byte		*state_storage, *storage;
	Locstate	*store_state, st, st1;
	size_t		sizest = wave->sizest;
	int		*iperm, iter, iii;

	len = 1;
	for(i=0; i<3; i++)
	{
	    imin[i] = -rgr->lbuf[i];
	    imax[i] = rgr->gmax[i] + rgr->ubuf[i];
	    len *= imax[i] - imin[i];
	}

	uni_array(&state_storage, len, sizest);
	uni_array(&store_state, len, sizeof(Locstate));
	storage = state_storage;

	for(i=0; i<len; i++)
	{
	    store_state[i] = (Locstate) storage;
	    storage += sizest;
	}

	printf("#normal bf store states\n");
	
	//store states.
	num = 0;
	for(k=imin[2]; k<imax[2]; k++)
	    for(j=imin[1]; j<imax[1]; j++)
		for(i=imin[0]; i<imax[0]; i++)
		{
		    icrds[0] = i;
		    icrds[1] = j;
		    icrds[2] = k;
		    ft_assign(store_state[num], Rect_state(icrds, wave), sizest);
		    num++;
		}

	printf("#normal bf\n");
	fflush(NULL);

	fr->interf->normal_unset = YES;
	fr->interf->curvature_unset = YES;

	//Vel(st) will be replaced with normal
	get_normal_from_front(fr, wave);

	//check_interior_val(wave);

	//scatter the normals	
	add_to_debug("gfmref");
	
	iperm = set_iperm(fr->step,3);
	for(i=0; i<3; i++)
	{
	    if(!scatter_states(wave,fr,iperm,i))
	    {
		screen("ERROR ghost: scatter_states() failed in direction\n",i);
		clean_up(ERROR);
	    }
	}

	remove_from_debug("gfmref");

	//check_interior_val(wave);
	//tecplot_interior_ghost("ghost0", wave);
	//tecplot_interior_normal("normal", wave);

	//extrapolate ghost states
	//Dens(st) will be used as a temp storage

if(!debugging("gfmnoit"))	
{
	for(iter=0; iter<5; iter++)
	{
	    //ghost for fluid
	    ghost_rk_forward(fr, wave, 2, 0, 0);    //dens
	    if(debugging("flgasexp"))
	    {
		ghost_rk_forward(fr, wave, 2, 0, 2);    //u fluid
		ghost_rk_forward(fr, wave, 2, 0, 3);    //v fluid
		ghost_rk_forward(fr, wave, 2, 0, 4);    //w fluid
	    }

	    //ghost for gas
	    ghost_rk_forward(fr, wave, 3, 1, 0);    //dens
	    if(debugging("flgasexp"))
	    {
		ghost_rk_forward(fr, wave, 3, 1, 1);    //pressure
	    }

	    printf("#ghost iter %d\n", iter);

	    for(i=0; i<3; i++)
	    {
		if(!scatter_states(wave,fr,iperm,i))
		{
		    screen("ERROR ghost iter: "
			   "scatter_states() failed in direction\n",i);
		    clean_up(ERROR);
		}
	    }
	}

	printf("#ghost iter af\n");
}

	//restore Dens, Vel fields for states.
	num = 0;
	for(k=imin[2]; k<imax[2]; k++)
	    for(j=imin[1]; j<imax[1]; j++)
		for(i=imin[0]; i<imax[0]; i++)
		{
		    icrds[0] = i;
		    icrds[1] = j;
		    icrds[2] = k;
		    //ft_assign(Rect_state(icrds, wave), store_state[num], sizest);

		    st = Rect_state(icrds, wave);
		    st1 = store_state[num];

		    Dens(st) = Dens(st1);
		    Press(st) = Press(st1);
                    if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                    {
                        for(iii = 0; iii < Params(st)->n_comps; iii++)
                            pdens(st)[iii] = pdens(st1)[iii];
                    }

		    ft_assign(Gn(st), Vel(st), 3*FLOAT);
		    ft_assign(Vel(st), Vel(st1), 3*FLOAT);

		    num++;
		}
	
	if(fr->step%10 == -1)
	{
	    char fname[200];

	    //sprintf(fname, "ghost_%d", fr->step);
	    //tecplot_interior_ghost(fname, wave);
	    tecplot_interface_geom("geom", fr->interf);
	    clean_up(0);
	}

	free_these(2, state_storage, store_state);

}

void	ghost_rk_forward(Front*,Wave*,int,int,int);

void	ghost_rk_forward(
	Front		*fr,
	Wave		*wave,
	int		comp,
	int		ind,
	int		flag)
{
	int		i, j, k, dir, imin[3], imax[3], icrds[3], ic[3];
	Locstate	st, st1;
	RECT_GRID	*rgr = wave->rect_grid;
	float		*h = rgr->h;
	float		tmp, an, dt;

	for(i=0; i<3; i++)
	{
	    imin[i] = -rgr->lbuf[i];
	    imax[i] = rgr->gmax[i] + rgr->ubuf[i];
	}
	dt = min3(h[0], h[1], h[2])/3.0;

	//half step forward and store the results in Dens(st)

	//init all the states as its initial value
	for(k=imin[2]; k<imax[2]; k++)
	    for(j=imin[1]; j<imax[1]; j++)
		for(i=imin[0]; i<imax[0]; i++)
		{
		    icrds[0] = i;
		    icrds[1] = j;
		    icrds[2] = k;
			
		    st = Rect_state(icrds, wave);

		    if(flag == 0)
			Dens(st) = Gd(st)[ind];
		    else if(flag == 1)
			Dens(st) = Gp(st)[ind];
		    else if(flag >= 2 && flag <=4)
			Dens(st) = Gvel(st)[ind][flag-2];
			
		    if(isnan(Dens(st)))
		    {
		        printf("#dens bf %15.8e\n", Dens(st));
			print_general_vector("nor=", Vel(st), 3, "\n");
			print_int_vector("icrds=", icrds, 3, "\n");
			clean_up(ERROR);
		    }
		}

	//calc half step values from Gd Gp or Gvel to Dens(st)
	for(k=imin[2]+1; k<imax[2]-1; k++)
	    for(j=imin[1]+1; j<imax[1]-1; j++)
		for(i=imin[0]+1; i<imax[0]-1; i++)
		{
		    icrds[0] = i;
		    icrds[1] = j;
		    icrds[2] = k;

		    st = Rect_state(icrds, wave);

		    if(Rect_comp(icrds, wave) == comp)
			continue;

		    //make sure the normal is calculated at the point.
		    if(Mag3d(Vel(st)) < 0.1)
			continue;

		    tmp = 0.0;
		    for(dir=0; dir<3; dir++)
		    {
			ft_assign(ic, icrds, 3*INT);
			an = Vel(st)[dir]*dt/h[dir];

			if((Vel(st)[dir] < 0.0 && comp == 2) ||
			   (Vel(st)[dir] > 0.0 && comp == 3))
			{
			    ic[dir]++;
			    st1 = Rect_state(ic, wave);

			    //if no normal ic calculated, no flux
			    if(Mag3d(Vel(st1)) < 0.1)
				st1 = st;

			    if(flag == 0)
				tmp += an*(Gd(st1)[ind]-Gd(st)[ind]);
			    else if(flag == 1)
				tmp += an*(Gp(st1)[ind]-Gp(st)[ind]);
			    else if(flag >= 2 && flag <=4)
				tmp += an*(Gvel(st1)[ind][flag-2]-Gvel(st)[ind][flag-2]);
			}
			else
			{
			    ic[dir]--;
			    st1 = Rect_state(ic, wave);
			    if(Mag3d(Vel(st1)) < 0.1)
				st1 = st;

			    if(flag == 0)
				tmp += an*(Gd(st)[ind]-Gd(st1)[ind]);
			    else if(flag == 1)
				tmp += an*(Gp(st)[ind]-Gp(st1)[ind]);
			    else if(flag >= 2 && flag <=4)
				tmp += an*(Gvel(st)[ind][flag-2]-Gvel(st1)[ind][flag-2]);
			}
			if(isnan(tmp))
			{
			    printf("Gd nan  %d %15.8e %15.8e %15.8e\n", dir, an, Gd(st1)[ind], Gd(st)[ind]);
			    printf("%15.8e  %15.8e  %15.8e \n", Vel(st)[dir], dt, h[dir]);
			}
		    }
		    if(comp == 2)
			tmp *= -1.0;

		    if(flag == 0)
			Dens(st) = 0.5*tmp + Gd(st)[ind];
		    else if(flag == 1)
			Dens(st) = 0.5*tmp + Gp(st)[ind];
		    else if(flag >= 2 && flag <=4)
			Dens(st) = 0.5*tmp + Gvel(st)[ind][flag-2];
		}

	//update ghost values from Dens(st)
	for(k=imin[2]+1; k<imax[2]-1; k++)
	    for(j=imin[1]+1; j<imax[1]-1; j++)
		for(i=imin[0]+1; i<imax[0]-1; i++)
		{
		    icrds[0] = i;
		    icrds[1] = j;
		    icrds[2] = k;
		    
		    st = Rect_state(icrds, wave);
		    if(Rect_comp(icrds, wave) == comp)
			continue;

		    //make sure the normal is calculated at the point.
		    if(Mag3d(Vel(st)) < 0.1)
			continue;

		    tmp = 0.0;
		    for(dir=0; dir<3; dir++)
		    {
			ft_assign(ic, icrds, 3*INT);
			an = Vel(st)[dir]*dt/h[dir];

			if((Vel(st)[dir] < 0.0 && comp == 2) ||
			   (Vel(st)[dir] > 0.0 && comp == 3))

			{
			    ic[dir]++;
			    st1 = Rect_state(ic, wave);

			    //if no normal ic calculated, no flux
			    if(Mag3d(Vel(st1)) < 0.1)
				st1 = st;
			    tmp += an*(Dens(st1)-Dens(st));
			}
			else
			{
			    ic[dir]--;
			    st1 = Rect_state(ic, wave);
			    if(Mag3d(Vel(st1)) < 0.1)
				st1 = st;
			    tmp += an*(Dens(st)-Dens(st1));
			}
		    }
		    if(comp == 2)
			tmp *= -1.0;

		    {
			float	*crds = Rect_coords(icrds, wave);
			if(the_pt(crds) && NO)
			{
			    printf("#iter af %d %d %15.8e \n", comp, ind, tmp);
			    print_general_vector("nor=", Vel(st), 3, "\n");
			    printf("#vel %15.8e  %15.8e  %15.8e \n", Gd(st)[ind], 
				Dens(st), Gd(st)[ind]+tmp);
			}
		    }
		    if(flag == 0)
			Gd(st)[ind] += tmp;
		    else if(flag == 1)
			Gp(st)[ind] += tmp;
		    else if(flag >= 2 && flag <=4)
			Gvel(st)[ind][flag-2] += tmp;
		}

}


void	get_ghost_state(Front*, Wave*, int, int);

void	solve_ghost_value(Front*, Wave*);

void	solve_ghost_value(
	Front		*fr,
	Wave		*wave)
{

	get_ghost_state(fr, wave, 2, 0);
	get_ghost_state(fr, wave, 3, 1);
	//tecplot_interior_ghost("ghost", wave);
	//clean_up(0);
}


bool	get_ave_normal(
	int		*imin,
	int		*imax,
	int		*ic,
	Wave		*wave,
	int		***norset)
{
	float		*coords, *coords1, len, f;
	int		i, j, k, ic1[3], dir, num;
	Locstate	st, st1;
	bool		found;

	found = NO;

	for(i=0; i<3; i++)
	    for(j=0; j<2; j++)
	    {
		dir = j == 0 ? -1 : 1;
		ft_assign(ic1, ic, 3*INT);
		ic1[i] = ic[i] + dir;

		if(ic1[i] < imin[i] || ic1[i] >= imax[i])
		    continue;
		if(norset[ic1[0]-imin[0]][ic1[1]-imin[1]][ic1[2]-imin[2]] == 1)
		    found = YES;
	    }
	if(!found)
	    return NO;

	st = Rect_state(ic, wave);
	Vel(st)[0] = 0.0;
	Vel(st)[1] = 0.0;
	Vel(st)[2] = 0.0;
	Kap(st) = 0.0;

	num = 0;
	coords = Rect_coords(ic, wave);
	for(i=ic[0]-1; i<=ic[0]+1; i++)
	    for(j=ic[1]-1; j<=ic[1]+1; j++)
		for(k=ic[2]-1; k<=ic[2]+1; k++)
		{
		    if(i < imin[0] || i >= imax[0] || 
		       j < imin[1] || j >= imax[1] ||
		       k < imin[2] || k >= imax[2])
			continue;
		    if(norset[i-imin[0]][j-imin[1]][k-imin[2]] != 1)
			continue;

		    ic1[0] = i;
		    ic1[1] = j;
		    ic1[2] = k;
		    st1 = Rect_state(ic1, wave);
		    
		    //coords1 = Rect_coords(ic1, wave);
		    //len = distance_between_positions(coords, coords1, 3);

		    //do not use length weighted normal direction
		    Vel(st)[0] += Vel(st1)[0];
		    Vel(st)[1] += Vel(st1)[1];
		    Vel(st)[2] += Vel(st1)[2];
		    Kap(st) += Kap(st1);
		    num++;
		}

	f = Mag3d(Vel(st));

	if(f < 1.0e-6)
	{
	    //print_int_vector("ica", ic, 3, "\n");
	    //print_general_vector("n", Vel(st), 3, "\n");

	    Vel(st)[0] = 0.0;
	    Vel(st)[1] = 0.0;
	    Vel(st)[2] = 0.0;
	    
	}
	else
	{
	    Vel(st)[0] /= f;
	    Vel(st)[1] /= f;
	    Vel(st)[2] /= f;
	}

	if(num > 0)
	    Kap(st) /= num;

	return YES;
}

void	check_interior_val(Wave*);
void	check_interior_val(
	Wave		*wave)
{
	RECT_GRID		*rgr = wave->rect_grid;
	Locstate		st;
	int			i,j,k, imin[3],imax[3], num;
	int			ic[3],ic1[3],ic2[3];

	for(i=0; i<3; i++)
	{
	    imin[i] = -rgr->lbuf[i];
	    imax[i] = rgr->gmax[i] + rgr->ubuf[i];
	}

	for(i=imin[0]; i<imax[0]; i++)
	    for(j=imin[1]; j<imax[1]; j++)
		for(k=imin[2]; k<imax[2]; k++)
		{
		    ic[0] = i;
		    ic[1] = j;
		    ic[2] = k;
		    st = Rect_state(ic, wave);

		    if(Mag3d(Vel(st)) < 1.0e-6)
		    {
		        printf("#check small found.\n");
			print_int_vector("ic=", ic, 3, "\n");
		    }
		}
}

void	get_normal_from_front(Front*, Wave*);

//it will fill Vel field of states by normal
void	get_normal_from_front(
	Front		*fr,
	Wave		*wave)
{
	INTERFACE               *intfc = fr->interf;
	RECT_GRID		*gr = &topological_grid(intfc);
	RECT_GRID		*rgr = wave->rect_grid;
	HYPER_SURF              *hs;
	HYPER_SURF_ELEMENT      *hse;
	POINT                   *p;
	Locstate		st;
	int			i,j,k, imin[3],imax[3], num;
	int			ic[3],ic1[3],ic2[3];
	float			curv, nor[3], d[3], pt1[3], f, d1, d2, d3, *pt, tol;
	int			***norset;
	bool			found;

	for(i=0; i<3; i++)
	{
	    imin[i] = -rgr->lbuf[i];
	    imax[i] = rgr->gmax[i] + rgr->ubuf[i];
	}

	tri_array(&norset,imax[0]-imin[0],imax[1]-imin[1],imax[2]-imin[2],INT);

	tol = min3(rgr->h[0], rgr->h[1], rgr->h[2])*1.0e-6;

	for(i=imin[0]; i<imax[0]; i++)
	    for(j=imin[1]; j<imax[1]; j++)
		for(k=imin[2]; k<imax[2]; k++)
		{
		    ic[0] = i;
		    ic[1] = j;
		    ic[2] = k;
		    st = Rect_state(ic, wave);
		    Vel(st)[0] = 0.0;
		    Vel(st)[1] = 0.0;
		    Vel(st)[2] = 0.0;
		    Kap(st) = 0.0;
		    Gn(st)[0] = 0.0;   //sum of weight
		    Gn(st)[1] = -1.0;  //flag
		    norset[i-imin[0]][j-imin[1]][k-imin[2]] = 0;
		}

	printf("#normal cal bf\n");

	(void) next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
	{
	    if(wave_type(hs) != FIRST_SCALAR_PHYSICS_WAVE_TYPE)
		continue;
	    if(Boundary_point(p))
	    {
		p->_nor[0] = 0.0;
		p->_nor[1] = 0.0;
		p->_nor[2] = 0.0;
		continue;
	    }
	    normal(p, hse, hs, nor, fr);
	    curv = p->curvature;

	    //printf("#curv %15.8e %15.8e \n", p->curvature, 
	    //	f_mean_curvature_at_point3d(p,hse,hs,fr));

	    //|nor| == 1
	    pt = Coords(p);
	   
	    //make sure pt is inside the topological grid.
	    for(i=0; i<3; i++)
		if(pt[i] <= gr->VL[i]+tol || pt[i] >= gr->VU[i]-tol)
		    break;
	    if(i < 3)
		continue;

	    //in which cell of top grid
	    if(!rect_in_which(pt, ic, gr))
		continue;

	    //corner
	    for(i=0; i<3; i++)
		pt1[i] = cell_edge(ic[i], i, gr);
	    //in which cell of comp grid
	    rect_in_which(pt1, ic1, rgr);

	    for(i=0; i<3; i++)
	    {
		d[i] = fabs(pt[i]-pt1[i])/rgr->h[i];
	        if(d[i] < -tol || d[i] > 1.0 + tol)
		{
		    printf("ERROR get_normal_from_front "
		           "d[%d] = %15.8e is out of range.\n",i,d[i]);
		    clean_up(ERROR);
		}
	    }

	    for(i=0; i<2; i++)
		for(j=0; j<2; j++)
		    for(k=0; k<2; k++)
		    {
			ic2[0] = ic1[0] + i;
			ic2[1] = ic1[1] + j;
			ic2[2] = ic1[2] + k;
			d1 = i == 0 ? fabs(1.0-d[0]) : d[0];
			d2 = j == 0 ? fabs(1.0-d[1]) : d[1];
			d3 = k == 0 ? fabs(1.0-d[2]) : d[2];
			st = Rect_state(ic2, wave);
			f = d1*d2*d3;

			if(ic2[0] == 5 && ic2[1] == 13 && ic2[2] == 2)
			{
			    print_general_vector("nor=", nor, 3, "\n");
			    print_general_vector("d=", d, 3, "\n");
			    printf("f= %15.8e\n", f);
			}
			Vel(st)[0] += nor[0]*f;
			Vel(st)[1] += nor[1]*f;
			Vel(st)[2] += nor[2]*f;
			Kap(st) += curv*f;
			Gn(st)[0] += f;
			Gn(st)[1] = 1.0;
		    
			norset[ic2[0]-imin[0]][ic2[1]-imin[1]][ic2[2]-imin[2]] = 1;
		    }
	}

	printf("#normal cal af\n");

	for(i=imin[0]; i<imax[0]; i++)
	    for(j=imin[1]; j<imax[1]; j++)
		for(k=imin[2]; k<imax[2]; k++)
		{
		    ic[0] = i;
		    ic[1] = j;
		    ic[2] = k;

		    //make sure Vel(st) is assigned 
		    if(norset[i-imin[0]][j-imin[1]][k-imin[2]] != 1)
			continue;

		    st = Rect_state(ic, wave);
		    f = Mag3d(Vel(st));

		    if(f < 1.0e-6)
		    {
		        //print_int_vector("zic", ic, 3, "\n");
			//printf("f=  %15.8e\n", f);
			//print_general_vector("n", Vel(st), 3, "\n");

			Vel(st)[0] = 0.0;
			Vel(st)[1] = 0.0;
			Vel(st)[2] = 0.0;

		    }
		    else
		    {
			Vel(st)[0] /= f;
			Vel(st)[1] /= f;
			Vel(st)[2] /= f;
		    }

		    if(Gn(st)[1] > 0.0 && Gn(st)[0] > 1.0e-6)
			Kap(st) /= Gn(st)[0];
		}

	found = YES;
	num = 1;
	while(found && num > 0)
	{
	    found = NO;
	    num = 0;

	    for(i=imin[0]; i<imax[0]; i++)
		for(j=imin[1]; j<imax[1]; j++)
		    for(k=imin[2]; k<imax[2]; k++)
		    {
			if(norset[i-imin[0]][j-imin[1]][k-imin[2]] != 0)
			    continue;

			found = YES;
			ic[0] = i;
			ic[1] = j;
			ic[2] = k;

			if(get_ave_normal(imin,imax,ic,wave,norset))
			{
			    num++;
			    norset[i-imin[0]][j-imin[1]][k-imin[2]] = 2;
			}
		    }

	    printf("found = %d  %d\n", found, num);

	    for(i=imin[0]; i<imax[0]; i++)
		for(j=imin[1]; j<imax[1]; j++)
		    for(k=imin[2]; k<imax[2]; k++)
			if(norset[i-imin[0]][j-imin[1]][k-imin[2]] == 2)
			    norset[i-imin[0]][j-imin[1]][k-imin[2]] = 1;

	}

	free(norset);
}

bool	get_ave_state(
	int		*imin,
	int		*imax,
	int		*ic,
	Wave		*wave,
	int		***norset,
	int		comp,
	int		ind)
{
	float		*coords, *coords1, len, f;
	int		i, j, k, num, ic1[3], dir, iii;
	Locstate	st, st1;
	float		gd, gp, gvel[3], gpd[10];
	bool		found;

	found = NO;

	for(i=0; i<3; i++)
	    for(j=0; j<2; j++)
	    {
		dir = j == 0 ? -1 : 1;
		ft_assign(ic1, ic, 3*INT);
		ic1[i] = ic[i] + dir;

		if(ic1[i] < imin[i] || ic1[i] >= imax[i])
		    continue;
		if(norset[ic1[0]-imin[0]][ic1[1]-imin[1]][ic1[2]-imin[2]] == 1)
		    found = YES;
	    }
	if(!found)
	    return NO;

	st = Rect_state(ic, wave);

	num = 0;
	gd = 0.0;
	gp = 0.0;
	gvel[0] = 0.0;
	gvel[1] = 0.0;
	gvel[2] = 0.0;

	coords = Rect_coords(ic, wave);
	for(i=ic[0]-1; i<=ic[0]+1; i++)
	    for(j=ic[1]-1; j<=ic[1]+1; j++)
		for(k=ic[2]-1; k<=ic[2]+1; k++)
		{
		    if(i < imin[0] || i >= imax[0] || 
		       j < imin[1] || j >= imax[1] ||
		       k < imin[2] || k >= imax[2])
			continue;
		    if(norset[i-imin[0]][j-imin[1]][k-imin[2]] != 1)
			continue;

		    ic1[0] = i;
		    ic1[1] = j;
		    ic1[2] = k;
		    st1 = Rect_state(ic1, wave);
		    coords1 = Rect_coords(ic1, wave);

		    if(Rect_comp(ic1, wave) == comp)
		    {
			gd += density(st1);
			gvel[0] += vel(0,st1);
			gvel[1] += vel(1,st1);
			gvel[2] += vel(2,st1);
			gp += pressure(st1);
                        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                        {
                            for(iii=0; iii<Params(st1)->n_comps; iii++)
                                gpd[iii] += pdens(st1)[iii];
                        }
		    }
		    else
		    {
			gd += Gd(st1)[ind];
			gvel[0] += Gvel(st1)[ind][0];
			gvel[1] += Gvel(st1)[ind][1];
			gvel[2] += Gvel(st1)[ind][2];
			gp += Gp(st1)[ind];
                        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                        {
                            for(iii=0; iii<Params(st1)->n_comps; iii++)
                                gpd[iii] += Gpd(st1)[ind][iii];
                        }
		    }
		    num++;
		}

	Gd(st)[ind] = gd/num;
	Gvel(st)[ind][0] = gvel[0]/num;
	Gvel(st)[ind][1] = gvel[1]/num;
	Gvel(st)[ind][2] = gvel[2]/num;
	Gp(st)[ind] = gp/num;
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            for(iii=0; iii<Params(st)->n_comps; iii++)
                Gpd(st)[ind][iii] = gpd[iii]/num;
	    
	    float sum = 0.0;
            for(iii=0; iii<Params(st)->n_comps; iii++)
                sum += Gpd(st)[ind][iii];
            for(iii=0; iii<Params(st)->n_comps; iii++)
                Gpd(st)[ind][iii] = Gd(st)[ind]/sum*Gpd(st)[ind][iii];
        }

	if(NO && fabs(Gvel(st)[ind][2]) > 10.0)
	{
	    printf("#Gvel %d  %15.8e\n", ind, Gvel(st)[ind][2]);
	    print_int_vector("ic= ", ic, 3, "\n");
	    clean_up(ERROR);
	}
	return YES;
}

void	copy_ghost_state(Wave*, Wave*);

void	copy_ghost_state(
	Wave		*newwave,
	Wave		*wave)
{
	RECT_GRID		*rgr = wave->rect_grid;
	int			i, j, k, ind, ic[3], imin[3], imax[3],iii;
	Locstate		st, st1;

	if(wave == newwave)
	    return;

	for(i=0; i<3; i++)
	{
	    imin[i] = -rgr->lbuf[i];
	    imax[i] = rgr->gmax[i] + rgr->ubuf[i];
	}
	
	for(i=imin[0]; i<imax[0]; i++)
	    for(j=imin[1]; j<imax[1]; j++)
		for(k=imin[2]; k<imax[2]; k++)
		{
		    ic[0] = i;
		    ic[1] = j;
		    ic[2] = k;
		    st = Rect_state(ic, newwave);
		    st1 = Rect_state(ic, wave);
		    for(ind=0; ind<2; ind++)
		    {
			Gd(st)[ind] = Gd(st1)[ind];
			Gp(st)[ind] = Gp(st1)[ind];
			Gvel(st)[ind][0] = Gvel(st1)[ind][0];
			Gvel(st)[ind][1] = Gvel(st1)[ind][1];
			Gvel(st)[ind][2] = Gvel(st1)[ind][2];
                        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
                        {
                            for(iii=0; iii<Params(st)->n_comps; iii++)
                                Gpd(st)[ind][iii] = Gpd(st1)[ind][iii];
                        }
		    }
		    Gn(st)[0] = Gn(st1)[0];
		    Gn(st)[1] = Gn(st1)[1];
		    Gn(st)[2] = Gn(st1)[2];
		    Kap(st) = Kap(st1);
		}

}

//init Gd Gvel Gp from averaging.

void	get_ghost_state(
	Front		*fr,
	Wave		*wave,
	int		comp,
	int		ind)
{
	INTERFACE               *intfc = fr->interf;
	RECT_GRID		*gr = &topological_grid(intfc);
	RECT_GRID		*rgr = wave->rect_grid;
	HYPER_SURF              *hs;
	HYPER_SURF_ELEMENT      *hse;
	POINT                   *p;
	Locstate		st;
	int			i,j,k, imin[3],imax[3];
	int			ic[3],ic1[3],ic2[3];
	float			nor[3], d[3], pt1[3], f, d1, d2, d3, *pt;
	int			***norset, c, num, iii;
	bool			found;

	for(i=0; i<3; i++)
	{
	    imin[i] = -rgr->lbuf[i];
	    imax[i] = rgr->gmax[i] + rgr->ubuf[i];
	}

	tri_array(&norset,imax[0]-imin[0],imax[1]-imin[1],imax[2]-imin[2],INT);
	
	for(i=imin[0]; i<imax[0]; i++)
	    for(j=imin[1]; j<imax[1]; j++)
		for(k=imin[2]; k<imax[2]; k++)
		{
		    ic[0] = i;
		    ic[1] = j;
		    ic[2] = k;
		    c = Rect_comp(ic, wave);
		    st = Rect_state(ic, wave);
		    if(c == comp)
		    {
			norset[i-imin[0]][j-imin[1]][k-imin[2]] = 1;
			Gd(st)[ind] = density(st);
			Gp(st)[ind] = pressure(st);
			Gvel(st)[ind][0] = vel(0, st);
			Gvel(st)[ind][1] = vel(1, st);
			Gvel(st)[ind][2] = vel(2, st);
			if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
			{
			    for(iii=0; iii<Params(st)->n_comps; iii++)
				Gpd(st)[ind][iii] = pdens(st)[iii];
			}
		    }
		    else
			norset[i-imin[0]][j-imin[1]][k-imin[2]] = 0;
		}

	printf("#ghost state bf\n");

	found = YES;
	num = 1;
	while(found && num > 0)
	{
	    found = NO;

	    num = 0;
	    for(i=imin[0]; i<imax[0]; i++)
		for(j=imin[1]; j<imax[1]; j++)
		    for(k=imin[2]; k<imax[2]; k++)
		    {
			if(norset[i-imin[0]][j-imin[1]][k-imin[2]] != 0)
			    continue;

			found = YES;
			ic[0] = i;
			ic[1] = j;
			ic[2] = k;

			if(get_ave_state(imin,imax,ic,wave,norset,comp,ind))
			{
			    num++;
			    norset[i-imin[0]][j-imin[1]][k-imin[2]] = 2;
			}
		    }

	    printf("found = %d, num = %d\n", found, num);

	    for(i=imin[0]; i<imax[0]; i++)
		for(j=imin[1]; j<imax[1]; j++)
		    for(k=imin[2]; k<imax[2]; k++)
			if(norset[i-imin[0]][j-imin[1]][k-imin[2]] == 2)
			    norset[i-imin[0]][j-imin[1]][k-imin[2]] = 1;

	}

	free(norset);
}

#define	corner_index(p,i,gr)	irint(floor(((p)-(gr)->L[i])/(gr)->h[i]-0.5))
bool	find_block_face(
	float		*f,
	int		*icrds,
	float		*p,
	int		dir,
	RECT_GRID	*gr)
{
	int	i;

	for(i=0; i<3; i++)
	{
	    if(dir == i)
		continue;
	    
	    icrds[i] = corner_index(p[i],i,gr);
	    if(icrds[i] < -gr->lbuf[i] || icrds[i] >= gr->gmax[i]+gr->ubuf[i]-1)
		return  NO;
	    f[i] = p[i] - (gr->L[i]+(0.5+icrds[i])*gr->h[i]);
	    f[i] /= gr->h[i];
	}

	icrds[dir] = (int)((p[dir] - gr->L[dir])/gr->h[dir]);
	if(icrds[dir] < -gr->lbuf[dir] || icrds[dir] >= gr->gmax[dir]+gr->ubuf[dir]-1)
		return  NO;

	f[dir] = 0.0;
	
	return  YES;
}


//make surface point pt is inside the top grid
bool	interp_st_in_face(
	Locstate	s,
	float		*pt,
	int		dir,
	int		comp,
	INTERFACE	*intfc,
	Wave		*wave)
{
	RECT_GRID	*gr = wave->rect_grid;
	size_t		sizest = wave->sizest;
	float		f[3], *crds1, *crds2;
	int		ic[3], icrds[3], dir1, dir2, i;
	Locstate	s1, s2;
	static Locstate	st[2];
	static bool	first = YES;

	if(first)
	{
	    for(i=0; i<2; i++)
		alloc_state(intfc, &(st[i]), sizest);
	    first = NO;
	}

	if(!find_block_face(f,icrds,pt,dir,gr))
	    return NO;

	dir1 = (dir + 1)%3;
	dir2 = (dir + 2)%3;


	if(NO)
	{
	    printf("#interp_face comp=%3d  dir=%3d\n", comp, dir);
	    print_general_vector("pt", pt, 3, "\n");
	    print_general_vector("f=", f, 3, "\n");
	
	    for(i=0; i<2; i++)
	    {
		ic[dir] = icrds[dir];
		ic[dir1] = icrds[dir1] + i;
		ic[dir2] = icrds[dir2];
	    
		crds1 = Rect_coords(ic,wave);
		printf("%d ", Rect_comp(ic,wave));
		print_general_vector("crds1", crds1, 3, "\n");

		ic[dir2]++;
	    
		crds2 = Rect_coords(ic,wave);
		printf("%d ", Rect_comp(ic,wave));
		print_general_vector("crds2", crds2, 3, "\n");
	    }
	    //clean_up(ERROR);
	}

	for(i=0; i<2; i++)
	{
	    ic[dir] = icrds[dir];
	    ic[dir1] = icrds[dir1] + i;
	    ic[dir2] = icrds[dir2];
	    
	    s1 = Rect_state(ic,wave);
	    crds1 = Rect_coords(ic,wave);
	    if(Rect_comp(ic,wave) != comp)
		return NO;

	    ic[dir2]++;
	    
	    s2 = Rect_state(ic,wave);
	    crds2 = Rect_coords(ic,wave);
	    if(Rect_comp(ic,wave) != comp)
		return NO;

	    bi_interpolate_intfc_states(intfc, 1.0-f[dir2], f[dir2], 
		crds1, s1, crds2, s2, st[i]);
	}

	bi_interpolate_intfc_states(intfc, 1.0-f[dir1], f[dir1], 
		crds1, st[0], crds2, st[1], s);

	return YES;
}

void	move_pt_in_face(
	float		*crds, 
	int		*icrds,
	int		dir,
	int		nb,
	INTERFACE	*intfc)
{
	RECT_GRID	*gr = &topological_grid(intfc);
	int		i;
	float		f1, f2, tol;
	float		*h = gr->h;

	tol = 1e-6*min3(h[0], h[1], h[2]);

	for(i=0; i<3; i++)
	{
	    if(i == dir)
		continue;

	    f1 = cell_edge(icrds[i], i, gr) + tol;
	    f2 = f1 + h[i] - tol;

	    if(crds[i] < f1)
		crds[i] = f1;
	    if(crds[i] > f2)
		crds[i] = f2;
	}
}


void	crx_on_cell_face(
	float		*crdsm,		//crx point in face
	int		*im,		//face in which direcion
	int		*tfm,		//face in left or right face
	float		*crds,		//the point
	float		*nor,		//the direction normalized
	int		*icrds,		//block index, top grid
	INTERFACE	*intfc)
{
	RECT_GRID	*gr = &topological_grid(intfc);
	int		i, tf;
	float		face[3][2], crx_f, v, tmax, t, tm, hmax;
	float		*h = gr->h;

	hmax = max3(h[0], h[1], h[2]);
	
	//the coords of the direction of 6 faces.
	for(i=0; i<3; i++)
	{
	    face[i][0] = cell_edge(icrds[i], i, gr);
	    face[i][1] = face[i][0] + h[i];
	}
	
	if(NO)
	{
	    printf("#crx_on_cell_face\n");
	    print_general_vector("crds=", crds, 3, "\n");
	    print_general_vector("nor=", nor, 3, "\n");
	    for(i=0; i<3; i++)
	        printf("%3d  %15.8e  %15.8e\n", i, face[i][0], face[i][1]);
	}

	tmax = HUGE_VAL;
	for(i=0; i<3; i++)
	{
	    //if crx in the i direction is on the top (tf = 1)
	    //or on the bottom (tf = 0)
	    if(nor[i] <= 0.0)
	    {
		crx_f = face[i][0];
		tf = 0;
	    }
	    else
	    {
		crx_f = face[i][1];
		tf = 1;
	    }

	    v = crx_f - crds[i];
	    if(fabs(v) > 3.0*fabs(nor[i])*hmax)
		continue;

	    t = v/nor[i];
	    if(fabs(t) < tmax)
	    {
		tmax = fabs(t);
		tm = t;
		*tfm = tf;
		*im = i;
	    }
	}

	for(i=0; i<3; i++)
	    crdsm[i] = crds[i] + tm*nor[i];
	move_pt_in_face(crdsm, icrds, *im, *tfm, intfc);
}

bool	crx_on_major_dir(
	float		*crdsm,		//crx point in face
	int		*im,		//face in which direcion
	float		*pt,		//the point
	float		*nor,		//the direction normalized
	int		dir,
	INTERFACE	*intfc)
{
	RECT_GRID	*gr = &topological_grid(intfc);
	int		i, k, icrds[3];
	float		crds[3], crx_f, tmax;
	float		*h = gr->h;

	DEBUG_ENTER(crx_on_major_dir)

	//in which cell 
	for(i=0; i<3; i++)
	    if(i != dir)
		crds[i] = nor[i]<0.0 ? pt[i] - 0.5*h[i] : pt[i] + 0.5*h[i];
	    else
		crds[i] = pt[i];

	if(!rect_in_which(crds, icrds, gr))
	{
	    DEBUG_LEAVE(crx_on_major_dir)
	    return NO;
	}
	
	tmax = -HUGE_VAL;
	for(i=0; i<3; i++)
	{
	    if(fabs(nor[i]) > tmax)
	    {
		tmax = fabs(nor[i]);
		k = i;
	    }
	}

	crx_f = cell_edge(icrds[k], k, gr);
	if(nor[k] > 0.0)
	    crx_f += h[k];
	
	for(i=0; i<3; i++)
	    crdsm[i] = pt[i] + (crx_f-pt[k])/nor[k]*nor[i];
	move_pt_in_face(crdsm, icrds, k, 0, intfc);
	
	*im = k;
	
	if(NO)
	{
	    printf("#crx_on_major_dir  *im=%3d  crx_f=%15.8e\n", *im, crx_f);
	    print_general_vector("pt=", pt, 3, "\n");
	    print_general_vector("h=", h, 3, "\n");
	    print_general_vector("nor=", nor, 3, "\n");
	    print_general_vector("crdsm=", crdsm, 3, "\n");
	}

	DEBUG_LEAVE(crx_on_major_dir)
	return YES;
}

bool	extrp_st_normal_edge(Locstate,float*,float*,int,int,Stencil*,Front*);

bool	extrp_st_normal_edge(
	Locstate	s,
	float		*pt,
	float		*nor,
	int		comp,
	int		dir_in,
	Stencil		*sten,
	Front		*fr)
{
	INTERFACE		*intfc;
	Wave			*wave;
	float			crdsm[3];
	int			dir;
	
	intfc = fr->interf;
	wave = sten->wave;

	if(!crx_on_major_dir(crdsm, &dir, pt, nor, dir_in, intfc))
	{
	    DEBUG_LEAVE(extrp_st_along_normal)
	    return NO;
	}
	if(!interp_st_in_face(s, crdsm, dir, comp, intfc, wave))
	{
	    DEBUG_LEAVE(extrp_st_along_normal)
	    return NO;
	}

	return YES;
}



bool	extrp_st_along_normal(
	Locstate	s,
	float		*pt,
	Stencil		*sten,
	Front		*fr)
{
	INTERFACE		*intfc;
	Wave			*wave;
	RECT_GRID		*gr;
	HYPER_SURF		*hs;
	HYPER_SURF_ELEMENT	*hse;
	TRI			*tri;
	float			t[3], coords_on[3], v[3], nor[3], crds[3], crdsm[3];
	float			*h, len;
	bool         		status;
	int			dir, nb, icrds[3], new_comp;
	int			i, j;

	DEBUG_ENTER(extrp_st_along_normal)

	intfc = fr->interf;
	wave = sten->wave;
	gr = &topological_grid(intfc);
	h = gr->h;
	new_comp = sten->newcomp;
	
	if ((status = nearest_interface_point(pt,new_comp,intfc,
					      NO_SUBDOMAIN,NULL,
					      coords_on,t,&hse,&hs)) != YES)
	{
	    printf("WARNING in extrp_st_along_normal(), "
	                  "nearest_interface_point failed\n");
	    
	    DEBUG_LEAVE(extrp_st_along_normal)
	    return NO;
	}

	tri = Tri_of_hse(hse);

	for(j=0; j<3; j++)
	    nor[j] = 0.0;

	for(i=0; i<3; i++)
	{
	    normal(Point_of_tri(tri)[i], hse, hs, v, fr);
	    for(j=0; j<3; j++)
		nor[j] += v[j]*t[i];
	}

if(!debugging("revjump"))
{
	//normal is pointing comp 2
	if(new_comp == 3)
	    for(j=0; j<3; j++)
		nor[j] *= -1.0;
}
else
{
	if(new_comp == 2)
	    for(j=0; j<3; j++)
		nor[j] *= -1.0;
}
	len = Mag3d(nor);
	for(j=0; j<3; j++)
	    nor[j] /= len;

	//pt is a cell center coords, so it is a node of top grid
	for(i=0; i<3; i++)
	    crds[i] = nor[i]<0.0 ? pt[i] - 0.5*h[i] : pt[i] + 0.5*h[i];

	if(NO)
	{
	    printf("#crx_on_cell_face bf\n");
	    print_general_vector("crds=", crds, 3, "\n");
	    print_general_vector("nor=", nor, 3, "\n");
	}

	if(!rect_in_which(crds, icrds, gr))
	{
	    DEBUG_LEAVE(extrp_st_along_normal)
	    return NO;
	}

	crx_on_cell_face(crdsm, &dir, &nb, pt, nor, icrds, intfc);
	if(!interp_st_in_face(s, crdsm, dir, new_comp, intfc, wave))
	{
	    DEBUG_LEAVE(extrp_st_along_normal)
	    return NO;
	}

	DEBUG_LEAVE(extrp_st_along_normal)
}

int	get_st_from_ghost(Locstate, int*, Stencil*);

//s must have its own memory allocated.
int	get_st_from_ghost(
	Locstate	s,
	int		*icoords,
	Stencil		*sten)
{
	Wave			*wave = sten->wave;
	INTERFACE		*intfc = sten->fr->interf;
	Locstate		st;
	int			comp, ind, i, iii;
	static Gas_param	*param_list[10];
	static bool		first = YES;

	if(first)
	{
	    first = NO;
	    param_list[2] = gas_params_for_comp(2, intfc);
	    param_list[3] = gas_params_for_comp(3, intfc);
	}

	if(debugging("gfm2nd"))
	{
	    get_st_from_interp(s, icoords, sten);
	    
	    //verbose_print_state("gfm2nd", s);
	    //clean_up(ERROR);
	    
	    return 0;
	}
	
	comp = sten->newcomp;

	//print_int_vector("icrds", icoords, 3, "\n");

	set_type_of_state(s, TGAS_STATE);
	Init_params(s, param_list[comp]);
	st = Rect_state(icoords, wave);
	ind = comp == 2 ? 0 : 1;

	//print_int_vector("#icrds=", icoords, 3, "\n");
	//printf("Gd = %15.9e ", Gd(st)[ind]);
	
	Dens(s) = Gd(st)[ind];
	for(i=0; i<3; i++)
	   Vel(s)[i] = Gvel(st)[ind][i];
	Press(s) = Gp(st)[ind];
        if(g_composition_type() == MULTI_COMP_NON_REACTIVE)
        {
            for(iii=0; iii<Params(s)->n_comps; iii++)
                pdens(s)[iii] = Gpd(st)[ind][iii];
        }

	set_state(s, GAS_STATE, s);
}

LOCAL   int      fixed_curve_propagate(
        Front     *fr,
        POINTER   wave,
        CURVE     *oldc,
        CURVE     *newc,
        float     dt)
{
BOND        *oldb, *newb;
POINT       *n_pos, *o_pos;
BOND_TRI    *oldbtri, *newbtri, *newbtri1;
O_BOND      *ordered_bond;
Locstate    sl, sr; 
float       V[MAXD];
size_t      sizest = fr->sizest; 
int	    i, sav_mov;
	
	DEBUG_ENTER(fixed_curve_propagate)

	sav_mov = fr->movingframe;
	fr->movingframe = NO;

	//printf("#curve type %d  %d\n", curve_type(oldc), curve_type(newc));
	if(curve_type(newc) == NEUMANN_CURVE_P)
	{
	    //printf("#prt curve states\n");
	    //fprint_curve_states(stdout, &oldc);
	}
	
	for(oldb = oldc->first, newb = newc->first; oldb && newb; 
	    oldb = oldb->next, newb = newb->next)
	{
            ordered_bond = make_o_bond(oldb);
	    for(i=0; i<ordered_bond->num_btris; i++)
	    {
	        oldbtri = ordered_bond->_btris[i];
		newbtri = find_correspond_btri(oldbtri, newb);
		
		//this is for moving curve cases, bond_tris are not 1-1 
		if(newbtri == NULL)
		    continue;
		if(wave_type(newbtri->surface) == PASSIVE_BOUNDARY)
		    continue;
		//TMP_TK
		//if(debugging("const_bdry"))
		{
		    //TMP inflow
		    sl = left_end_btri_state(newbtri);
		    sr = right_end_btri_state(newbtri);
		    ft_assign(sl, left_end_btri_state(oldbtri), sizest);
		    ft_assign(sr, right_end_btri_state(oldbtri), sizest);
		    
                    //printf("#stype=%d %d\n", state_type(sl), state_type(sr));		    
		    if(state_type(sl)==0)
			set_type_of_state(sl, GAS_STATE);
		    
		    if(state_type(sr)==0)
			set_type_of_state(sr, GAS_STATE);

		    get_constant_state(sl, negative_component(newbtri->surface), 
				    Coords(newb->start),fr->time);
		    get_constant_state(sr, positive_component(newbtri->surface), 
				    Coords(newb->start), fr->time);
		    
		    continue;
		}

		o_pos = oldb->end;
		n_pos = newb->end;
		normal_at_point(o_pos)[0] = HUGE_VAL;
	    	
	        ft_assign(left_state(o_pos),left_end_btri_state(oldbtri), sizest);
		ft_assign(right_state(o_pos),right_end_btri_state(oldbtri), sizest);
	
	        //printf("#o_pos %d  %d n_pos %d  %d  %d", Params(left_state(o_pos)), Params(right_state(o_pos)),
		//       Params(left_state(n_pos)), Params(right_state(n_pos)), n_pos);
		
	    	if(curve_type(newc) != NEUMANN_CURVE_P)
		{
		    point_propagate(fr,wave,o_pos,n_pos,oldbtri->tri,oldbtri->surface,dt,V);
		    Coords(n_pos)[0] = Coords(o_pos)[0];
		    Coords(n_pos)[1] = Coords(o_pos)[1];
		    Coords(n_pos)[2] = Coords(o_pos)[2];
		}
		else
		{
		    if(oldb == oldc->last && !is_closed_curve(oldc))
		    {
		        //printf("#not closed ed\n");
			ft_assign(left_state(n_pos),left_state(o_pos), sizest);
		        ft_assign(right_state(n_pos),right_state(o_pos), sizest);
		    }
		    else
		    {
		        point_propagate(fr,wave,o_pos,n_pos,oldbtri->tri,oldbtri->surface,dt,V);
			if(NO && the_point(n_pos))
			{
			    print_tri(oldbtri->tri, fr->interf);
			    printf("#neg  pos %d  %d\n", 
			        negative_component(oldbtri->surface),
				positive_component(oldbtri->surface));
			}
			//printf("#cur prop comp %d\n", component(Coords(n_pos), fr->interf));
			//inside the wall, should be moved ouside
			if(gas_params_for_comp(component(Coords(n_pos), fr->interf), fr->interf) != NULL)
		            point_propagate_along_wall(fr,wave,o_pos,oldb,oldc, 
			         oldbtri->tri,oldbtri->surface,n_pos,dt,V);
		
			//printf("#addr newp  %d\n", n_pos );
			//print_general_vector("n_pos", Coords(n_pos), 3, "\n");
			
			if(NO && the_point(n_pos))
			{
			    printf("#n_pos found\n");
			    print_general_vector("oldp", Coords(o_pos), 3, "\n");
			}
		    }
		}

		sl = left_end_btri_state(newbtri);
		sr = right_end_btri_state(newbtri);
		ft_assign(sl, left_state(n_pos), sizest);
		ft_assign(sr, right_state(n_pos), sizest);

	    }
	    if(curve_type(newc) != NEUMANN_CURVE_W && 
	       curve_type(newc) != NEUMANN_CURVE_P)
	        reset_fixed_btri_states(newb, ordered_bond, fr);
	}

	//for closed curve, the start state of the first bond comes from the end 
	//state of the last bond. 
	//for none closed curve, ASSUMING the node is always on the subdomain bdry
	//so it will be removed after scatter_front, just fill the state with the state
	//from last time step.

        ordered_bond = make_o_bond(oldc->first);
	for(i=0; i<ordered_bond->num_btris; i++)
	{
	    oldbtri = ordered_bond->_btris[i];
	    newbtri = find_correspond_btri(oldbtri, newc->first);
	    
	    //this is for moving curve cases, bond_tris are not 1-1 
	    if(newbtri == NULL)
	        continue;
	    if(wave_type(newbtri->surface) == PASSIVE_BOUNDARY)
		continue;

	    newbtri1 = find_correspond_btri(oldbtri, newc->last);

	    //newbtri: btri of the first bond.
	    //newbtri1: btri of the last bond.

	    if(newbtri1 == NULL)
	    {
	        printf("ERROR: fixed_curve_propagate, interface inconsistent for btri.\n");
		clean_up(ERROR);
	    }
	    if(newbtri1->surface != newbtri->surface)
	    {
	        printf("ERROR: fixed_curve_propagate, interface inconsistent for surface.\n");
		clean_up(ERROR);
	    }

	    if(is_closed_curve(oldc))
	    {
	        sl = left_end_btri_state(newbtri1);
	        sr = right_end_btri_state(newbtri1);
	        ft_assign(left_start_btri_state(newbtri), sl, sizest);
		ft_assign(right_start_btri_state(newbtri), sr, sizest);
	    }
	    else
	    {
	        sl = left_start_btri_state(oldbtri);
	        sr = right_start_btri_state(oldbtri);
	        ft_assign(left_start_btri_state(newbtri), sl, sizest);
		ft_assign(right_start_btri_state(newbtri), sr, sizest);
	    }
	}

	DEBUG_LEAVE(fixed_curve_propagate)

	fr->movingframe = sav_mov;
        
	return GOOD_CURVE; 
}   /*end fixed_curve_propagate */


IMPORT  void	smooth_curve(CURVE *, float);

//#bjet2  
LOCAL   int      subdomain_curve_propagate(
        Front     *fr,
        POINTER   wave,
        CURVE     *oldc,
        CURVE     *newc,
        float     dt)
{
	BOND        *oldb, *newb;
	BOND_TRI    *oldbtri, *newbtri;
	size_t      sizest = fr->sizest; 
	int	    i, k;
	INTERFACE   *intfc;
	RECT_GRID   *gr;
	POINT	    *oldp, *newp;
	Locstate    sl, sr;
	float	    V[3], nor[4], hmin, bdry_fac = 0.5;

	DEBUG_ENTER(subdomain_curve_propagate)
	
	//printf("#curve type subdomain %d  %d\n", curve_type(oldc), curve_type(newc));
	intfc = fr->interf;
	gr = computational_grid(intfc);
	hmin = min3(gr->h[0], gr->h[1], gr->h[2]);

	for(oldb = oldc->first, newb = newc->first; 
	    oldb && newb; 
	    oldb = oldb->next, newb = newb->next)
	{
	    if(Btris(oldb) == NULL || Btris(newb) == NULL)
	    {
	        printf("ERROR: subdomain_curve_propagate. NULL Btris\n");
		clean_up(ERROR);
	    }
	    //subdomain curves have only one bond tri.
	    oldbtri = Btris(oldb)[0];
	    newbtri = Btris(newb)[0];
	    if(oldbtri == NULL || newbtri == NULL)
	    {
	        printf("ERROR: subdomain_curve_propagate. NULL bond_tri\n");
		clean_up(ERROR);
	    }
	    //subdomain curves states will be discarded after scatter_front
	    ft_assign(left_start_btri_state(newbtri), 
	    	left_start_btri_state(oldbtri), sizest);
	    ft_assign(right_start_btri_state(newbtri), 
	    	right_start_btri_state(oldbtri), sizest);
	    ft_assign(left_end_btri_state(newbtri), 
	    	left_end_btri_state(oldbtri), sizest);
	    ft_assign(right_end_btri_state(newbtri), 
	    	right_end_btri_state(oldbtri), sizest);

	    //moving the starting point.
	    oldp = oldb->start;
	    newp = newb->start;
	    if(point_outside_open_bdry(&k, nor, oldp, intfc))
	    {
		point_propagate(fr,wave,oldp,newp,
		    Hyper_surf_element(oldbtri->tri), 
		    Hyper_surf(oldbtri->surface),dt,V);
		slsr(newp, Hyper_surf_element(newbtri->tri), 
			   Hyper_surf(newbtri->surface), &sl, &sr);
		ft_assign(sl,left_state(newp),sizest);
		ft_assign(sr,right_state(newp),sizest);
		
		if(V[k]*nor[k] < 0.0 || nor[3] > bdry_fac)
		    ft_assign(Coords(newp), Coords(oldp), 3*FLOAT);
		
		//printf("#reach open bdry.\n");
		
		if(the_point(newp))
		{
		    verbose_print_state("sub_sl", sl);
		    verbose_print_state("sub_sr", sr);

		    printf("dt = %24.15e p = %d op = %d\n",dt, newp, oldp);
		    printf("cnum %d  slsr %d  %d  %d  %d\n", 
		    	curve_number(oldc),left_start_btri_state(oldbtri), 
			right_start_btri_state(oldbtri), 
			left_end_btri_state(oldbtri),
			right_end_btri_state(oldbtri));

		    print_general_vector("oldp", Coords(oldp), 3, "\n");
		    print_tri(newbtri->tri, intfc);
		}
	    }

	    if(newb->next != NULL)
		continue;
	    if(is_closed_curve(newc))
		continue;
	    
	    //moving the ending point.
	    newp = newb->end;
	    oldp = oldb->end;
	    if(point_outside_open_bdry(&k, nor, oldp, intfc))
	    {
		point_propagate(fr,wave,oldp,newp,
		    Hyper_surf_element(oldbtri->tri), 
		    Hyper_surf(oldbtri->surface),dt,V);
		slsr(newp, Hyper_surf_element(newbtri->tri), 
			   Hyper_surf(newbtri->surface), &sl, &sr);
		ft_assign(sl,left_state(newp),sizest);
		ft_assign(sr,right_state(newp),sizest);
		
		if(V[k]*nor[k] < 0.0 || nor[3] > bdry_fac)
		    ft_assign(Coords(newp), Coords(oldp), 3*FLOAT);
	    }
	}

	//if(fr->step%20 == 0)
	//    smooth_curve(newc, hmin*0.01);

	DEBUG_LEAVE(subdomain_curve_propagate)
		
        return GOOD_CURVE; 
}   /*end fixed_curve_propagate */

LOCAL   SURFACE  *find_physical_surface_at_curve(
	CURVE       *c,
	ORIENTATION *orient)
{
	SURFACE     **s, **surf;
	int         i;
	ORIENTATION s_or;


	for (i = 0, s = c->pos_surfaces, s_or = POSITIVE_ORIENTATION;
	     i < 2;
	     ++i,   s = c->neg_surfaces, s_or = NEGATIVE_ORIENTATION)
	{
    	    for (surf = s; surf && *surf; ++surf)
            {
                 if (wave_type(*surf) >= FIRST_PHYSICS_WAVE_TYPE)
                 {
                      *orient = s_or;
                      return *surf;
                 }
             }
 	
	}
        return NULL; 		
}

LOCAL   SURFACE  *adjacent_surface(
        SURFACE          *surf,
        ORIENTATION      s_orient,
        CURVE            *curve,
        ANGLE_DIRECTION  angle_dir,
        ORIENTATION      *adj_s_orient)
{
        INTERFACE   *intfc = surf->interface;
        BOND        *bond = curve->first;
        CURVE       *c;
        SURFACE     **s, *ans = NULL; 
        BOND_TRI    **btris;
        const TRI   *tri,*tri1,*tri2;
        float       sin12,cos12,oldsin12,oldcos12;
	const float *t1,*t2;
	int         i, dim = curve->interface->dim;
	COMPONENT   test_comp;
        static int nfail = 0;
 
        if (debugging("adjacent_surface"))
        {
            print_surface(surf);
            print_orientation("s_orient = ",s_orient,", ");
            print_angle_direction("angle_dir = ",angle_dir,"\n");
        }
        printf("Entered adjacent_surface()\n");
        test_comp = (surf_ang_oriented_l_to_r(angle_dir,s_orient)) ?
                                positive_component(surf) :
                                negative_component(surf);

	for(btris =  Btris(bond); btris && *btris;++btris)
	{
	   tri = (*btris)->tri;
	   if((*btris)->surface == surf)
	   {
               tri = (*btris)->tri;
	       t1 = Tri_normal(tri);	   
	   }
	}

        for(s == curve->pos_surfaces; s && *s; ++s)
	{
            if (*s == surf && s_orient == NEGATIVE_ORIENTATION)
                continue;

                /* Test for consistent component */
                                                                                
            if (((angle_dir == CLOCKWISE) &&
                    !matchable_comps(positive_component(*s),test_comp,intfc))
                                 ||
                ((angle_dir == COUNTER_CLOCK) &&
                    !matchable_comps(negative_component(*s),test_comp,intfc)))
                     continue;
            for(btris =  Btris(bond); btris && *btris;++btris)
            {
                if((*btris)->surface == *s)
                {
                    tri = (*btris)->tri; 
                    t2 = Tri_normal(tri);
                }
            }
            (void) vector_product(t1,t2,&sin12,dim);
            cos12 = scalar_product(t1,t2,dim);
            if (ans == NULL)
            {
                oldsin12 = sin12;
                oldcos12 = cos12;
                ans = *s;
                *adj_s_orient = NEGATIVE_ORIENTATION;
                printf("in adjacent_surface(), test1\n" );
                print_surface(ans);
                continue;
            }
            if (is_new_angle_smaller(sin12,cos12,oldsin12,oldcos12,angle_dir))
            {
                oldsin12 = sin12;
                oldcos12 = cos12;
                ans = *s;
                *adj_s_orient = NEGATIVE_ORIENTATION;
                printf("in adjacent_surface(), test2\n" );
                print_surface(ans);
            }
        }
        
        for(s == curve->neg_surfaces; s && *s; ++s)
        {
            if (*s == surf && s_orient == POSITIVE_ORIENTATION)
                continue;
                                                                                                                            
                /* Test for consistent component */
                                                                                                                             
            if (((angle_dir == CLOCKWISE) &&
                    !matchable_comps(negative_component(*s),test_comp,intfc))
                                 ||
                ((angle_dir == COUNTER_CLOCK) &&
                    !matchable_comps(positive_component(*s),test_comp,intfc)))
                     continue;
            for(btris =  Btris(bond); btris && *btris;++btris)
            {
                if((*btris)->surface == *s)
                {
                    tri = (*btris)->tri;
                    t2 = Tri_normal(tri);
                }
            }
            (void) vector_product(t1,t2,&sin12,dim);
            cos12 = scalar_product(t1,t2,dim);
            if (ans == NULL)
            {
                oldsin12 = sin12;
                oldcos12 = cos12;
                ans = *s;
                *adj_s_orient = POSITIVE_ORIENTATION;
                continue;
                printf("in adjacent_surface(), test3\n" );
                print_surface(ans);

            }
            if (is_new_angle_smaller(sin12,cos12,oldsin12,oldcos12,angle_dir))
            {
                oldsin12 = sin12;
                oldcos12 = cos12;
                ans = *s;
                *adj_s_orient = POSITIVE_ORIENTATION;
                printf("in adjacent_surface(), test4\n" );
                print_surface(ans);
            }
        }
        
        if (ans == NULL)
        {
            if (nfail++ < 10) /* TOLERANCE */
                (void) printf("WARNING in adjacent_surface(), returning null\n");
            else
            {
                screen("ERROR in adjacent_surface(), "
                       "can't find adjacent surface\n");
                clean_up(ERROR);
            }
        }
        else
            nfail = 0;
        debug_print("adjacent_surface","Leaving adjacent_surface(), ans = %d\n",ans);

        return  ans;           
} 

LOCAL   bool    matchable_comps(
        COMPONENT comp1,
        COMPONENT comp2,
        INTERFACE *intfc)
{
        if (equivalent_comps(comp1,comp2,intfc))
            return YES;
        else if (is_exterior_comp(comp1,intfc) && is_excluded_comp(comp2,intfc))            return YES;
        else if (is_exterior_comp(comp2,intfc) && is_excluded_comp(comp1,intfc))            return YES;
        else
            return NO;
}               /*end matchable_comps*/

LOCAL  	O_BOND 	*make_o_bond(
	BOND		*b)
{
	float		x_axis[3], y_axis[3], z_axis[3];
	float		vector_in_xz_plane[3], uni_array[3];
	float		*origin = Coords(b->start);
	float		magx, magy;
	float		x_proj, y_proj;
	int		i, j_start, j_end, other;
	int		num_btris = 0;
	int		not_ordered;
	TRI		*tri;
	O_BOND		*ordered_bond;
	static int  	third_point[3][3] = 
	{
	    {-1, 2, 1},
	    { 2,-1, 0},
	    { 1, 0,-1}
	};

	static int 	orient_map[3][3] =
	{
	    { 0, 1,-1},
	    {-1, 0, 1},
	    { 1,-1, 0}
	};

	DEBUG_ENTER(make_o_bond)
	while (Btris(b)[num_btris])
	    num_btris++;

	ordered_bond = (O_BOND*) Store(sizeof(O_BOND));
	ordered_bond->_btris = (BOND_TRI**) Store(num_btris*sizeof(BOND_TRI*));
	ordered_bond->angle = (float*) Store(num_btris*sizeof(float));
	ordered_bond->orientation = (int*) Store(num_btris*sizeof(int));

	for (i = 0; i < num_btris; i++)
	    ordered_bond->_btris[i] = Btris(b)[i];
	ordered_bond->num_btris = num_btris;
	ordered_bond->old_bond = b;

	/* Construct local coordinate system. *
	 * x_axis defined to be zero radians  *
	 *  and in plane of first tri,        *
	 * z_axis defined to be bond b,       *
	 * y_axis = z_axis X x_axis           */
	
	tri = Btris(b)[0]->tri;  
	j_start = j_end = 0;
	
	while (Point_of_tri(tri)[j_start] != b->start && j_start < 3)
	    j_start++;
	while (Point_of_tri(tri)[j_end] != b->end && j_start < 3)
	    j_end++;
	other = third_point[j_start][j_end];

	difference3d(Coords(b->end),origin,z_axis);
	difference3d(Coords(Point_of_tri(tri)[other]),origin,
		     vector_in_xz_plane);
	magy = vector_product(z_axis,vector_in_xz_plane,y_axis,3);
	magx = vector_product(y_axis,z_axis,x_axis,3);
	vector_scale3d(y_axis,(1.0/magy));
	vector_scale3d(x_axis,(1.0/magx));


	for (i = 0; i < num_btris; i++)
	{
	    tri = Btris(b)[i]->tri;
	    j_start = j_end = 0;
	    
	    while (Point_of_tri(tri)[j_start] != b->start && j_start < 3)
	        j_start++;
	    while (Point_of_tri(tri)[j_end] != b->end && j_start < 3)
	        j_end++;
	    ordered_bond->orientation[i] = orient_map[j_start][j_end];
	    other = third_point[j_start][j_end];
	    
	    if (ordered_bond->orientation[i] == 0)
	    {
	        screen("unable to orient tri wrt bond in make_o_bond()\n");
		clean_up(ERROR);
	    }

	    difference3d(Coords(Point_of_tri(tri)[other]),origin,uni_array);

	    x_proj = scalar_product(uni_array,x_axis,3);
	    y_proj = scalar_product(uni_array,y_axis,3);
      
/*	    if (sqr(x_proj) + sqr(y_proj) <= .000001) //TOLERANCE
	    {
	        screen("ERROR in make_o_bond(), degenerate TRI\n");
		print_tri(tri,Btris(b)[i]->surface->interface);
		clean_up(ERROR);
	    }
*/
	    ordered_bond->angle[i] = normalized_angle(atan2(y_proj,x_proj));
	}
        
	ordered_bond->angle[0] = 0.0;

	if(fixed_curve_debug)
	{
	    printf("IN make_o_bond(),before order the bond\n");
	    print_bond(b);
	    for(i = 0; i < num_btris; i++)
	    {
		printf("btris[%d] on the surface %llu,angle = %f\n",
			i,surface_number(Btris(b)[i]->surface),ordered_bond->angle[i]);    
		//print_tri_coords(Btris(b)[i]->tri);
	    }
	}
	not_ordered = YES;
	while (not_ordered)  /* TODO: this should be rewritten using qsort() */
        {
	    not_ordered = NO;
	    for (i = 0; i < num_btris - 1; i++)
	        if (ordered_bond->angle[i+1] < ordered_bond->angle[i])
		{
		    float ftmp;
		    int   itmp;
		    BOND_TRI *btmp;

		    ftmp = ordered_bond->angle[i];
		    btmp = ordered_bond->_btris[i];
		    itmp = ordered_bond->orientation[i];

		    ordered_bond->angle[i] = ordered_bond->angle[i+1];
		    ordered_bond->_btris[i] = ordered_bond->_btris[i+1];
		    ordered_bond->orientation[i] = 
		        ordered_bond->orientation[i+1];

		    ordered_bond->angle[i+1] = ftmp;
		    ordered_bond->_btris[i+1] = btmp;
		    ordered_bond->orientation[i+1] = itmp;

		    not_ordered = YES;
		} 
	}
        DEBUG_LEAVE(make_o_bond)
	return ordered_bond;
}
		/*end make_o_bond*/

LOCAL	void	difference3d(
	float *a,
	float *b,
	float *c)
{
	c[0] = a[0] - b[0];
	c[1] = a[1] - b[1];
	c[2] = a[2] - b[2];
}		/*end difference3d*/

LOCAL	void	vector_scale3d(
	float *a,
	float s)
{
	a[0] *= s;
	a[1] *= s;
	a[2] *= s;
}		/*end vector_scale3d*/

LOCAL   void set_curve_vel_to_zero(CURVE *c)
{
        BOND         *b;
        BOND_TRI     **btris;
        int          i;
        Locstate     start_left,start_right,end_left,end_right;
        INTERFACE    *intfc = c->interface;
                                                                                           
        for(b = c->first; b; b=b->next)
        {
            for(btris = Btris(b); btris && *btris; btris++)
            {
                start_left = left_start_btri_state(*btris);
                start_right = right_start_btri_state(*btris);
                end_left = left_end_btri_state(*btris);
                end_right = right_end_btri_state(*btris);
                for(i = 0; i < 3;i++)
                {
                    Vel(start_left)[i] = 0.0;
                    Vel(start_right)[i] = 0.0;
                    Vel(end_left)[i] = 0.0;
                    Vel(end_right)[i] = 0.0;
                }
            }
        }
}

#endif /* defined(THREED) */
