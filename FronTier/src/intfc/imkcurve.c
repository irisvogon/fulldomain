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
*				imkcurve.c:
*
*       Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Containing function of rebuilding interface within a mesh block.
*
*/


#define DEBUG_STRING "i_make_curve"

#include <intfc/int.h>

struct  _EG_CRX_2D {
        BBI_POINT ***x_crx;
        BBI_POINT ***y_crx;
        BBI_POINT ***x_curve_crx;
        BBI_POINT ***y_curve_crx;
        BBI_POINT ***node_crx;
        BBI_POINT *crx_store;
        COMPONENT **comp;
        CURVE     **nodes;
        int       num_nodes;
};
typedef struct _EG_CRX_2D EG_CRX_2D;

struct _BLK_CRX_2D {
        COMPONENT       **comp;
        int             **ix;
        int             **iy;
        int             num_comps;
        int             nv[4];
        COMPONENT       comps[4];
        BBI_POINT       ***crx;
        BBI_POINT       *node_crx;
        BBI_POINT       *crx_store;
        BLK_INFO        *blk_info;
        BLK_TYPE        blk_type;
};
typedef struct _BLK_CRX_2D BLK_CRX_2D;

LOCAL	bool onfront_block2d(int,int,const EG_CRX_2D*);
LOCAL	bool increment_theta(float*,float,float,float,ELLIP_PARAMS*,
			float*,float);
LOCAL	void reset_domain_comp2d(COMPONENT**,RECT_GRID*);
LOCAL   void assign_two_comp_domain2d(float (*func)(POINTER,float*),POINTER,
                        COMPONENT**,RECT_GRID*,COMPONENT,COMPONENT);
LOCAL	void make_grid_curves(BLK_CRX_2D*,const EG_CRX_2D*,int*,BOND**,int*);
LOCAL 	void assign_blk_crx2d(BLK_CRX_2D*,int,int,const EG_CRX_2D*);
LOCAL	void connect_adj_blk(BLK_BOND*,BLK_BOND*);
LOCAL	int  count_crx_through_comp2d(int*,COMPONENT**);
LOCAL   int  install_grid_crx2d(float (*func)(POINTER,float*),POINTER,
                        EG_CRX_2D*,RECT_GRID*);
LOCAL	int  construct_comp2_blk2d(BLK_CRX_2D*,BLK_BOND*);
LOCAL	BLK_CRX_2D *alloc_blk_crx2d(bool);

/*	Initialization for level function parameters */
LOCAL	POINTER init_line_params(RECT_GRID*);
LOCAL	POINTER init_multi_circle_params(RECT_GRID*);
LOCAL	POINTER init_sine_wave_params(RECT_GRID*);
LOCAL	POINTER init_test_disk_params(RECT_GRID*);

/*	Level functions for curves */

LOCAL   float line_func(POINTER,float*);
LOCAL   float multi_circle_func(POINTER,float*);
LOCAL   float test_disk_func(POINTER,float*);

#define		MAX_NUM_SEGMENTS		100


EXPORT CURVE *make_elliptic_curve(
	ELLIP_PARAMS *ellip,
	COMPONENT compin,
	COMPONENT compout,
	float space)
{
	int	    i;
	bool	    done;
	ORIENTATION nor_or = ellip->nor_orient;
	int	    closed = ellip->closed;
	float	    theta[MAXD-1];
	RECT_GRID   *gr = ellip->gr;
	float	    *h = gr->h, *L = gr->L, *U = gr->U;
	float	    coords[MAXD];
	float	    sgn;
	float	    ThetaS = ellip->ThetaS[0], ThetaE = ellip->ThetaE[0];
	CURVE	    *cur;
	NODE	    *ns, *ne;
	COMPONENT   l_comp,r_comp;

	if (nor_or == NEGATIVE_ORIENTATION)
	{
	    sgn = -1.0;
	    l_comp = compout;
	    r_comp = compin;
	}
	else
	{
	    sgn = 1.0;
	    l_comp = compin;
	    r_comp = compout;
	}
	if (closed == YES)
	{
	    if (nor_or == NEGATIVE_ORIENTATION)
	    {
		ThetaS = normalized_angle(ThetaS);
		ThetaE = ThetaS - 2.0*PI;
	    }
	    else
	    {
		ThetaS = normalized_angle(ThetaS);
		ThetaE = ThetaS + 2.0*PI;
	    }
	    coords_on_ellips(ellip->ThetaS,coords,ellip);
	    ns = make_node(Point(coords));
	    ne = ns;
	}
	else
	{
	}
	cur = make_curve(l_comp,r_comp,ns,ne);

	done = NO;
	theta[0] = ThetaS;
	while (done == NO)
	{
	    done = increment_theta(theta,ThetaS,ThetaE,sgn,ellip,h,space);
	    if (done == YES)
		break;
	    coords_on_ellips(theta,coords,ellip);
	    if (insert_point_in_bond(Point(coords),cur->last,cur) !=
		FUNCTION_SUCCEEDED)
	    {
	        screen("ERROR in g_make_ellipse(), "
		       "insert_point_in_bond() failed\n");
	        clean_up(ERROR);
	    }
	}
	return cur;
}	/* end make_elliptic_curve */


EXPORT	void coords_on_ellips(
	float		*theta,
	float		*coords,
	ELLIP_PARAMS	*ellip)
{
	float		*cen = ellip->cen, *rad = ellip->rad;
	FOURIER_POLY	*fpoly = ellip->fpoly;
	LEGENDRE_POLY	*lpoly = ellip->lpoly;
	int		dim = ellip->dim;
	float		Dr[MAXD];
	float		er[3];
	float		r;
	int		i, j;

	switch (dim)
	{
	case 2:
	    er[0] = cos(theta[0]);
	    er[1] = sin(theta[0]);
	    r = 1.0/hypot(er[0]/rad[0],er[1]/rad[1]);
	    break;
#if defined(THREED)
	case 3:
	    er[0] = cos(theta[0])*sin(theta[1]);
	    er[1] = sin(theta[0])*sin(theta[1]);
	    er[2] = cos(theta[1]);
	    r = 1.0/sqrt(sqr(er[0]/rad[0])+sqr(er[1]/rad[1])+sqr(er[2]/rad[2]));
	    break;
#endif /* defined(THREED) */
	default:
	    r = ERROR_FLOAT;
	    screen("ERROR in coords_on_ellips(), "
	           "invalid dimension %d\n",dim);
	    clean_up(ERROR);
	}
	if (fpoly != NULL)
	    r += fourier_poly(theta,fpoly);
	if ((dim == 2) && (lpoly != NULL))
	    r += legendre_poly(er[1],lpoly);

	for (i = 0; i < dim; i++)
            Dr[i] = 1.0*r*er[i];
	for (i = 0; i < dim; i++)
	    coords[i] = cen[i] + Dr[i];
}	/*end coords_on_ellips */


/*
*			increment_theta():
*
*	Computes the incremented theta in the polar coordinated
*	of the z = constant cross section of an ellipsoid.
*	The difference new_theta - theta is determined so that
*	the displacement of the two consecutive points on the
*	ellipsoid are approximately separated by the distance
*	space in the scaled metric with scale vector h.
*	This function assumes that if sgn and theta_e - theta_s
*	have the same algebraic sign.
*/

LOCAL	bool increment_theta(
	float		*theta,
	float		theta_s,
	float		theta_e,
	float		sgn,
	ELLIP_PARAMS	*ellip,
	float		*h,
	float		space)
{
	float		d_theta, dl;
	float		new_theta;
	float		rad[MAXD], r;
	float		er[2];
	bool		done;
	rad[0] = ellip->rad[0];	rad[1] = ellip->rad[1];
	dl = space*hypot(h[0],h[1]);
	er[0] = cos(*theta);	er[1] = sin(*theta);
	r = 1.0/hypot(er[0]/rad[0],er[1]/rad[1]);
	d_theta = dl/(r*r*r*hypot(er[0]/sqr(rad[0]),er[1]/sqr(rad[1])));

	new_theta = theta[0] + sgn*d_theta;
	done = NO;
	if (!Between(new_theta,theta_s,theta_e))
	{
	    new_theta = theta_e;
	    done = YES;
	}

	*theta = new_theta;
	return done;
}		/*end increment_theta*/

EXPORT float multi_sine_mode_func(
	POINTER func_params, 
        float *coords)
{
	FOURIER_POLY *sine_params = (FOURIER_POLY*)func_params;
	float z0,dist;
	float *nu,A,phase;	
	float *L = sine_params->L;
	float *U = sine_params->U;
	int i,j;
	int dim = sine_params->dim;

	
	z0 = sine_params->z0;
	for (i = 0; i < sine_params->num_modes; ++i)
	{
	    A = sine_params->A[i];
	    nu = sine_params->nu[i];
	    phase = sine_params->phase[i];
	    for (j = 0; j < dim-1; ++j)
	    {
		phase += nu[j]*(coords[j] - L[j])*2.0*PI/(U[j] - L[j]);
	    }
	    z0 += A*sin(phase);
	}
	dist = coords[dim-1] - z0;
	return dist;
}	/* end multi_sine_mode_func */


/*******************************************************************
*       This function make a curve described by the function       *
*       func = 0. The negative side of the curve has right_comp    *
*       and the positive side of the curve has left_comp.          *
*******************************************************************/

EXPORT CURVE *make_level_curve(
	RECT_GRID   *rgr,
	INTERFACE   *intfc,
	COMPONENT   left_comp,
	COMPONENT   right_comp,
	float       (*func)(POINTER,float*),
        POINTER     func_params,
	bool	    use_dual_grid)
{
	int		i,num_crx,num_segs,np,*gmax;
	RECT_GRID	dual_gr,*grid;
	CURVE		*curve;
	EG_CRX_2D	Eg_crx;
	BLK_INFO	blk_info;
	static BLK_CRX_2D  *blk_crx;
	float		coords[2] = {0.0,0.0};
	NODE		*ns,*ne;
	BOND 		*b,*first[MAX_NUM_SEGMENTS];
        INTERFACE       *save_intfc;

        save_intfc = current_interface();
        set_current_interface(intfc);

	if (blk_crx == NULL)
	    blk_crx = alloc_blk_crx2d(NO);

	zero_scalar(&Eg_crx,sizeof(EG_CRX_2D));
	if (use_dual_grid)
	{
	    set_grid_for_surface_construction(&dual_gr,rgr);
	    grid = &dual_gr;
	}
	else
	    grid = rgr;
	gmax = grid->gmax;
	bi_array(&Eg_crx.comp,gmax[0]+1,gmax[1]+1,sizeof(COMPONENT));

	reset_domain_comp2d(Eg_crx.comp,grid);
	assign_two_comp_domain2d(func,func_params,Eg_crx.comp,
				grid,right_comp,left_comp);
	num_crx = count_crx_through_comp2d(gmax,Eg_crx.comp);
	if (num_crx == 0)
	{
	    intfc->default_comp = Eg_crx.comp[0][0];
	    free_these(1,Eg_crx.comp);
	    return NULL;
	}

	bi_array(&Eg_crx.x_crx,gmax[0],gmax[1]+1,sizeof(BBI_POINT*));
        bi_array(&Eg_crx.y_crx,gmax[0]+1,gmax[1],sizeof(BBI_POINT*));
        uni_array(&Eg_crx.crx_store,num_crx,sizeof(BBI_POINT));

	num_crx = install_grid_crx2d(func,func_params,&Eg_crx,grid);

	ns = make_node(Point(coords));
	ne = make_node(Point(coords));
	curve = make_curve(left_comp,right_comp,ns,ne);

	curve->first = curve->last = NULL;

	blk_info.num_surfs = 1;
	uni_array(&blk_info.curves,1,sizeof(CURVE*));
	uni_array(&blk_info.cur_bonds,1,sizeof(BOND*));

	blk_info.cur_bonds[0] = NULL;
	blk_info.curves[0] = curve;
	curve->num_points = 0;

	blk_crx->comps[0] = right_comp;
	blk_crx->comps[1] = left_comp;
	blk_crx->blk_info = &blk_info;

	for (i = 0; i < num_crx; ++i)
	    Eg_crx.crx_store[i].c = curve;

	make_grid_curves(blk_crx,&Eg_crx,gmax,first,&num_segs);
	curve->num_points = 2;
	for (b = curve->first; b->next != NULL; b = b->next)
	    curve->num_points++;
	curve->last = b;
	if (curve->start->posn == curve->end->posn)
	{
	    change_node_of_curve(curve,NEGATIVE_ORIENTATION,curve->start);
	    delete_node(ne);
	}
	if (num_segs != 1)
	{
	    CURVE *c;
	    for (i = 0; i < num_segs; ++i)
	    {
	    	if (first[i]->start == curve->start->posn)
		    continue;
		np = 2;
		for (b = first[i]; b->next != NULL; b = b->next) ++np;
		if (b->end == first[i]->start)
		{
		    ns = ne = make_node(first[i]->start);
		}
		else
		{
		    ns = make_node(first[i]->start);
		    ne = make_node(b->end);
		}
		c = make_curve(left_comp,right_comp,ns,ne);
		c->first = first[i];
		c->last = b;
		c->num_points = np;
	    }
	}

	free_these(3,Eg_crx.x_crx,Eg_crx.y_crx,Eg_crx.crx_store);
	free_these(3,Eg_crx.comp,blk_info.curves,blk_info.cur_bonds);
	set_current_interface(save_intfc);
	if (num_segs == 1) return curve;
	else return NULL;
}	/* end make_level_surface */


/**********************************************************************
*	This function reset components of the domain to NO_COMP       *
**********************************************************************/

LOCAL	void reset_domain_comp2d(
	COMPONENT **comp,
	RECT_GRID *gr)
{
	int i,j;
	int *gmax = gr->gmax;

	for (i = 0; i <= gmax[0]; ++i)
	{
	    for (j = 0; j <= gmax[1]; ++j)
	    {
	    	comp[i][j] = NO_COMP;
	    }
	}
}	/* end reset_domain_comp */


LOCAL	int count_crx_through_comp2d(
	int *gmax,
	COMPONENT **comp)
{
	int i,j,num_crx,num_comp;

	num_crx = 0;
	for (i = 0; i <= gmax[0]; ++i)
	{
	    for (j = 0; j <= gmax[1]; ++j)
	    {
		if (i != gmax[0])
		{
		    if (comp[i][j] != comp[i+1][j])
			    ++num_crx;
		}
		if (j != gmax[1])
		{
		    if (comp[i][j] != comp[i][j+1])
			    ++num_crx;
		}
	    }
	}
	return num_crx;
}	/* end count_crx_through_comp */


LOCAL	void assign_two_comp_domain2d(
	float (*func)(POINTER,float*),
	POINTER func_params,
	COMPONENT **comp,
	RECT_GRID *gr,
	COMPONENT left_comp,
	COMPONENT right_comp)
{
	int i,j,k;
	int *gmax = gr->gmax;
	float *L = gr->L;
	float *h = gr->h;
	float coords[MAXD];

	for (i = 0; i <= gmax[0]; ++i)
	{
	    coords[0] = L[0] + i*h[0];
	    for (j = 0; j <= gmax[1]; ++j)
	    {
	    	coords[1] = L[1] + j*h[1];
		if (comp[i][j] != NO_COMP) continue;
		if ((*func)(func_params,coords) > 0)
		    comp[i][j] = left_comp;
		else 
		    comp[i][j] = right_comp;
	    }
	}
}	/* end assign_two_comp_domain */


LOCAL	int install_grid_crx2d(
	float (*func)(POINTER,float*),
        POINTER func_params,
	EG_CRX_2D *eg_crx,
	RECT_GRID *grid)
{
	float coords1[MAXD];
	float coords2[MAXD];
	float crds_crx[MAXD];
	float *L = grid->L;
	float *h = grid->h;
	int *gmax = grid->gmax;
	int dim = grid->dim;
	int i,j,n_crx;
	BBI_POINT ***x_crx = eg_crx->x_crx;
	BBI_POINT ***y_crx = eg_crx->y_crx;
	BBI_POINT *crx_store = eg_crx->crx_store;
	COMPONENT **comp = eg_crx->comp;

	n_crx = 0;

	/* install x-crossings */

	for (j = 0; j <= gmax[1]; ++j)
	{
	    coords1[1] = coords2[1] = L[1] + j*h[1];
	    for (i = 0; i < gmax[0]; ++i)
	    {
		x_crx[i][j] = NULL;
		if (comp[i][j] != comp[i+1][j])
		{
		    coords1[0] = L[0] + i*h[0];
		    coords2[0] = L[0] + (i+1)*h[0];
		    if (! grid_line_crx_in_dir(func,func_params,
				dim,coords1,coords2,crds_crx,0))
		    {
			screen("ERROR: in install_grid_crx(), no x-crxing!");
			clean_up(ERROR);
		    }
		    if (crds_crx[0] - coords1[0] < 0.004*h[0])
			crds_crx[0] = coords1[0] + 0.004*h[0];
		    if (coords2[0] - crds_crx[0] < 0.004*h[0])
			crds_crx[0] = coords2[0] - 0.004*h[0];
		    crx_store[n_crx].p = Point(crds_crx);
		    x_crx[i][j] = &crx_store[n_crx++];
		}
	    }
	}

	/* install y-crossings */

	for (i = 0; i <= gmax[0]; ++i)
	{
	    coords1[0] = coords2[0] = L[0] + i*h[0];
	    for (j = 0; j < gmax[1]; ++j)
	    {
		y_crx[i][j] = NULL;
		if (comp[i][j] != comp[i][j+1])
		{
		    coords1[1] = L[1] + j*h[1];
		    coords2[1] = L[1] + (j+1)*h[1];
		    if (!grid_line_crx_in_dir(func,func_params,
			dim,coords1,coords2,crds_crx,1))
		    {
			screen("ERROR: in install_grid_crx(), no y-crxing!");
			clean_up(ERROR);
		    }
		    if (crds_crx[1] - coords1[1] < 0.004*h[1])
			crds_crx[1] = coords1[1] + 0.004*h[1];
		    if (coords2[1] - crds_crx[1] < 0.004*h[1])
			crds_crx[1] = coords2[1] - 0.004*h[1];
		    crx_store[n_crx].p = Point(crds_crx);
		    y_crx[i][j] = &crx_store[n_crx++];
		}
	    }
	}

	return n_crx;
}	/* end install_grid_crx */


LOCAL	void make_grid_curves(
	BLK_CRX_2D      *blk_crx,
	const EG_CRX_2D *eg_crx,
	int             *gmax,
	BOND		**first,
	int		*num_segs)
{
	int i,j,k,num_blk;
	BLK_BOND *bm,***blk_mem,*blk_mem_store;
	BOND *b;

	num_blk = 0;
	for (i = 0; i < gmax[0]; ++i)
	{
	    for (j = 0; j < gmax[1]; ++j)
	    {
		    if (onfront_block2d(i,j,eg_crx))
			++num_blk;
	    }
	}
	bi_array(&blk_mem,gmax[0],gmax[1],sizeof(BLK_BOND*));
        uni_array(&blk_mem_store,num_blk,sizeof(BLK_BOND));

	num_blk = 0;
	for (i = 0; i < gmax[0]; ++i)
	{
	    for (j = 0; j < gmax[1]; ++j)
	    {
		if (onfront_block2d(i,j,eg_crx))
		{
		    bm = blk_mem[i][j] = &blk_mem_store[num_blk++];
		    bm->blk_info = blk_crx->blk_info;
		    assign_blk_crx2d(blk_crx,i,j,eg_crx);
		    switch (blk_crx->blk_type)
		    {
		    case COMP2_BLOCK:
			construct_comp2_blk2d(blk_crx,bm);
			break;
		    default:
			screen("UNKNOWN BLOCK: code needed!\n");
			clean_up(ERROR);
		    }
		    if (i != 0 && blk_mem[i-1][j] != NULL)
			connect_adj_blk(blk_mem[i-1][j],bm);
		    if (j != 0 && blk_mem[i][j-1] != NULL)
			connect_adj_blk(blk_mem[i][j-1],bm);
		}
	    }
	}
	*num_segs = 0;
	for (i = 0; i < num_blk; ++i)
	{
	    int nc = blk_mem_store[i].num_curves;
	    for (j = 0; j < nc; ++j)
	    {
	    	int nb = blk_mem_store[i].num_bonds[j];
		for (k = 0; k < nb; ++k)
		{
		    if (blk_mem_store[i].bonds[j][k]->prev == NULL)
		    	first[(*num_segs)++] = blk_mem_store[i].bonds[j][k];
		}
	    }
	}

	free_these(2,blk_mem,blk_mem_store);
}	/* end make_grid_curves */


LOCAL	bool onfront_block2d(
	int          i,
	int          j,
	const EG_CRX_2D *eg_crx)
{
	int ii,jj;

	for (ii = 0; ii < 2; ++ii)
	{
	    for (jj = 0; jj < 2; ++jj)
	    {
		if (eg_crx->x_crx[i][j+ii] != NULL)
		    return YES;
		if (eg_crx->y_crx[i+jj][j] != NULL)
		    return YES;
	    }
	}
	return NO;
}	/* end onfront_block */


LOCAL	BLK_CRX_2D *alloc_blk_crx2d(
	bool alloc_BBI_POINT)
{
	BLK_CRX_2D *blk_crx;

	scalar(&blk_crx,sizeof(BLK_CRX));
	bi_array(&blk_crx->comp,2,2,sizeof(COMPONENT));
	bi_array(&blk_crx->ix,2,2,sizeof(int));
	bi_array(&blk_crx->iy,2,2,sizeof(int));
	bi_array(&blk_crx->crx,2,2,sizeof(BBI_POINT*));
	if (alloc_BBI_POINT)
	{
	    int i,j,k,num_crx = 0;
	    uni_array(&blk_crx->crx_store,5,sizeof(BBI_POINT));
	    for (i = 0; i < 2; ++i)
	    for (j = 0; j < 2; ++j)
	    	blk_crx->crx[i][j] = &blk_crx->crx_store[num_crx++];
	    blk_crx->node_crx = &blk_crx->crx_store[num_crx++];
	}
	return blk_crx;
}	/*end alloc_blk_crx*/


LOCAL void assign_blk_crx2d(
	BLK_CRX_2D   *blk_crx,
	int          i,
	int          j,
	const EG_CRX_2D *eg_crx)
{
	int ic,num_curve_crx,ii,jj;
	BBI_POINT ***x_crx = eg_crx->x_crx;
	BBI_POINT ***y_crx = eg_crx->y_crx;
	COMPONENT c,**comp = eg_crx->comp;

	blk_crx->num_comps = 0;
	for (ii = 0; ii < 2; ++ii)
	{
	    for (jj = 0; jj < 2; ++jj)
	    {
		c = comp[i+ii][j+jj];
		blk_crx->ix[ii][jj] = ii;
		blk_crx->iy[ii][jj] = jj;
		for (ic = 0; ic < blk_crx->num_comps; ++ic)
		{
		    if (c == blk_crx->comps[ic])
		    {
			++blk_crx->nv[ic];
			break;
		    }
		}
		if (ic == blk_crx->num_comps)
		{
		    blk_crx->comps[ic] = c;
		    blk_crx->nv[ic] = 1;
		    ++blk_crx->num_comps;
		}
		blk_crx->comp[ii][jj] = c;
	    }
	    blk_crx->crx[0][ii] = x_crx[i][j+ii];
	    blk_crx->crx[1][ii] = y_crx[i+ii][j];
	}
	for (ii = 0; ii < blk_crx->num_comps-1; ++ii)
	{
	    for (jj = 1; jj < blk_crx->num_comps; ++jj)
	    {
	    	if (blk_crx->nv[ii] > blk_crx->nv[jj])
		{
		    int nv_tmp;
		    COMPONENT c_tmp;
		    nv_tmp = blk_crx->nv[ii];
		    blk_crx->nv[ii] = blk_crx->nv[jj];
		    blk_crx->nv[jj] = nv_tmp;
		    c_tmp = blk_crx->comps[ii];
		    blk_crx->comps[ii] = blk_crx->comps[jj];
		    blk_crx->comps[jj] = c_tmp;
		}
	    }
	}
	blk_crx->blk_type = COMP2_BLOCK;
	return;
}	/* end assign_blk_crx */


LOCAL	int construct_comp2_blk2d(
	BLK_CRX_2D *blk_crx,
	BLK_BOND *blk_mem)
{
	int       i,j;
	COMPONENT **comp = blk_crx->comp;
	int       num_crx;
	CURVE 	  *c = blk_crx->blk_info->curves[0];
	COMPONENT left_c  = negative_component(c);
	COMPONENT right_c = positive_component(c);
	POINT	  *p1,*p2;

	/*
	printf("Entering construct_comp2_blk2d()\n");
	printf("nv = %d %d\n",blk_crx->nv[0],blk_crx->nv[1]);
	printf("comps = %d %d\n",blk_crx->comps[0],blk_crx->comps[1]);
	printf("comp[%d][%d] = %d\n",0,0,comp[0][0]);
	printf("comp[%d][%d] = %d\n",1,0,comp[1][0]);
	printf("comp[%d][%d] = %d\n",0,1,comp[0][1]);
	printf("comp[%d][%d] = %d\n",1,1,comp[1][1]);
	printf("crx[%d][%d] = %d\n",0,0,blk_crx->crx[0][0]);
	printf("crx[%d][%d] = %d\n",0,1,blk_crx->crx[0][1]);
	printf("crx[%d][%d] = %d\n",1,0,blk_crx->crx[1][0]);
	printf("crx[%d][%d] = %d\n",1,1,blk_crx->crx[1][1]);
	*/
	if (blk_crx->nv[0] == 1)
	{
	    blk_mem->num_bonds[0] = 1;
	    if (comp[0][0] == blk_crx->comps[0])
	    {
	    	p1 = blk_crx->crx[0][0]->p;
		p2 = blk_crx->crx[1][0]->p;
		blk_mem->bonds[0][0] = (comp[0][0] == left_c) ?
			Bond(p1,p2) : Bond(p2,p1);
	    }
	    else if (comp[0][1] == blk_crx->comps[0])
	    {
	    	p1 = blk_crx->crx[1][0]->p;
		p2 = blk_crx->crx[0][1]->p;
		blk_mem->bonds[0][0] = (comp[0][1] == left_c) ?
			Bond(p1,p2) : Bond(p2,p1);
	    }
	    else if (comp[1][0] == blk_crx->comps[0])
	    {
	    	p1 = blk_crx->crx[1][1]->p;
		p2 = blk_crx->crx[0][0]->p;
		blk_mem->bonds[0][0] = (comp[1][0] == left_c) ?
			Bond(p1,p2) : Bond(p2,p1);
	    }
	    else if (comp[1][1] == blk_crx->comps[0])
	    {
	    	p1 = blk_crx->crx[0][1]->p;
		p2 = blk_crx->crx[1][1]->p;
		blk_mem->bonds[0][0] = (comp[1][1] == left_c) ?
			Bond(p1,p2) : Bond(p2,p1);
	    }
	}
	else if (blk_crx->nv[0] == 2)
	{
	    if (comp[0][0] == comp[1][0])
	    {
	    	blk_mem->num_bonds[0] = 1;
	    	p1 = blk_crx->crx[1][1]->p;
	    	p2 = blk_crx->crx[1][0]->p;
		blk_mem->bonds[0][0] = (comp[0][0] == left_c) ?
			Bond(p1,p2) : Bond(p2,p1);
	    }
	    else if (comp[0][0] == comp[0][1])
	    {
	    	blk_mem->num_bonds[0] = 1;
	    	p1 = blk_crx->crx[0][0]->p;
	    	p2 = blk_crx->crx[0][1]->p;
		blk_mem->bonds[0][0] = (comp[0][0] == left_c) ?
			Bond(p1,p2) : Bond(p2,p1);
	    }
	    else if (comp[0][0] == comp[1][1])
	    {
	    	blk_mem->num_bonds[0] = 2;
	    	p1 = blk_crx->crx[0][0]->p;
		p2 = blk_crx->crx[1][0]->p;
		blk_mem->bonds[0][0] = (comp[0][0] == left_c) ?
			Bond(p1,p2) : Bond(p2,p1);
	    	p1 = blk_crx->crx[0][1]->p;
		p2 = blk_crx->crx[1][1]->p;
		blk_mem->bonds[0][1] = (comp[1][0] == left_c) ?
			Bond(p1,p2) : Bond(p2,p1);
	    }
	}
	blk_mem->num_curves = 1;
	blk_mem->curves[0] = c;
	for (i = 0; i < blk_mem->num_curves; ++i)
	for (j = 0; j < blk_mem->num_bonds[i]; ++j)
	{
	    blk_mem->bonds[i][j]->prev = blk_mem->bonds[i][j]->next = NULL;
	}

	return FUNCTION_SUCCEEDED;
}	/* end reconstruct_blk_intfc */


LOCAL	void connect_adj_blk(
	BLK_BOND *bm1,
	BLK_BOND *bm2)
{
	int i1,j1,nc1 = bm1->num_curves;
	int i2,j2,nc2 = bm2->num_curves;
	BOND *b1,*b2;
	CURVE *c;

	for (i1 = 0; i1 < nc1; ++i1)
	{
	    for (i2 = 0; i2 < nc2; ++i2)
	    {
	    	for (j1 = 0; j1 < bm1->num_bonds[i1]; ++j1)
		{
		    for (j2 = 0; j2 < bm2->num_bonds[i2]; ++j2)
		    {
		    	b1 = bm1->bonds[i1][j1];
		    	b2 = bm2->bonds[i2][j2];
			c = bm1->curves[i1];
			if (b1->end == b2->start)
			{
			    b1->next = b2;
			    b2->prev = b1;
			    c = bm1->curves[i1];
			    while (b1 != NULL)
			    {
			    	c->first = b1;
				c->start->posn = b1->start;
				if (b1 == b2) /* closed curve */
				{
				    c->last = b1->prev;
				    c->last->next = NULL;
				    c->first->prev = NULL;
				    c->end->posn = b1->start;
				    break;
				}
				b1 = b1->prev;
			    }
			    while (b2 != NULL && b1 != b2)
			    {
			    	c->last = b2;
				c->end->posn = b2->end;
				b2 = b2->next;
			    }
			}
			else if (b2->end == b1->start)
			{
			    b2->next = b1;
			    b1->prev = b2;
			    while (b2 != NULL)
			    {
			    	c->first = b2;
				c->start->posn = b2->start;
				if (b1 == b2) /* closed curve */
				{
				    c->last = b2->prev;
				    c->last->next = NULL;
				    c->first->prev = NULL;
				    c->end->posn = b2->start;
				    break;
				}
				b2 = b2->prev;
			    }
			    while (b1 != NULL && b1 != b2)
			    {
			    	c->last = b1;
				c->end->posn = b1->end;
				b1 = b1->next;
			    }
			}
		    }
		}
	    }
	}
}	/* end connect_adj_blk */

EXPORT 	CURVE *prompt_make_linear_curve(
	INTERFACE *intfc,
	RECT_GRID *gr)
{
	INTERFACE       *infc;
	float coords1[2],coords2[2],coords[2];
	POINT *p1,*p2;
	NODE *ns,*ne;
	CURVE *curve;
	COMPONENT left_c,right_c;
	int i,num_pts;
	float dx,dy;

	screen("Enter two integers as the left and right components: ");
	Scanf("%d %d\n",&left_c,&right_c);
	screen("Enter the coordinates of the start node: ");
	Scanf("%f %f\n",&coords1[0],&coords1[1]);
	screen("Enter the coordinates of the end node: ");
	Scanf("%f %f\n",&coords2[0],&coords2[1]);

	p1 = Point(coords1);
	ns = make_node(p1);
	p2 = Point(coords2);
	ne = make_node(p2);
	curve = make_curve(left_c,right_c,ns,ne);

	screen("Enter number of interior points of the line: ");
	Scanf("%d\n",&num_pts);
	dx = (coords2[0] - coords1[0])/(num_pts+1);
	dy = (coords2[1] - coords1[1])/(num_pts+1);
	for (i = 0; i < num_pts; ++i)
	{
	    coords[0] = coords1[0] + (i+1)*dx;
	    coords[1] = coords1[1] + (i+1)*dy;
	    if (insert_point_in_bond(Point(coords),curve->last,curve) !=
	                    FUNCTION_SUCCEEDED)
	    {
	    	screen("ERROR in prompt_make_linear_curve(), "
		       "insert_point_in_bond() failed\n");
		clean_up(ERROR);
	    }
	}
	return curve;
}	/* end prompt_make_linear_curve */

EXPORT 	CURVE *prompt_make_elliptic_curve(
	INTERFACE *intfc,
	RECT_GRID *gr)
{
	COMPONENT compin,compout;
	ELLIP_PARAMS ellip;
	float cen[MAXD],rad[MAXD];
	CURVE *curve;

	screen("Enter two integers as the in and out components: ");
	Scanf("%d %d\n",&compin,&compout);
	screen("Enter the coordinates of the elliptic center: ");
	Scanf("%f %f\n",&cen[0],&cen[1]);
	screen("Enter the radii of the ellipse: ");
	Scanf("%f %f\n",&rad[0],&rad[1]);

	ellip.cen = cen;	ellip.rad = rad;
	ellip.closed = YES;	
	ellip.nor_orient = POSITIVE_ORIENTATION;
	ellip.dim = 2;		
	ellip.gr = &topological_grid(intfc);
	ellip.fpoly = NULL;
	ellip.lpoly = NULL;
	ellip.ThetaS[0] = 0.0;
	curve = make_elliptic_curve(&ellip,compin,compout,0.5);
	return curve;

}	/* end prompt_make_elliptic_curve */

EXPORT void prompt_make_level_curves(
	INTERFACE *intfc,
	RECT_GRID *gr)
{
	COMPONENT left_c,right_c;
	POINTER func_params;
	float (*func)(POINTER,float*);
	char s[20];

	screen("Enter two integers as the left and right components: ");
	Scanf("%d %d\n",&left_c,&right_c);

	screen("Enter the level curve options\n");
	screen("Supported level curve types are \n"
	       "\tStraight line (l)\n"
	       "\tEllipse (e)\n"
	       "\tMultiple circles (m)\n"
	       "\tSine waves (s)\n"
	       "\tTest disk (t)\n"
	       "Enter choice: ");
	(void) Scanf("%s\n",s);
	switch (s[0])
	{
	case 'l':
	case 'L':
	    func_params = init_line_params(gr);
	    func = line_func;
	    break;
	case 'e':
	case 'E':
	    func_params = init_ellipse_params(gr);
	    func = ellipse_func;
	    break;
	case 'm':
	case 'M':
	    func_params = init_multi_circle_params(gr);
	    func = multi_circle_func;
	    break;
	case 's':
	case 'S':
	    func_params = init_sine_wave_params(gr);
	    func = multi_sine_mode_func;
	    break;
	case 't':
	case 'T':
	    func_params = init_test_disk_params(gr);
	    func = test_disk_func;
	    break;
	}
	if (!make_level_curve(gr,intfc,left_c,right_c,func,func_params,NO))
        {
            screen("make_level_curve() failed!\n");
	    clean_up(ERROR);
        }
}	/* end prompt_make_level_interface */


/*	Initialization of curve parameters
*/

LOCAL	POINTER init_line_params(RECT_GRID *gr)
{
	static LINE_PARAMS lparams;
	screen("The equation of line is a*x + b*y = c\n");
	screen("Enter a, b, and c: ");
	Scanf("%f %f %f\n",&lparams.a,&lparams.b,&lparams.c);
	return (POINTER)&lparams;
}	/* end init_line_params */

EXPORT	POINTER init_ellipse_params(RECT_GRID *gr)
{
	static ELLIP2D_PARAMS eparams;
	screen("The equation of ellipse is "
	    	   "(x-x0)^2/a^2 + (y-y0)^2/b^2 = 1\n");
	screen("Enter center x0, and y0: ");
	Scanf("%f %f\n",&eparams.x0,&eparams.y0);
	screen("Enter radii a, and b: ");
	Scanf("%f %f\n",&eparams.a,&eparams.b);
	return (POINTER)&eparams;
}	/* end init_ellipse_params */

LOCAL	POINTER init_multi_circle_params(RECT_GRID *gr)
{
	static MC_PARAMS mc_params;
	int i;
	screen("Enter number of circles: ");
	Scanf("%d\n",&mc_params.num_cir);
	uni_array(&mc_params.rad,mc_params.num_cir,FLOAT);
	bi_array(&mc_params.cen,mc_params.num_cir,2,FLOAT);
	for (i = 0; i < mc_params.num_cir; ++i)
	{
	    screen("Enter the center coords of the %d-th circle: ",i+1);
	    Scanf("%f %f\n",&mc_params.cen[i][0],&mc_params.cen[i][1]);
	    screen("Enter the radius of the %d-th circle: ",i+1);
	    Scanf("%f\n",&mc_params.rad[i]);
	}
	return (POINTER)&mc_params;
}	/* end init_multi_circle_params */

LOCAL	POINTER init_sine_wave_params(RECT_GRID *gr)
{
	static FOURIER_POLY sine_params;
	int i;
	screen("Enter mean position of the curve (in z-direction): ");
	Scanf("%f\n",&sine_params.z0);
	screen("Enter number of modes: ");
	Scanf("%d\n",&sine_params.num_modes);
	bi_array(&sine_params.nu,sine_params.num_modes,1,FLOAT);
	uni_array(&sine_params.A,sine_params.num_modes,FLOAT);
	uni_array(&sine_params.phase,sine_params.num_modes,FLOAT);
	for (i = 0; i < sine_params.num_modes; ++i)
	{
	    screen("Enter the frequency of the %d-th mode: ",i);
	    Scanf("%f\n",&sine_params.nu[i][0]);
	    screen("Enter the amplitude of the %d-th mode: ",i);
	    Scanf("%f\n",&sine_params.A[i]);
	    screen("Enter the phase of the %d-th mode: ",i);
	    Scanf("%f\n",&sine_params.phase[i]);
	    sine_params.phase[i] *= PI/180.0;
	}
	sine_params.dim = 2;
	sine_params.U = gr->U;
	sine_params.L = gr->L;
	return (POINTER)&sine_params;
}	/* end init_sine_wave_params */


LOCAL	POINTER init_test_disk_params(RECT_GRID *gr)
{
	static TDISK_PARAMS td_params;
	screen("Enter center of the disk x0, and y0: ");
	Scanf("%f %f\n",&td_params.x0,&td_params.y0);
	screen("Enter radius of the disk: ");
	Scanf("%f\n",&td_params.r);
	screen("Concave well faces up (along the y-axis)\n");
	screen("Enter the width and depth of the well: ");
	Scanf("%f %f\n",&td_params.w,&td_params.h);
	return (POINTER)&td_params;
}	/* end init_test_disk_params */

/*	Level functions for curves 
*/

LOCAL float line_func(
	POINTER func_params, 
        float *coords)
{
	LINE_PARAMS *lparams = (LINE_PARAMS*)func_params;
	float dist;
	float a = lparams->a;
	float b = lparams->b;
	float c = lparams->c;

	dist = (a*coords[0] + b*coords[1] - c)/sqrt(a*a + b*b);
	return dist;
}	/* end line_func */

LOCAL float multi_circle_func(
	POINTER func_params, 
        float *coords)
{
	MC_PARAMS *mc_params = (MC_PARAMS*)func_params;
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
}	/* end multi_circle_func */

EXPORT float ellipse_func(
	POINTER func_params, 
        float *coords)
{
	ELLIP2D_PARAMS *eparams = (ELLIP2D_PARAMS*)func_params;
	float dist;

	dist = sqr(coords[0] - eparams->x0)/sqr(eparams->a) +
	       sqr(coords[1] - eparams->y0)/sqr(eparams->b) - 1.0;
	return dist;
}	/* end ellipse_func */

LOCAL float test_disk_func(
	POINTER func_params, 
        float *coords)
{
	TDISK_PARAMS *td_params = (TDISK_PARAMS*)func_params;
	float xl,xr,yh;
	float dist1,dist2,dist3;

	xl = td_params->x0 - td_params->w/2.0;
	xr = td_params->x0 + td_params->w/2.0;
	yh = td_params->y0 + td_params->r - td_params->h;
	dist1 = sqrt(sqr(coords[0] - td_params->x0) + 
                sqr(coords[1] - td_params->y0)) - td_params->r;
	dist2 = (coords[0] < td_params->x0) ? coords[0] - xl :
		xr - coords[0];
	dist3 = coords[1] - yh;
	if (dist1 > 0.0) return dist1;
	else if (dist2 > 0.0 && dist3 > 0.0)
	    return min(dist2,dist3);
	else if (dist2 > 0 && dist3 < 0.0)
	    return max(dist1,dist3);
	else if (dist2 < 0 && dist3 > 0.0)
	    return max(dist1,dist2);
	else
	    return max(dist1,max(dist2,dist3));
}	/* end test_disk_func */


EXPORT CURVE *make_array_curve(
	INTERFACE   *intfc,
	COMPONENT   left_comp,
	COMPONENT   right_comp,
	int	    num_points,
	float	    **point_array,
	bool	    is_closed_curve)
{
	int i;
	CURVE *c;
	NODE *ns,*ne;
	INTERFACE *save_intfc;
	float coords[MAXD];

	save_intfc = current_interface();
        set_current_interface(intfc);

	if (is_closed_curve)
	{
	    coords[0] = point_array[0][0];
	    coords[1] = point_array[0][1];
	    ns = ne = make_node(Point(coords));
	}
	else
	{
	    coords[0] = point_array[0][0];
	    coords[1] = point_array[0][1];
	    ns = make_node(Point(coords));
	    coords[0] = point_array[num_points-1][0];
	    coords[1] = point_array[num_points-1][1];
	    ne = make_node(Point(coords));
	    num_points--;
	}
	c = make_curve(left_comp,right_comp,ns,ne);
	for (i = 1; i < num_points; ++i)
	{
	    coords[0] = point_array[i][0];
	    coords[1] = point_array[i][1];
	    insert_point_in_bond(Point(coords),c->last,c);
	}

	set_current_interface(save_intfc);
	return c;
}	/* end make_array_curve */
