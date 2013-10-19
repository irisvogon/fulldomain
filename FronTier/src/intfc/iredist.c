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
*				iredist.c:
*
*       Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Containing function of redistribute interface.
*
*/

#if defined(TWOD)

#include <intfc/int.h>

#define DEBUG_STRING "i_redistribute"

LOCAL	int rb_cell_index(float*,int,RECT_GRID*);
LOCAL   void  replace_curve_seg_by_bond(CURVE*,BOND*,BOND*);

/*
*			equi_redist_curve_seg():
*
*	Redistributes a curve segment from bs->start to be->end by
*	deleting all points between, and inserting N equally spaced points.
*	NOTE that the resultant point spacing is equal only if measured
*	with respect to the arc length of the ORIGINAL, input curve
*	segment and NOT with respect to the final redistributed curve
*	segment.
*	If nbds <= 0, the routine computes N such that the spacing between
*	points on the redistributed curve is <= space.
*	If nbds > 0, the routine uses N == nbds.
*/

EXPORT	void equi_redist_curve_seg(
	CURVE		*c,
	BOND		*bs,
	BOND		*be,
	int		nbds,
	float		seg_len,
	float		space,
	RECT_GRID	*rgr)
{
	BOND		*b, *bstart, *bend;
	float		b_len, sc_len, offset, s, oms;
	float		coords[MAXD];
	float		*h = rgr->h;
	int		dim = rgr->dim;
	int		new_nbds;

	DEBUG_ENTER(equi_redist_curve_seg)
	if (nbds <= 0)
	{
	    new_nbds = (int)(ceil((double)seg_len/(double)space));
	    b_len = seg_len/(float)new_nbds;
	}
	else
	{
	    new_nbds = nbds;
	    b_len = seg_len/(float)new_nbds;
	}
	if (new_nbds <= 1)
	{
	    if (is_closed_curve(c) && bs == c->first && be == c->last)
	    {
	    	new_nbds = c->num_points - 1;
	    	if (new_nbds > 1)
	            equi_redist_curve_seg(c,bs,be,new_nbds,
					  seg_len,space,rgr);
	    }
	    else
	        replace_curve_seg_by_bond(c,bs,be);

	    DEBUG_LEAVE(equi_redist_curve_seg)
	    return;
	}

	offset = b_len;
	bstart = bs;		bend = be->next;
	while (bstart != bend)
	{
	    b = bstart;
	    while ((sc_len = scaled_bond_length(b,h,dim)) < offset)
	    {
	    	if (b->next == bend)
	    	{
	    	    replace_curve_seg_by_bond(c,bstart,b);
	    	    goto leave;
	    	}
	    	offset -= sc_len;
	    	b = b->next;
	    }
	    if ((b->next != bend) ||
	        (sc_len >= offset + MIN_SC_SEP(c->interface)))
	    {
	    	s = offset/sc_len;	oms = 1.0 - s;
	    	coords[0] = oms * Coords(b->start)[0] + s * Coords(b->end)[0];
	    	coords[1] = oms * Coords(b->start)[1] + s * Coords(b->end)[1];
	    	if (insert_point_in_bond(Point(coords),b,c) !=
		    FUNCTION_SUCCEEDED)
	    	{
	    	    screen("ERROR in equi_redist_curve_seg(), "
	    	           "insert_point_in_bond failed\n");
	    	    clean_up(ERROR);
	    	}
	    }
	    replace_curve_seg_by_bond(c,bstart,b);
	    bstart = bstart->next;
	    offset = b_len;
	}

leave:
	DEBUG_LEAVE(equi_redist_curve_seg)
	return;
}		/*end equi_redist_curve_seg*/


EXPORT	void	rect_bdry_curve_redist(
	CURVE		*c,
	ORIENTATION	orient,
	RECT_GRID	*gr,
	float		*tol)
{
	BOND		*b, *bb;
	float		tolx = tol[0], toly = tol[1];
	float		*ps, *pe;
	float		x, y, coords[MAXD];
	int		ixs, ixe, iys, iye, ix, iy;
	int		ibegin, ifinish;
	ORIENTATION	opor = Opposite_orient(orient);
	int		di;

	DEBUG_ENTER(rect_bdry_curve_redist)
	ps = Coords(c->start->posn);		pe = Coords(c->end->posn);

	if (fabs(ps[1]-pe[1]) < fabs(ps[0]-pe[0]))
	{

	    /* horizontal boundary */

	    ixs = rb_cell_index(ps,0,gr);
	    ixe = rb_cell_index(pe,0,gr);

	    if (ixe == ixs)
	    {
	    	while (c->first != c->last)
		    i_delete_point_adjacent_to_node(c,POSITIVE_ORIENTATION);
	    	DEBUG_LEAVE(rect_bdry_curve_redist)
	    	return;
	    }
	    if (orient == POSITIVE_ORIENTATION)
	    {
	    	di = (ps[0] < pe[0]) ? 1 : -1;
	    	ibegin = ixs;
	    	ifinish = ixe;
	    }
	    else
	    {
	    	di = (ps[0] < pe[0]) ? -1 : 1;
	    	ibegin = ixe;
	    	ifinish = ixs;
	    }

			/* Check positions at nodes */
		
	    x = cell_center(ibegin,0,gr);
	    if ((fabs(x - Coords(Node_of(c,orient)->posn)[0]) < tolx) ||
	    				(!Between(x,ps[0],pe[0])))
	    	ibegin += di;

	    x = cell_center(ifinish,0,gr);
	    if ((fabs(x - Coords(Node_of(c,opor)->posn)[0]) < tolx) ||
	    				(!Between(x,ps[0],pe[0])))
	    	ifinish -= di;

	    if (di*(ibegin - ifinish) > 0)
	    {
	    	while (c->first != c->last)
		    i_delete_point_adjacent_to_node(c,POSITIVE_ORIENTATION);
	    	DEBUG_LEAVE(rect_bdry_curve_redist)
	    	return;
	    }

	    coords[1] = y = Coords(c->start->posn)[1];
	    b = Bond_at_node(c,orient);
	    for (ix = ibegin; di*(ifinish - ix) >= 0; ix += di)
	    {
	        x = cell_center(ix,0,gr);
	        for (bb = b; bb != NULL; bb = Following_bond(bb,orient))
	        {
	    	    if (fabs(x-Coords(Point_of_bond(bb,opor))[0]) < tolx)
	    	    {
	    		Coords(Point_of_bond(bb,opor))[0] = x;
	    		break;
	    	    }
	    	    else if (Between(x,Coords(bb->start)[0],Coords(bb->end)[0]))
	    	    {
	    	        coords[0] = x;
	    	        if (insert_point_in_bond(Point(coords),bb,c) !=
			    FUNCTION_SUCCEEDED)
	    	        {
	    	            screen("ERROR in rect_bdry_curve_redist(), "
	    	                   "insert_point_in_bond failed\n");
	    	            clean_up(ERROR);
	    	        }
	    	        if (orient == NEGATIVE_ORIENTATION) 
	    	        {
	    		    if (bb == b)
				b = bb->next;
	    		    bb = bb->next;
	    	        }
	    	        break;
	            }
	        }
		if (bb == NULL)
		{
		    bb = Bond_at_node(c,opor);
		    if (bb == b) break;
		}
		if (orient == POSITIVE_ORIENTATION)
		{
		    replace_curve_seg_by_bond(c,b,bb);
		    b = b->next;
		}
		else
		{
		    replace_curve_seg_by_bond(c,bb,b);
		    b = bb->prev;
		}
		if (b == NULL) break;
	    }
	    if (orient == POSITIVE_ORIENTATION)
	    	replace_curve_seg_by_bond(c,b,Bond_at_node(c,opor));
	    else
	    	replace_curve_seg_by_bond(c,Bond_at_node(c,opor),b);
	}
	else if (fabs(ps[0]-pe[0]) < fabs(ps[1]-pe[1]))
	{	
	    /* vertical boundary */

	    iys = rb_cell_index(ps,1,gr);
	    iye = rb_cell_index(pe,1,gr);
	    if (iye == iys)
	    {
	    	while (c->first != c->last)
		    i_delete_point_adjacent_to_node(c,POSITIVE_ORIENTATION);
	    	DEBUG_LEAVE(rect_bdry_curve_redist)
	    	return;
	    }
	    if (orient == POSITIVE_ORIENTATION)
	    {
	    	di = (ps[1] < pe[1]) ? 1 : -1;
	    	ibegin = iys;
	    	ifinish = iye;
	    }
	    else
	    {
	    	di = (ps[1] < pe[1]) ? -1 : 1;
	    	ibegin = iye;
	    	ifinish = iys;
	    }

	    /* Check positions at nodes */

	    y = cell_center(ibegin,1,gr);
	    if ((fabs(y - Coords(Node_of(c,orient)->posn)[1]) < toly) ||
	    				(!Between(y,ps[1],pe[1])))
	    	ibegin += di;

	    y = cell_center(ifinish,1,gr);
	    if ((fabs(y - Coords(Node_of(c,opor)->posn)[1]) < toly) ||
	    				(!Between(y,ps[1],pe[1])))
	        ifinish -= di;

	    if (di*(ibegin - ifinish) > 0)
	    {
	    	while (c->first != c->last)
		    i_delete_point_adjacent_to_node(c,POSITIVE_ORIENTATION);
	    	DEBUG_LEAVE(rect_bdry_curve_redist)
	    	return;
	    }

	    coords[0] = x = Coords(c->start->posn)[0];
	    b = Bond_at_node(c,orient);
	    for (iy = ibegin; di*(ifinish - iy) >= 0; iy += di)
	    {
	        y = cell_center(iy,1,gr);
	        for (bb = b; bb != NULL; bb = Following_bond(bb,orient))
	        {
	    	    if (fabs(y-Coords(Point_of_bond(bb,opor))[1]) < toly)
	    	    {
	    		Coords(Point_of_bond(bb,opor))[1] = y;
	    		break;
	    	    }
	    	    else if (Between(y,Coords(bb->start)[1],Coords(bb->end)[1]))
	    	    {
	    	        coords[1] = y;

	    	        if (insert_point_in_bond(Point(coords),bb,c) !=
			    FUNCTION_SUCCEEDED)
	    	        {
	    	            screen("ERROR in rect_bdry_curve_redist(), "
	    	                   "insert_point_in_bond failed\n");
	    	            clean_up(ERROR);
	    	        }
	    	        if (orient == NEGATIVE_ORIENTATION)
	    	        {
	    	            if (bb == b) b = bb->next;
	    	            bb = bb->next;
	    	        }
	    	        break;
	    	    }
	        }  
		if (bb == NULL)
		{
		    bb = Bond_at_node(c,opor);
		    if (bb == b) break;
		}
		if (orient == POSITIVE_ORIENTATION)
		{
		    replace_curve_seg_by_bond(c,b,bb);
		    b = b->next;
		}
		else
		{
		    replace_curve_seg_by_bond(c,bb,b);
		    b = bb->prev;
		}
		if (b == NULL)
		    break;
	    }
	    if (orient == POSITIVE_ORIENTATION)
	    	replace_curve_seg_by_bond(c,b,Bond_at_node(c,opor));
	    else
	    	replace_curve_seg_by_bond(c,Bond_at_node(c,opor),b);
	}
	DEBUG_LEAVE(rect_bdry_curve_redist)
}		/*end rect_bdry_curve_redist*/


EXPORT bool i_delete_point_adjacent_to_node(
	CURVE           *c,
	ORIENTATION     orient)
{
	BOND		*b = Bond_at_node(c,orient);
	bool		stat;

	if (orient == POSITIVE_ORIENTATION) 
	{
		if (b->next == NULL) return FUNCTION_FAILED;
		if (debugging("delete_point_adjacent_to_node"))
		{
			(void) printf("deleting start of bond - ");
			print_bond(b->next);
		}

		stat = delete_start_of_bond(b->next,c);
	}
	else 
	{
		if (b->prev == NULL) return FUNCTION_FAILED;
		if (debugging("delete_point_adjacent_to_node"))
		{
			(void) printf("deleting start of bond - ");
			print_bond(b);
		}

		stat = delete_start_of_bond(b,c);
	}

	return stat;
}	/* end i_delete_point_adjacent_to_node */


LOCAL	int rb_cell_index(
	float		*p,
	int		i,
	RECT_GRID	*gr)
{
	int		indx;
	int		*gmax = gr->gmax;
	int		*lbuf = gr->lbuf, *ubuf = gr->ubuf;

	indx = cell_index(p[i],i,gr);
	if (indx < -lbuf[i])
	    indx = -lbuf[i];
	if (indx >= gmax[i] + ubuf[i])
	    indx = gmax[i] + ubuf[i] - 1;
	return indx;
}		/*end rb_cell_index*/

/*
*			replace_curve_seg_by_bond():
*
*	Replaces the curve segment from the start of bond bs to the end
*	of bond be by a single bond. The resultant single bond has the
*	same pointer value as bs.
*/

LOCAL   void replace_curve_seg_by_bond(
	CURVE		*c,
	BOND		*bs,
	BOND		*be)
{
	BOND		*bnext;

	if ((bs == NULL) || (be == NULL))
	    return;
	if (bs == be)
	    return;
	bnext = be->next;
	while (bs->next != bnext)
		(void) delete_start_of_bond(bs->next,c);
}		/*end replace_curve_seg_by_bond*/

EXPORT  void rect_bdry_redist2d(
	INTERFACE	*intfc,
	RECT_GRID	*rgr,
	int		iperm)
{
	CURVE		**c;
	float		tol[MAXD];
	int		i;
	ORIENTATION	orient;

	DEBUG_ENTER(rect_bdry_redist2d)
	for (i = 0; i < rgr->dim; ++i)
	    tol[i] = MIN_SC_SEP(intfc) * rgr->h[i];
	orient = (iperm % 2) ? POSITIVE_ORIENTATION : NEGATIVE_ORIENTATION;
	for (c = intfc->curves; c && *c; ++c) 
	{
	    if (!is_bdry(*c))
	        continue;

	    rect_bdry_curve_redist(*c,orient,rgr,tol);
	}
	DEBUG_LEAVE(rect_bdry_redist2d)
}		/*end rect_bdry_redist2d*/


/*
*			closed_curve_node_redistribute():
*
*	A node on a closed curve never gets redistributed. With the
*	assumption that a node on a closed curve is more of a convenience
*	than an important point, this routine randomly changes the
*	position of the node on a closed curve each time redistribute
*	is called.
*/

EXPORT	bool closed_curve_node_redistribute(
	INTERFACE	*intfc,
	bool		status)
{
	BOND		*b;
	CURVE		**c;

#if defined(DEBUG_NODE_PROPAGATE)
        if (debugging("noCCNR"))
	    return status;
#endif /* defined(DEBUG_NODE_PROPAGATE) */
	for (c = intfc->curves; *c ; ++c)
	{
	    if (is_bdry(*c))
		continue;
	    if (!is_closed_curve(*c))
		continue;

			/* pick a random bond - any bond */

	    b = random_bond_on_curve(*c);
	    if (!move_closed_loop_node(*c,b))
		status = NO;
	}
	return status;
}		/*end closed_curve_node_redistribute*/
#endif /* defined(TWOD) */
