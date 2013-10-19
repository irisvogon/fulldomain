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
*
*			          ideal.h	
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Template file for the implementation of new equation of state models
*
*/

#if !defined(_IDEAL_H)
#define _IDEAL_H

#include <geos/geosdecs.h>

	/* Structures for accessing equation of state information */

struct _IDEAL_EOS_Component
{
	char	_name[10];			/* Name of Component */
	double	_coeff_small[7]; 		/* Coeffs for T<T0 */
	double	_coeff_large[7];		/* Coeffs for T>=T0 */
};
typedef struct _IDEAL_EOS_Component IDEAL_EOS_Component;

#define name(component)    			((IDEAL_EOS_Component*)component)->_name
#define M(component)    			((IDEAL_EOS_Component*)component)->_M
#define coeff_small(component,n)    		((IDEAL_EOS_Component*)component)->_coeff_small[n]
#define coeff_large(component,n)    		((IDEAL_EOS_Component*)component)->_coeff_large[n]

struct _IDEAL_EOS {
	EOS	eos;
	
	/* Constants */
	float  R;	/* Gas-Constant */
	float  T;	/* Critical Temperature */

	float _M[MAX_NUM_GAS_COMPS];	/* Molecular weights of components */

	/* Properties from NASA File, that will be used in this simulation */
	IDEAL_EOS_Component* comp;	// Size is params->n_comps
};
typedef struct _IDEAL_EOS IDEAL_EOS;

#define	IDEAL_Eos(state)			((IDEAL_EOS*)Params(state)->eos)
#define	IDEAL_Eos_Component(state,n)		(&IDEAL_Eos(state)->comp[n])
#define	IDEAL_Eos_Component_Count(state)	(IDEAL_Eos(state)->nComps)
#define	IDEAL_Eos_Component_Alloc(state,N)	(IDEAL_Eos(state)->comp=(IDEAL_EOS_Component*)malloc(N*sizeof(IDEAL_EOS_Component)))

typedef struct {
	float *c_small;
	float *c_large;
	float E;
	float R;
} IDEAL_RT_PRMS;

#endif /* !defined(_IDEAL_H) */
