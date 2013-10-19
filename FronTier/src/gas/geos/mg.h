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
*				mg.h:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Template file for the implementation of new equation of state models
*
*/

#if !defined(_MG_H)
#define _MG_H

#include <geos/geosdecs.h>

struct _Tref_Spline {
    float *rho; /* Independent variable array, monotone nondecreasing */
    float *T;   /* Dependent variable array */
    int    m;    /* Number of data points */

    int    k;    /* Degree of spline polynomial */
    int    n;    /* Number of knots in the spline approximation */
    float *t;   /* Knot array */
    float *c;   /* B-spline coefficients */
    int   nest; /* Allocated length of the arrays t and c */

    float     Tref_atol;
    float     Tref_rtol;
    float     Tref_terr;
};
typedef struct _Tref_Spline Tref_Spline;

struct _MG_params {
    const char   *name;	/*Name of material,  may be NULL*/
    float        Rho0;	/* Reference state density*/
    float        P0;	/* Reference state pressure*/
    float        E0;	/* Reference state energy*/
    float        T0;	/* Reference state temperature*/
    float        C0;	/* Reference state sound speed*/
    float        S1,    /* U_s = C0 + (S1 + S2*U_s/U_p + S3*(U_s/U_p)^2)*U_p */
                 S2,    /* along reference Hugoniot */
	         S3;     
    float        gamma0;/* Gruneisen exponent at reference state */
    float        b;	/* Gruneisent exponent at infinite compression */
    float        CV;	/* Specific heat at constant volume */

    float        Rho_max, /* V_min = 1/Rho_max,  Rho_max is the maximum  */
	         V_min;	  /* allowed density,  for rho > Rho_max either  */
			  /* the sound speed on the reference curve      */
			  /* is imaginary or the reference curve         */
			  /* blows up somewhere between Rho_max and rho  */
    float       P_inf;
    float       E_inf;
    Tref_Spline *TrefSpline;
};
typedef struct _MG_params MG_params;

	/* Structures for accessing equation of state information */

struct _MG_EOS {
	EOS	eos;

	MG_params MGparams;
};
typedef struct _MG_EOS MG_EOS;

#define	MG_Eos(state)	((MG_EOS *)Params(state)->eos)

#endif /* !defined(_MG_H) */
