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
*				fvelo.h
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*
*	Contains structures related to velocity functions
*/
#if !defined(_FVELO_H)
#define _FVELO_H

struct _TRANS_PARAMS
{
	int dim;
	float vel[MAXD];
};
typedef struct _TRANS_PARAMS TRANS_PARAMS;

struct _RADIAL_MOTION_PARAMS
{
	int dim;
	float cen[MAXD];
	float speed;
};
typedef struct _RADIAL_MOTION_PARAMS RADIAL_MOTION_PARAMS;

struct _SHEAR_MOTION_PARAMS
{
	int dim;
	float dir[MAXD];
	float coords[MAXD];
	float dvdh;
};
typedef struct _SHEAR_MOTION_PARAMS SHEAR_MOTION_PARAMS;

struct _SINE_MOTION_PARAMS
{
	int dim;
	float wave_length;
	float vmax;
	float phase;
};
typedef struct _SINE_MOTION_PARAMS SINE_MOTION_PARAMS;

struct _CIRCULAR_ROTATION_PARAMS
{
	int dim;
	float x_0,y_0;
	float omega_0;
	float grad;
};
typedef struct _CIRCULAR_ROTATION_PARAMS CIRCULAR_ROTATION_PARAMS;

struct _NORV_PARAMS
{
	int dim;
	float coeff;
	float epsilon;
};
typedef struct _NORV_PARAMS NORV_PARAMS;

struct _FLAME_PARAMS
{
        POINTER stuff;
        int problem_type;
        int dim;
        float alpha;
        float U;
        float m;
        float c;
        float epsilon;
        float wind[MAXD];
        float a,b;
};
typedef struct _FLAME_PARAMS FLAME_PARAMS;

struct _BURGERS_PARAMS
{
	int dim;
	float cen[MAXD];
	float speed;
};
typedef struct _BURGERS_PARAMS BURGERS_PARAMS;

struct _BIPOLAR_PARAMS
{
	int dim;
	float i1,i2;
	float cen1[MAXD],cen2[MAXD];
	float reverse_time;
};
typedef struct _BIPOLAR_PARAMS BIPOLAR_PARAMS;

struct _VORTEX_PARAMS
{
    int dim;
    float coeff;
    float time;
    float cos_time;
    char type[10];
    float cen[MAXD];
    float rad;
};
typedef struct _VORTEX_PARAMS VORTEX_PARAMS;

#endif	/* !defined(_FVELO_H) */
