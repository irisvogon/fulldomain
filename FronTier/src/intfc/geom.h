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
*			geom.h
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#if !defined(_IGEOM_H)
#define _IGEOM_H

#include <cdecs.h>
#include <vmalloc.h>
#include <intfc/triangledefs.h>

/* The default is to compile all dimensional versions of the code */
#if !defined(ONED) && !defined(TWOD) && !defined(THREED)
#   define ONED
#   define TWOD
#   define THREED
#endif /* !defined(ONED) && !defined(TWOD) && !defined(THREED) */

enum {
#if defined(THREED)
	MAXD = 3
#elif defined(TWOD)
	MAXD = 2
#elif defined(ONED)
	MAXD = 1
#endif /* defined(THREED) */
};

enum {
#if defined(THREED)
	MAX_LSQ_PTS = 100
#elif defined(TWOD)
	MAX_LSQ_PTS = 40
#elif defined(ONED)
	MAX_LSQ_PTS = 5
#endif /* defined(THREED) */
};

#define	NOBUF	NULL

	/* Possible values for the variable remap */
typedef	enum {
	INVALID_REMAP     = 0,
	IDENTITY_REMAP    = 1,	/* identity conformal map */
	CYLINDRICAL_REMAP = 2,	/* cylindrical geometry */
	SPHERICAL_REMAP   = 3	/* spherical geometry (oned) */
} GEOMETRY_REMAP;

typedef	enum {
	COMP_GRID = 1,
	DUAL_GRID = 2,
	EXPANDED_DUAL_GRID = 3
} GRID_TYPE;

#define d_index1d(ix,gmax)       \
        (ix)

#define d_index2d(ix,iy,gmax)                                        \
        ((iy)*((gmax)[0]+1)+(ix))

#define d_index3d(ix,iy,iz,gmax)                                        \
        (((iz)*((gmax)[1]+1)+(iy))*((gmax)[0]+1)+(ix))


		/* Structure Defining a Regular Rectangular Grid: */
/*
*	Support has been added to support virtual domain buffer zones
*	as well as possibly nonconstant grid spacings.  The storage
*	for the grid lines and spacings is taken from the single
*	array glstore.  The allocated size of glstore should be
*
*				dim-1
*				-----
*				\
*				 \
*	size_of(glstore) =        \	(3*(gmax[i]+ub[i]+lb[i]) + 1)
*				  /
*				 /
*				/
*				-----
*				i = 0
*
*	Which is just the combined length of the edge and dl uni_arrays.
*	This storage should be mapped out as follows
*	
*	edge[0] = glstore + lb[0]
*	edge[i]   = edge[i-1]   + gmax[i-1]   + ub[i-1]   + 1 + lb[i]	(i > 0)
*	dl[0]     = edge[dim-1] + gmax[dim-1] + ub[dim-1] + 1 + lb[0]
*	dl[i]     = dl[i-1]     + gmax[i-1]   + ub[i-1]       + lb[i]	(i > 0)
*	center[0] = dl[dim-1]   + gmax[dim-1] + ub[dim-1]     + lb[0]
*	center[i] = center[i-1] + gmax[i-1]   + ub[i-1]       + lb[i]	(i > 0)
*
*	This arrangement has the following benefits,
*
*	1. Storage need only be allocated once for grid line storage
*	2. The addresses edge[i] - glstore,  center[i] - slstore,
*	   and dl[i] - glstore are absolute numbers independent of processor.
*	   This is exploited in the pointer reconstruction for the topological
*	   grid that is transferred across multiple processors with the
*	   interface.
*
*	   Note that in the the case of the topological grid,  consistency
*	   requires that glstore be allocated by store().
*
*	All allocation of the edge, center, and dl fields should be done through
*	provided functions.
*
*	Besides supporting parallel virtual domains,  the precomputation
*	of the grid cell boundaries and centers should provide increased
*	efficiency in many cases.  Also note that the cell centers, edges
*	and spacings are related by the formulas:
*
*	dl[i][j] = edge[i][j+1] - edge[i][j]
*	center[i][j] = 0.5 * (edge[i][j+1] + edge[i][j])
*/

struct _RECT_GRID {
	float L[3];	  /* Lower corner of rectangle containing grid */
	float U[3];	  /* Upper corner of rectangle containing grid */
	float h[3];	  /* Average grid spacings in the grid	       */
	int   gmax[3];	  /* Number of grid blocks 		       */
	int   dim;	  /* Dimension of Grid */

		/* Specifications for virtual domains and variable grids */

	float GL[3];	  /* Lower corner of global grid */
	float GU[3];	  /* Upper corner of global grid */
	float VL[3];	  /* Lower corner of virtual domain */
	float VU[3];	  /* Upper corner of virtual domain */
	int   lbuf[3];	  /* Lower buffer zone width */
	int   ubuf[3];	  /* Upper buffer zone width */

		/* Specifications for variable mesh grids */

	float *edges[3];	/* Coordinate cell edges */
	float *centers[3];	/* Coordinate cell centers */
	float *dh[3];	   	/* Coordindate cell widths */
	float *glstore;	   	/* Storage for edges, centers and dh arrays */
	int   variable_mesh[3]; /* YES for variable dh in ith direction */

	struct _REMAP {
	    GEOMETRY_REMAP remap;
	    const char     *remap_name;
	    const char     *dnm[3],   /* Symbolic names for coord. directions*/
	                   *Dnm[3],
	                   *dirs[3];
	    float          area;      /* Total (Remappped) Computational Area */
	    /*grid cell area*/
	    float          (*Area)(const float*,const struct _RECT_GRID*);
	} Remap;
}; 
typedef struct _RECT_GRID RECT_GRID;

struct _COMM_BOX {
        int lmin[4][2][3];
        int lmax[4][2][3];
        int umin[4][2][3];
        int umax[4][2][3];
        float L0;
        float U0;
        float L1;
        float U1;
        int lx;
        int ly;
        int ux;
        int uy;
        int flag;
};
typedef struct _COMM_BOX COMM_BOX;


struct _RECT_BOX {
	struct _RECT_BOX  *prev;
	struct _RECT_BOX  *next;
	int bmin[MAXD];	  /* Lower bound of box indices */
	int bmax[MAXD];	  /* Upper bound of box indices */
	RECT_GRID *grid;
	int smin[MAXD];
	int smax[MAXD];
	int num_cross;
	struct _Cross *cross[100];
};
typedef struct _RECT_BOX RECT_BOX;
#if defined(__cplusplus)
typedef struct _RECT_GRID::_REMAP REMAP;
#else /* defined(__cplusplus) */
typedef	struct _REMAP REMAP;
#endif /* defined(__cplusplus) */


/*
*		Macros related to RECT_GRID
*	(regular grids only are currently supported)
*/

#define	grid_center_coord(i,gr)	(0.5 * ((gr)->L[i] + (gr)->U[i]))
#define	cell_index(p,i,gr)	irint(floor(((p)-(gr)->L[i])/(gr)->h[i]))
#define	cell_floor(p,i,gr)	irint(floor(((p)-(gr)->L[i])/(gr)->h[i]))
#define	cell_ceil(p,i,gr)	irint(ceil(((p)-(gr)->L[i])/(gr)->h[i]))
#define cell_center(indx,i,gr)	((gr)->L[i] + ((indx) + 0.5)*(gr)->h[i])
#define cell_edge(indx,i,gr)	((gr)->L[i] + (indx)*(gr)->h[i])
#define vd_cell_edge(indx,i,gr) ((gr)->VL[i] + (indx)*(gr)->h[i])
#define cell_width(index,i,gr)	((gr)->h[i])

/*
 * Note on SIDE.  Indicates the side of a hypersurface,  where the positive
 * side is the side of the hypersurface into which the normal points.
 * See the note regarding postive and negative components.
 */

enum _SIDE {
    UNKNOWN_SIDE  = -3,
    POSITIVE_SIDE =  0,
    NEGATIVE_SIDE =  1,
    ONEDGE	  =  2,
    ONVERTEX      =  3,
    COPLANAR      =  4
};
typedef enum _SIDE SIDE;

enum {
	PLUS_CUT  =  1,
	MINUS_CUT = -1
};

IMPORT	double	EPSILON;

 /* Direction names for oriented hypersurfaces */

typedef	enum {
	POSITIVE_ORIENTATION=1,
	NEGATIVE_ORIENTATION=2,
	ORIENTATION_NOT_SET=-1
} ORIENTATION;
typedef enum {
	COUNTER_CLOCK=1,
	CLOCKWISE=2,
	ANGLE_DIRECTION_NOT_SET=-1
} ANGLE_DIRECTION;

#define Opposite_orient(orient)						\
	(((orient)==POSITIVE_ORIENTATION)				\
	 ? NEGATIVE_ORIENTATION : POSITIVE_ORIENTATION)


#define Opposite_ang_dir(dir) ((dir) == CLOCKWISE ? COUNTER_CLOCK : CLOCKWISE)

 /* Angle orientation of curves */

#define curve_ang_oriented_l_to_r(ang_dir,orient)			      \
	((((ang_dir) == CLOCKWISE && (orient) == POSITIVE_ORIENTATION) ||     \
	((ang_dir) == COUNTER_CLOCK && (orient) == NEGATIVE_ORIENTATION)) ?  \
	YES : NO)

#define	Next_m3(n)		(((n) + 1) % 3)

#define	Prev_m3(n)		(((n) + 2) % 3)

/* Debugging printing of angles */

#define print_angle(mesg,ang,end)					\
	(void) printf("%s %g (%g deg)%s",mesg,ang,degrees(ang),end)

/*
*			Fourier Polynomial
*
*	           num_modes-1
*		     -----
*		     \                                 
*	z(p) = z0 +   \     A[k]*sin(<nu[k],p>) + phase[k])
*		      /                                
*		     /
*		     -----
*	             k = 0
*
*/

typedef struct {
	int   num_modes;
	int   dim;
	float **nu;
	float z0;
	float *A, *phase;
	float *L,*U; 		/* domain length in each direction */
} FOURIER_POLY;

/*
*			Legendre Polynomial
*
*		        max_degree
*			 -----
*		         \                                 
*		r(x) =    \     A[n]*P (x)
*			  /           n                    
*			 /
*			 -----
*		         n = 0
*
*	P(x) = Legendre polynomial of degree n.
*	 n
*
*	P(x) = 1
*        0
*
*	P(x) = x
*	 1
*
*	P(x) = 2*x*P(x) - P(x) - [x*P(x) - P(x)]/(n+1)   n >= 2
*	 n+1        n      n-1       n      n-1
*/

typedef struct {
	int   max_degree;
	float *A;
} LEGENDRE_POLY;

typedef struct {
        float *L;
        float *U;
} BDRY_BOX_PARAMS;

typedef struct {
	float 	      *cen;		/* Center of Ellipsoid */
	float 	      *rad;		/* Lengths of radii */
	float 	      ThetaS[2],
		      ThetaE[2];	/* Spherical coords of start and end */
	bool          closed;		/* Closed ellipsoid if YES */
	ORIENTATION   nor_orient;	/* Specifies inward or outward normal */
	FOURIER_POLY  *fpoly;		/* Fourier Perturbation factors */
	LEGENDRE_POLY *lpoly;		/* Legendre Perturbation factors */
	int 	      dim;		/* Dimension of embedding space */
	RECT_GRID     *gr;
} ELLIP_PARAMS;

typedef struct {
        float N[MAXD];		/* normal of the plane */
        float P[MAXD];		/* a point on the plane */
} PLANE_PARAMS;

typedef struct {
        float x0;
        float x1;
        float y;
        float z;
        float R;
        float rr;
} DUMBBELL_PARAMS;

typedef struct {
	/* equation for line is a*x + b*y = c */
        float a;		
        float b;		
        float c;		
} LINE_PARAMS;

typedef struct {
	/* equation for line is x^2/a^2 + y^2/b^2 = 1 */
        float x0;		
        float y0;		
        float a;		
        float b;		
} ELLIP2D_PARAMS;

typedef struct {
	/* equation for line is x^2/a^2 + y^2/b^2 = 1 */
        float x0;		
        float y0;		
        float r;		
        float w;		
        float h;		
} TDISK_PARAMS;

typedef struct {
	int num_cir;
	float **cen;
	float *rad;
} MC_PARAMS;

typedef struct {
        int           num_ellip;
        float         **cen;            /* Center of Ellipsoid 1 */
        float         **rad;            /* Lengths of radii 1 */
        ORIENTATION   nor_orient;       /* Specifies inward or outward normal */        int           dim;              /* Dimension of embedding space */
        RECT_GRID     *gr;
} M_ELLIP_PARAMS;

 /* Geometry EXPORTED Function Declarations*/

/*	geomutils.c*/
IMPORT	ANGLE_DIRECTION	fread_angle_direction(FILE*);
IMPORT	SIDE    Opposite_side(const SIDE);
IMPORT	double	dscaled_hypot(const double*,const double*,int);
IMPORT	const char *angle_direction_name(ANGLE_DIRECTION);
IMPORT	const char *orientation_name(ORIENTATION);
IMPORT	const char *side_name(SIDE);
IMPORT	float	_scaled_separation(const float*,const float*,const float*,int);
IMPORT	float	angle(float,float);
IMPORT	float	avg_angle_and_normalize(float,float);
IMPORT	float	cal_angle(const float*,const float*,const float*,
			  int,int,const float*);
IMPORT	float	distance_between_positions(const float*,const float*,int);
IMPORT	float	dscalar_product(const double*,const double*,int);
IMPORT	float	grid_size_in_direction(const float*,const float*,int);
IMPORT	float	mag_vector(const float*,int);
IMPORT	float	normalized_angle(float);
IMPORT	float	random_gaussian(float,float,unsigned short int [3]);
IMPORT	float	scalar_product(const float*,const float*,const int);
IMPORT	float	scaled_hypot(const float*,const float*,int);
IMPORT	float	triple_product(const float*,const float*,const float*,int);
IMPORT	float	vector_product(const float*,const float*,float*,int);
IMPORT	float	vector_product_on_points(const float*,const float*,
					 const float*,int,float*);
IMPORT	int	is_new_angle_smaller(float,float,float,float,int);
IMPORT	int	outside_point(const float*,const float*,const float*,int);
IMPORT	void	affine_fit(const float* const*,int,int,const float*,
			   float*,float**,float*);
IMPORT	void	fprint_angle_direction(FILE*,const char*,ANGLE_DIRECTION,
				       const char*);
IMPORT	void	fprint_general_vector(FILE*,const char*,const float*,
				      int,const char*);
IMPORT	void	fprint_orientation(FILE*,const char*,ORIENTATION,const char*);
IMPORT	void	sprint_general_vector(char*,const char*,const float*,
				      int,const char*);
IMPORT	void	print_angle_direction(const char*,ANGLE_DIRECTION,const char*);
IMPORT	void	print_general_vector(const char*,const float*,int,const char*);
IMPORT	void	print_side(const char*,SIDE,const char*);
IMPORT	void	print_orientation(const char*,ORIENTATION,const char*);

/*	igrid.c*/
IMPORT	GEOMETRY_REMAP	read_remap_from_string(const char*);
IMPORT	REMAP	*remap_info(void);
IMPORT	bool	adjust_top_grid_for_square(RECT_GRID*,const RECT_GRID*);
IMPORT	bool 	point_in_buffer(const float*,const RECT_GRID*);
IMPORT	bool 	is_rotational_symmetry(void);
IMPORT	bool	rect_in_which(const float*,int*,const RECT_GRID*);
IMPORT	float	grid_tolerance(const RECT_GRID*);
IMPORT	float	ident_Area(const float*,const RECT_GRID*);
IMPORT	int	set_grid_lines(RECT_GRID*);
IMPORT	void	set_rotational_symmetry(bool);
IMPORT	void	copy_rect_grid(RECT_GRID*,const RECT_GRID*);
IMPORT	void	fprint_rectangular_grid(FILE*,const RECT_GRID*);
IMPORT	void	free_grid_lines(RECT_GRID*);
IMPORT	void	i_init_remap_and_rect_grid(RECT_GRID*);
IMPORT	void	i_print_remap_values(void);
IMPORT	void	print_RECT_GRID_structure(const RECT_GRID*);
IMPORT	void	print_rectangular_grid(const RECT_GRID*);
IMPORT	void	read_rectangular_grid(const IO_TYPE*,RECT_GRID*,bool,REMAP*);
IMPORT	void	set_dual_grid(RECT_GRID*,const RECT_GRID*);
IMPORT	void	set_rect_grid(const float*,const float*,const float*,
			      const float*,const int*,const int*,const int*,
			      int,const REMAP*,RECT_GRID*);
IMPORT	void	set_remap(int,GEOMETRY_REMAP,REMAP*);
IMPORT	void	set_remap_identity(int,GEOMETRY_REMAP);
IMPORT	void	zoom_rect_grid(RECT_GRID*,const RECT_GRID*);
IMPORT  void	rect_grid_corner(const RECT_GRID*,const int*,float*);
IMPORT  void	rect_grid_center(const RECT_GRID*,const int*,float*);
IMPORT	void	init_topological_grid(RECT_GRID*,const RECT_GRID*);
IMPORT  void    set_remap_and_rect_grid(float*, float*,int*,GEOMETRY_REMAP,
					RECT_GRID*);

#endif /* !defined(_IGEOM_H) */
