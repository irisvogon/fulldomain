/*******************************************************************
 * 		CARTESIAN.c
 *******************************************************************/
#include "liquid.h"

//-------------------------------------------------------------------
//		INC_STATE
//-------------------------------------------------------------------
INC_STATE::INC_STATE(): m_P(0)
{
	m_U[0] = 0;
	m_U[1] = 0;
}
void INC_STATE::setZero(void)
{
	m_U[0] = 0;
	m_U[1] = 0;
	m_P = 0;
}


//------------------------------------------------------------------
// 		RECT_EDGE
// note: m_rectangle[] is uninitialized.
//------------------------------------------------------------------
RECT_EDGE::RECT_EDGE():m_comp(-1)
{
}

void RECT_EDGE::setRectangle(int i, int j)
{
	m_rectangle[0] = i;
	m_rectangle[1] = j;
}


//----------------------------------------------------------------
//		RECTANGLE
//----------------------------------------------------------------

//RECTANGLE::RECTANGLE()
RECTANGLE::RECTANGLE(): m_index(-1), m_comp(-1), m_visited(false)
{
}

void RECTANGLE::setCoords(
	double *coords,
	int dim)
{
	int i;
	for (i = 0; i < dim; ++i)
	    m_coords[i] = coords[i];
}
//--------------------------------------------------------------------------
// 		CARTESIAN
//--------------------------------------------------------------------------

CARTESIAN::~CARTESIAN()
{
}

//		CARTESIAN::initRectGrid()
// 1) there is a buffer zone outside the domain;
// 2) gmax[] means the number of cells;
void CARTESIAN::initRectGrid(Front *front)
{
	int i;
	
	makeGridIntfc();

	m_advection.initRectGrid(top_L,top_U,top_h,top_gmax);

	deleteGridIntfc();
}
//---------------------------------------------------------------
//	initMesh
// include the following parts
// 1) setup edges:		
// 2) setup cell_center
//---------------------------------------------------------------
void CARTESIAN::initMesh(void)
{
	int i,j,k, index;
	double coords[2];
	int num_cells,num_edges[MAXD],total_num_edges;
	int icoords[MAXD];
	int cell_index;

	// init vertices,edges & cell_center
	RECT_EDGE   	edge;
	RECTANGLE       rectangle;

	num_cells = 1;
	total_num_edges = 0;
	for (i = 0; i < dim; ++i)
	{
	    num_cells *= (top_gmax[i] + 1);
	    num_edges[i] = top_gmax[i] + 1;
	    for (j = 1; j < dim; ++j)
	    	num_edges[i] *= (top_gmax[(i+j)%dim] + 2);
	    total_num_edges += num_edges[i];
	}
	m_cell_center.insert(m_cell_center.end(),num_cells,rectangle);
	m_edges.insert(m_edges.end(),total_num_edges,edge);
	
	// setup vertices
	// left to right, down to up
	switch (dim)
	{
	case 2:
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	coords[0] = top_L[0] + top_h[0]*i;
	    	coords[1] = top_L[1] + top_h[1]*j;
		index = d_index2d(i,j,top_gmax);
	    	m_cell_center[index].setCoords(coords,dim);
	    	m_cell_center[index].icoords[0] = i;
	    	m_cell_center[index].icoords[1] = j;
	    }
	    break;
	case 3:
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	coords[0] = top_L[0] + top_h[0]*i;
	    	coords[1] = top_L[1] + top_h[1]*j;
	    	coords[2] = top_L[2] + top_h[2]*k;
		index = d_index3d(i,j,k,top_gmax);
	    	m_cell_center[index].setCoords(coords,dim);
	    	m_cell_center[index].icoords[0] = i;
	    	m_cell_center[index].icoords[1] = j;
	    	m_cell_center[index].icoords[2] = k;
	    }
	}
	
	// setup edge
	// first rows edges, then column edges
	// left to right, down to up
	// for cell, edges are ordered counter-clock wise, start from SOUTH
	index = 0;
	switch (dim)
	{
	case 2:
	    for (i = 0; i <= top_gmax[0]; i++)
	    for (j = 0; j <= top_gmax[1]+1; j++)
	    {
	    	m_edges[index].m_dir = RECT_EDGE_U;
		if (j == 0) 
		    m_edges[index].m_rectangle[0] = -1;
		else 
		{
		    cell_index = d_index2d(i,j-1,top_gmax);
		    m_edges[index].m_rectangle[0] = cell_index; //SOUTH
		    m_cell_center[cell_index].m_edge[2] = index; //NORTH
		}
		if (j == top_gmax[1]+1)
		    m_edges[index].m_rectangle[1] = -1;
		else 
		{
		    cell_index = d_index2d(i,j,top_gmax);
		    m_edges[index].m_rectangle[1] = cell_index; //NORTH
		    m_cell_center[cell_index].m_edge[0] = index; //SOUTH
		}
	    	index++;
	    }
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]+1; i++)
	    {
	    	m_edges[index].m_dir = RECT_EDGE_V;
		if (i == 0) 
		    m_edges[index].m_rectangle[0] = -1;
		else 
		{
		    cell_index = d_index2d(i-1,j,top_gmax);
		    m_edges[index].m_rectangle[0] = cell_index; //WEST
		    m_cell_center[cell_index].m_edge[1] = index; //EAST
		}
		if (i == top_gmax[0]+1)
		    m_edges[index].m_rectangle[1] = -1;
		else 
		{
		    cell_index = d_index2d(i,j,top_gmax);
		    m_edges[index].m_rectangle[1] = cell_index; //EAST
		    m_cell_center[cell_index].m_edge[3] = index; //WEST
		}
	    	index++;
	    }
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; i++)
	    for (j = 0; j <= top_gmax[1]+1; j++)
	    for (k = 0; k <= top_gmax[2]+1; k++)
	    {
	    	m_edges[index].m_dir = RECT_EDGE_U;
	    	index++;
	    }
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (k = 0; k <= top_gmax[2]+1; k++)
	    for (i = 0; i <= top_gmax[0]+1; i++)
	    {
	    	m_edges[index].m_dir = RECT_EDGE_V;
	    	index++;
	    }
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (i = 0; i <= top_gmax[0]+1; i++)
	    for (j = 0; j <= top_gmax[1]+1; j++)
	    {
	    	m_edges[index].m_dir = RECT_EDGE_W;
	    	index++;
	    }
	    break;
	}

	setComponent();
}

void CARTESIAN::setComponent(void)
{
	int i;
	double coords[2];
	
	// cell center components
	for (i = 0; i < m_cell_center.size(); i++)
	{
	    m_cell_center[i].m_comp = 
	    		getComponent(m_cell_center[i].icoords);
	}
}

// smooth discontinuous properties, such as the density, mu, force.
void CARTESIAN::setSmoothedProperties(void)
{
	int i,j,k; 
	int index, index0, index1, index2, index3;
	int comp, c[4];
	double center[2], point[2], phi, H;
	bool is_near_intfc_pt;

	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index  = d_index2d(i,j,top_gmax);			
	    comp  = m_cell_center[index].m_comp;
	    index0 = d_index2d(i,j-1,top_gmax);			
	    c[0] = m_cell_center[index0].m_comp;
	    index1 = d_index2d(i+1,j,top_gmax);			
	    c[1] = m_cell_center[index1].m_comp;
	    index2 = d_index2d(i,j+1,top_gmax);			
	    c[2] = m_cell_center[index2].m_comp;
	    index3 = d_index2d(i-1,j,top_gmax);			
	    c[3] = m_cell_center[index3].m_comp;
			
	    is_near_intfc_pt = NO;
	    for (k = 0; k < 4; ++k)
	    {
	    	if (comp != c[k] && !is_exterior_comp(c[k],front->interf))
		    is_near_intfc_pt = YES;
	    }
	    if (is_near_intfc_pt)
	    {
		getRectangleCenter(index, center);
		getNearestInterfacePoint(center, point);	

		phi = getDistance(center,point);
		if (comp == m_comp[0])
		    phi *= -1;
			
		H = getSmoothingFunction(phi);
		m_cell_center[index].m_state.m_mu = m_mu[0]  + 
					(m_mu[1]-m_mu[0])*H;
		m_cell_center[index].m_state.m_rho = m_rho[0] + 
					(m_rho[1]-m_rho[0])*H; 
	    }
	    else
	    {
		m_cell_center[index].m_state.m_mu  = m_mu[comp];
		m_cell_center[index].m_state.m_rho = m_rho[comp];
	    }

	}
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{	
	    index  = d_index2d(i,j,top_gmax);
	    array[index] = m_cell_center[index].m_state.m_mu;
	}
	scatMeshArray();
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{	
	    	index  = d_index2d(i,j,top_gmax);
	    	m_cell_center[index].m_state.m_mu = array[index];
	}
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    array[index] = m_cell_center[index].m_state.m_rho;
	}
	scatMeshArray();
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    m_cell_center[index].m_state.m_rho = array[index];
	}
}

// [1] M Kang, R Fedkiw and X.-D Liu, A boundary Condition capturing method for
//     muliphase incompressible flow.
//
double CARTESIAN::getSmoothingFunction(double phi)
{
	/*if (phi<0)
		return 0;
	else
		return 1;*/
	// Heaviside function [1]
	if (phi < -m_smoothing_radius)	
		return 0;
	else if (phi > m_smoothing_radius)
		return 1;
	else
		return 1.0/2 + phi/(2*m_smoothing_radius) + 1/(2*PI)*sin(PI*phi/m_smoothing_radius);
}

void CARTESIAN::setInitialCondition(void)
{
	int i;
	double coords[2];
	INC_STATE state;

	m_rho[0] = 1.0;		m_rho[1] = 0.1;
	m_mu[0]  = PI/50.0;	m_mu[1]  = PI/50.0;
	m_smoothing_radius = top_h[0] < top_h[1] ? top_h[1] : top_h[0];
	m_smoothing_radius *= 2.0;

	// cell_center
	for (i = 0; i < m_cell_center.size(); i++)
	{
	    getRectangleCenter(i, coords);
	    computeExactSolution(coords,m_cell_center[i].m_state);
	}
}

void CARTESIAN::setIndexMap(void)
{
	static bool first = YES;
	int i,j,k,index;
	int llbuf[MAXD],uubuf[MAXD];

	if (first)
	{
	    first = NO;
	    switch (dim)
	    {
	    case 2:
	    	bi_array(&ij_to_I,top_gmax[0]+1,top_gmax[1]+1,INT);
	    	break;
	    case 3:
	    	tri_array(&ijk_to_I,top_gmax[0]+1,top_gmax[1]+1,top_gmax[2]+1,
					INT);
	    	break;
	    }
	}

	index = 0;
	for (i = 0; i < dim; ++i)
	{
	    llbuf[i] = lbuf[i] != 0 ? lbuf[0] : 1;
	    uubuf[i] = ubuf[i] != 0 ? ubuf[0] : 1;
	}
	switch (dim)
	{
	case 2:
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
	    	ij_to_I[i][j] = index + ilower;
	    	index++;
	    }
	    scatter_cell_index(front,llbuf,uubuf,(POINTER)ij_to_I);
	    break;
	case 3:
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
	    	ijk_to_I[i][j][k] = index + ilower;
	    	index++;
	    }
	    scatter_cell_index(front,llbuf,uubuf,(POINTER)ijk_to_I);
	    break;
	}
}

// no need to consider the jump condition since the velocity is continuous.
void CARTESIAN::interpolateEdgeVelocity(void)
{
	int i, j;
	INTERFACE *intfc = front->interf;

	int index0, index1;

	switch (rect_boundary_type(intfc,1,0))
	{
	case NEUMANN_BOUNDARY:
	    j = 0;
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	m_cell_center[d_index2d(i,j,top_gmax)].m_state.m_U[0] = - 
	    	    m_cell_center[d_index2d(i,1,top_gmax)].m_state.m_U[0];
	    	m_cell_center[d_index2d(i,j,top_gmax)].m_state.m_U[1] = - 
	    	    m_cell_center[d_index2d(i,1,top_gmax)].m_state.m_U[1];
	    }
	    break;
	}
	switch (rect_boundary_type(intfc,1,1))
	{
	case NEUMANN_BOUNDARY:
	    j = top_gmax[1];
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	m_cell_center[d_index2d(i,j,top_gmax)].m_state.m_U[0] = 
	    	    - m_cell_center[d_index2d(i,j-1,top_gmax)].m_state.m_U[0];
	    	m_cell_center[d_index2d(i,j,top_gmax)].m_state.m_U[1] = 
	    	    - m_cell_center[d_index2d(i,j-1,top_gmax)].m_state.m_U[1];
	    }
	}

	for (i=0; i<m_edges.size(); i++)
	{
	    index0 = m_edges[i].m_rectangle[0];
	    index1 = m_edges[i].m_rectangle[1];

	    if (index0 == -1 || index1 == -1) continue;

	    m_edges[i].m_state.m_U[0] = 
			0.5*(m_cell_center[index0].m_state.m_U[0] + 
			m_cell_center[index1].m_state.m_U[0]);
	    m_edges[i].m_state.m_U[1] = 
			0.5*(m_cell_center[index0].m_state.m_U[1] + 
			m_cell_center[index1].m_state.m_U[1]);
	}
}

void CARTESIAN::computeAdvectionSimple(void)
{
	int i, j, k, index, index0, index1, index2, index3, size;
	INC_STATE state;

	double *u, *v;
	size = (top_gmax[0]+1)*(top_gmax[1]+1);
	uni_array(&u,size,sizeof(double));
	uni_array(&v,size,sizeof(double));

	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    u[index] = m_cell_center[index].m_state.m_U[0];
	    v[index] = m_cell_center[index].m_state.m_U[1];
	}

	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{	
	    index  = d_index2d(i,j,top_gmax);
	    index0 = d_index2d(i,j-1,top_gmax);
	    index1 = d_index2d(i+1,j,top_gmax);
	    index2 = d_index2d(i,j+1,top_gmax);
	    index3 = d_index2d(i-1,j,top_gmax);
	    m_cell_center[index].m_state.m_U[0] += - m_dt* 
			  (u[index]*(u[index1] - u[index3])/(2*top_h[0])
			   +v[index]*(u[index2] - u[index0])/(2*top_h[1]));
	    m_cell_center[index].m_state.m_U[1] += - m_dt* 
			  (u[index]*(v[index1] - v[index3])/(2*top_h[0])
			   +v[index]*(v[index2] - v[index0])/(2*top_h[1]));
	}
	for (k = 0; k < 2; ++k)
	{
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {	
	    	index  = d_index2d(i,j,top_gmax);
	    	array[index] = m_cell_center[index].m_state.m_U[k];
	    }
	    scatMeshArray();
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {	
	    	index  = d_index2d(i,j,top_gmax);
	    	m_cell_center[index].m_state.m_U[k] = array[index];
	    }
	}
	free_these(2,u,v);
}


double CARTESIAN::computeVelocityDivergence(int indexRectangle)
{
	int i, j, index0, index1, index2, index3, comp;
	double coords0[2], coords1[2], coeff0, coeff1, rhs[2][2], I, J, ret;

	getRectangleIndex(indexRectangle, i, j);
	ret = 0;
	ret += 1/(2*top_h[0]) * 
		(m_cell_center[d_index2d(i+1,j,top_gmax)].m_state.m_U[0] - 
		m_cell_center[d_index2d(i-1,j,top_gmax)].m_state.m_U[0]);
	ret += 1/(2*top_h[1]) * 
		(m_cell_center[d_index2d(i,j+1,top_gmax)].m_state.m_U[1] - 
		m_cell_center[d_index2d(i,j-1,top_gmax)].m_state.m_U[1]);
	return ret;
}

void CARTESIAN::compPressGradWithSmoothProperty(int index,double *gradP)
{
	int index0,index1,index2,index3,i,j;

	getRectangleIndex(index, i, j);
	index0 = d_index2d(i,j-1,top_gmax);
	index1 = d_index2d(i+1,j,top_gmax);
	index2 = d_index2d(i,j+1,top_gmax);
	index3 = d_index2d(i-1,j,top_gmax);
	// Px	
	gradP[0] = 1/(2*top_h[0])*(m_cell_center[index1].m_state.m_P - 
			m_cell_center[index3].m_state.m_P);
	// Py
	gradP[1] = 1/(2*top_h[1])*(m_cell_center[index2].m_state.m_P - 
			m_cell_center[index0].m_state.m_P);
}

void CARTESIAN::compDiffWithSmoothProperty(void)
{
	int comp,index,index0,index1,index2,index3,
		index4,index5,index6,index7,size;
	int I,I0,I1,I2,I3,I4,I5,I6,I7;
	double coords[2],coeff[8],mu[4],mu0,rho,corner[8],rhs;
	INC_STATE state;
	int i,j,k;
	INTERFACE *intfc = front->interf;

	PETSc solver;

	setIndexMap();

	size = 1;
	for (i = 0; i < dim; ++i) size *= gmax[i];
	solver.Create(2*ilower, 2*iupper-1, 9, 0);	// two buffer zones

	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    index0 = d_index2d(i,j-1,top_gmax);
	    index1 = d_index2d(i+1,j,top_gmax);
	    index2 = d_index2d(i,j+1,top_gmax);
	    index3 = d_index2d(i-1,j,top_gmax);
	    index4 = d_index2d(i-1,j-1,top_gmax);
	    index5 = d_index2d(i+1,j-1,top_gmax);
	    index6 = d_index2d(i+1,j+1,top_gmax);
	    index7 = d_index2d(i-1,j+1,top_gmax);
	    I  = ij_to_I[i][j];
	    I0 = ij_to_I[i][j-1];
	    I1 = ij_to_I[i+1][j];
	    I2 = ij_to_I[i][j+1];
	    I3 = ij_to_I[i-1][j];
	    I4 = ij_to_I[i-1][j-1];
	    I5 = ij_to_I[i+1][j-1];
	    I6 = ij_to_I[i+1][j+1];
	    I7 = ij_to_I[i-1][j+1];

	    for (k = 0; k < 8; k++)
		corner[k] = 1;

	    if (j == 1 && 
	    	rect_boundary_type(intfc,1,0) == NEUMANN_BOUNDARY)
	    {
		index0 = index;		corner[0] = -1;
		index5 = index1;	corner[5] = -1;
		index4 = index3;	corner[4] = -1;
		I0 = I;
		I5 = I1;
		I4 = I3;
	    }
	    if (j == (top_gmax[1]-1) &&
	    	rect_boundary_type(intfc,1,1) == NEUMANN_BOUNDARY)
	    {
		index2 = index;		corner[2] = -1;
		index6 = index1;	corner[6] = -1;
		index7 = index3;	corner[7] = -1;
		I2 = I;
		I6 = I1;
		I7 = I3;
	    }

	    mu0   = m_cell_center[index].m_state.m_mu;
	    rho   = m_cell_center[index].m_state.m_rho;
	    mu[0] = 1.0/2 * (mu0 + m_cell_center[index0].m_state.m_mu);// down
	    mu[1] = 1.0/2 * (mu0 + m_cell_center[index1].m_state.m_mu);// right
	    mu[2] = 1.0/2 * (mu0 + m_cell_center[index2].m_state.m_mu);// up
	    mu[3] = 1.0/2 * (mu0 + m_cell_center[index3].m_state.m_mu);// left

	    coeff[0] = 0.5*m_dt/rho * mu[0]/(top_h[1]*top_h[1]);	// down
	    coeff[1] = 0.5*m_dt/rho * mu[1]/(top_h[0]*top_h[0]);	// right
	    coeff[2] = 0.5*m_dt/rho * mu[2]/(top_h[1]*top_h[1]);	// up	
	    coeff[3] = 0.5*m_dt/rho * mu[3]/(top_h[0]*top_h[0]);	// left
	    coeff[4] = 0.5*m_dt/rho * m_cell_center[index0].m_state.m_mu/
	    				(4*top_h[0]*top_h[1]);
	    coeff[5] = 0.5*m_dt/rho * m_cell_center[index1].m_state.m_mu/
	    				(4*top_h[0]*top_h[1]);
	    coeff[6] = 0.5*m_dt/rho * m_cell_center[index2].m_state.m_mu/
	    				(4*top_h[0]*top_h[1]);
	    coeff[7] = 0.5*m_dt/rho * m_cell_center[index3].m_state.m_mu/
	    				(4*top_h[0]*top_h[1]);


	    // source term
	    getRectangleCenter(index, coords);
	    computeSourceTerm(coords, state);

	    solver.Add_A(I*2,I*2,1+(2*coeff[1]+2*coeff[3]+coeff[0]+coeff[2])); 
	    solver.Add_A(I*2,I0*2,-coeff[0]*corner[0]);
	    solver.Add_A(I*2,I1*2,-2*coeff[1]*corner[1]);
	    solver.Add_A(I*2,I2*2,-coeff[2]*corner[2]);
	    solver.Add_A(I*2,I3*2,-2*coeff[3]*corner[3]);
	    solver.Add_A(I*2,I4*2+1,-coeff[4]*corner[4]); 
	    solver.Add_A(I*2,I5*2+1,coeff[4]*corner[5]); 
	    solver.Add_A(I*2,I6*2+1,-coeff[6]*corner[6]); 
	    solver.Add_A(I*2,I7*2+1,coeff[6]*corner[7]); 
	    rhs = m_dt*state.m_U[0];
	    rhs += (1 - (2*coeff[1]+2*coeff[3]+coeff[0]+coeff[2]))*
	    		m_cell_center[index].m_state.m_U[0];
	    rhs +=   coeff[0]*m_cell_center[index0].m_state.m_U[0]*corner[0];
	    rhs += 2*coeff[1]*m_cell_center[index1].m_state.m_U[0]*corner[1];
	    rhs +=   coeff[2]*m_cell_center[index2].m_state.m_U[0]*corner[2];
	    rhs += 2*coeff[3]*m_cell_center[index3].m_state.m_U[0]*corner[3];
	    rhs +=   coeff[4]*m_cell_center[index4].m_state.m_U[1]*corner[4];
	    rhs += - coeff[4]*m_cell_center[index5].m_state.m_U[1]*corner[5];
	    rhs +=   coeff[6]*m_cell_center[index6].m_state.m_U[1]*corner[6];
	    rhs += - coeff[6]*m_cell_center[index7].m_state.m_U[1]*corner[7];
	    solver.Add_b(I*2, rhs);
	    solver.Add_A(I*2+1,I*2+1,     
			1 + (coeff[1]+coeff[3]+2*coeff[0]+2*coeff[2])); 
	    solver.Add_A(I*2+1,I0*2+1,-2*coeff[0]*corner[0]);
	    solver.Add_A(I*2+1,I1*2+1,-coeff[1]*corner[1]);
	    solver.Add_A(I*2+1,I2*2+1,-2*coeff[2]*corner[2]);
	    solver.Add_A(I*2+1,I3*2+1,-coeff[3]*corner[3]);
	    solver.Add_A(I*2+1,I4*2,-coeff[7]*corner[4]); 
	    solver.Add_A(I*2+1,I5*2,coeff[5]*corner[5]); 
	    solver.Add_A(I*2+1,I6*2,-coeff[5]*corner[6]); 
	    solver.Add_A(I*2+1,I7*2,coeff[7]*corner[7]); 
	    rhs = m_dt * state.m_U[1];
	    rhs += (1 - (coeff[1]+coeff[3]+2*coeff[0]+2*coeff[2]))*
	    		m_cell_center[index].m_state.m_U[1];
	    rhs += 2*coeff[0]*m_cell_center[index0].m_state.m_U[1]*corner[0];
	    rhs +=   coeff[1]*m_cell_center[index1].m_state.m_U[1]*corner[1];
	    rhs += 2*coeff[2]*m_cell_center[index2].m_state.m_U[1]*corner[2];
	    rhs +=   coeff[3]*m_cell_center[index3].m_state.m_U[1]*corner[3];
	    rhs +=   coeff[7]*m_cell_center[index4].m_state.m_U[0]*corner[4];
	    rhs += - coeff[5]*m_cell_center[index5].m_state.m_U[0]*corner[5];
	    rhs +=   coeff[5]*m_cell_center[index6].m_state.m_U[0]*corner[6];
	    rhs += - coeff[7]*m_cell_center[index7].m_state.m_U[0]*corner[7];
	    solver.Add_b(I*2+1, rhs);
	}

	solver.SetMaxIter(500);
	solver.SetTol(1e-6);
	solver.Solve();

	// get back the solution
	double *x;
	uni_array(&x,2*size,sizeof(double));
	solver.Get_x(x);

	int num_iter;
	double rel_residual;
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);

	if (debugging("PETSc"))
	    (void) printf("CARTESIAN::compDiffWithSmoothProperty: "
	       		"num_iter = %d, rel_residual = %le. \n", 
			num_iter,rel_residual); 

	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    I = ij_to_I[i][j];
	    index = d_index2d(i,j,top_gmax);
	    if (I>=0)
	    {
		m_cell_center[index].m_state.m_U[0] = 
				x[I*2-ilower*2];
		m_cell_center[index].m_state.m_U[1] = 
				x[I*2+1-ilower*2];
	    }
	    else
		printf("CARTESIAN::compDiffWithSmoothProperty: ???\n");
	}
	for (k = 0; k < 2; ++k)
	{
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {	
	    	index  = d_index2d(i,j,top_gmax);
	    	array[index] = m_cell_center[index].m_state.m_U[k];
	    }
	    scatMeshArray();
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {	
	    	index  = d_index2d(i,j,top_gmax);
	    	m_cell_center[index].m_state.m_U[k] = array[index];
	    }
	}

	free_these(1,x);
}

void CARTESIAN::compProjWithSmoothProperty(void)
{
	int index, index0, index1, index2, index3, size;
	double rhs, coeff[4], rho[4], rho0;
	int I,I0,I1,I2,I3;
	int i,j;
	INTERFACE *intfc = front->interf;

	PETSc solver;
	solver.Create(ilower, iupper-1, 5, 0);

	setIndexMap();

	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    index0 = d_index2d(i,j-1,top_gmax);
	    index1 = d_index2d(i+1,j,top_gmax);
	    index2 = d_index2d(i,j+1,top_gmax);
	    index3 = d_index2d(i-1,j,top_gmax);
	    I = ij_to_I[i][j];
	    I0 = ij_to_I[i][j-1];
	    I1 = ij_to_I[i+1][j];
	    I2 = ij_to_I[i][j+1];
	    I3 = ij_to_I[i-1][j];
	
	    if (j == 1 &&
	    	rect_boundary_type(intfc,1,0) == NEUMANN_BOUNDARY)
	    {
		index0 = index;
		I0 = I;
	    }
	    if (j == (top_gmax[1]-1) &&
	    	rect_boundary_type(intfc,1,1) == NEUMANN_BOUNDARY)
	    {
		index2 = index;
		I2 = I;
	    }


	    rho0   = m_cell_center[index].m_state.m_rho;
	    rho[0] = 1.0/2*(rho0 + m_cell_center[index0].m_state.m_rho);// down
	    rho[1] = 1.0/2*(rho0 + m_cell_center[index1].m_state.m_rho);// right
	    rho[2] = 1.0/2*(rho0 + m_cell_center[index2].m_state.m_rho);// up
	    rho[3] = 1.0/2*(rho0 + m_cell_center[index3].m_state.m_rho);// left
	    coeff[0] = 1/rho[0]/(top_h[1]*top_h[1]); 
	    coeff[1] = 1/rho[1]/(top_h[0]*top_h[0]); 
	    coeff[2] = 1/rho[2]/(top_h[1]*top_h[1]); 
	    coeff[3] = 1/rho[3]/(top_h[0]*top_h[0]); 

	    solver.Add_A(I,I,-coeff[0]-coeff[1]-coeff[2]-coeff[3]);
	    solver.Add_A(I,I0,coeff[0]);
	    solver.Add_A(I,I1,coeff[1]);
	    solver.Add_A(I,I2,coeff[2]);
	    solver.Add_A(I,I3,coeff[3]);

	    rhs = 1/m_dt*computeVelocityDivergence(index);
	    solver.Add_b(I,rhs);
	}
	
	double gradP[2];

	solver.SetMaxIter(500);
	solver.SetTol(1e-6);
	solver.Solve_withPureNeumann();

	double *x;
	uni_array(&x,m_cell_center.size(),sizeof(double));
	solver.Get_x(x);

	int num_iter = 0;
	double rel_residual = 0;
	solver.GetNumIterations(&num_iter);
	solver.GetFinalRelativeResidualNorm(&rel_residual);

	if (debugging("PETSc"))
	    (void) printf("CARTESIAN::compProjWithSmoothProperty: "
	       		"num_iter = %d, rel_residual = %le \n", 
			num_iter, rel_residual);
	
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    I = ij_to_I[i][j];

	    if (j == 0 &&
	    	rect_boundary_type(intfc,1,0) == NEUMANN_BOUNDARY)
	    {
		I0 = ij_to_I[i][1];
		m_cell_center[index].m_state.m_P = x[I0-ilower];
	    }
	    else if (j == (top_gmax[1]) &&
	    	rect_boundary_type(intfc,1,1) == NEUMANN_BOUNDARY)
	    {
		I0 = ij_to_I[i][top_gmax[1]-1];
		m_cell_center[index].m_state.m_P = x[I0-ilower];
	    }
	    else
	    	m_cell_center[index].m_state.m_P = x[I-ilower];
	}
	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index  = d_index2d(i,j,top_gmax);
	    array[index] = m_cell_center[index].m_state.m_P;
	}
	scatMeshArray();
	for (j = 0; j <= top_gmax[1]; j++)
	for (i = 0; i <= top_gmax[0]; i++)
	{
	    	index  = d_index2d(i,j,top_gmax);
	    	m_cell_center[index].m_state.m_P = array[index];
	}
	free_these(1,x);
}


void CARTESIAN::compNewVeloWithSmoothProperty(void)
{
	int i, j, k, index0, index;
	double gradP[2], rho;

	for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
	{
	    index = d_index2d(i,j,top_gmax);
	    rho = m_cell_center[index].m_state.m_rho;
	    compPressGradWithSmoothProperty(index,gradP);
	    m_cell_center[index].m_state.m_U[0] -= m_dt*1/rho*gradP[0];
	    m_cell_center[index].m_state.m_U[1] -= m_dt*1/rho*gradP[1];
	}
	for (k = 0; k < 2; ++k)
	{
	    for (j = jmin; j <= jmax; j++)
            for (i = imin; i <= imax; i++)
	    {	
	    	index  = d_index2d(i,j,top_gmax);
	    	array[index] = m_cell_center[index].m_state.m_U[k];
	    }
	    scatMeshArray();
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {	
	    	index  = d_index2d(i,j,top_gmax);
	    	m_cell_center[index].m_state.m_U[k] = array[index];
	    }
	}
}


void CARTESIAN::computeSourceTerm(double *coords, INC_STATE &state) 
{
	state.m_U[0] = 0;
	state.m_U[1] = -100;

	state.m_P = HUGE_VAL;
}
void CARTESIAN::computeSourceTerm(double *coords, double t, INC_STATE &state) 
{
	computeSourceTerm(coords, state);
}

void CARTESIAN::computeExactSolution(double *coords, INC_STATE &state) 
{
	state.setZero();
}

// for initial condition: 
// 		setInitialCondition();	
// this function should be called before solve()
// for the source term of the momentum equation: 	
// 		computeSourceTerm();
void CARTESIAN::solve(double dt)
{
	m_dt = dt;

	start_clock("solve");
	makeGridIntfc();

	setGlobalIndex();
	setComponent();
	setSmoothedProperties();
	
	interpolateEdgeVelocity();

	// 1) solve for intermediate velocity
	computeAdvectionSimple();
	
	start_clock("compDiffWithSmoothProperty");
	compDiffWithSmoothProperty();
	stop_clock("compDiffWithSmoothProperty");
	interpolateEdgeVelocity();

	// 2) projection step
	start_clock("compProjWithSmoothProperty");
	compProjWithSmoothProperty();
	stop_clock("compProjWithSmoothProperty");
	compNewVeloWithSmoothProperty();
	deleteGridIntfc();
	stop_clock("solve");
}


// check http://en.wikipedia.org/wiki/Bilinear_interpolation
void CARTESIAN::getVelocity(double *p, double *U)
{
	// locate the point
	int icoords[MAXD];
	int i,j,k;
	double c1[MAXD], c2[MAXD], denominator;

	if (!rect_in_which(p,icoords,top_grid))
	{
	    for (i=0; i<2; i++)
	    {
	    	U[i] = 0.0;
	    }
	    return;
	}

	switch (dim)
	{
	case 2:
	    int index11, index21, index12, index22;
	    i = icoords[0];
	    j = icoords[1];
	    index11 = d_index2d(i,j,top_gmax);
	    index12 = d_index2d(i,j+1,top_gmax);
	    index21 = d_index2d(i+1,j,top_gmax);
	    index22 = d_index2d(i+1,j+1,top_gmax);

	    getRectangleCenter(index11, c1);
	    getRectangleCenter(index22, c2);

	    denominator = (c2[0]-c1[0])*(c2[1]-c1[1]);
	    for (i=0; i<2; i++)
	    {
		U[i] =  m_cell_center[index11].m_state.m_U[i]/denominator
				*(c2[0]-p[0])*(c2[1]-p[1]) 
			  + m_cell_center[index21].m_state.m_U[i]/denominator 
			  	*(p[0]-c1[0])*(c2[1]-p[1])
			  + m_cell_center[index12].m_state.m_U[i]/denominator
			  	*(c2[0]-p[0])*(p[1]-c1[1])
			  + m_cell_center[index22].m_state.m_U[i]/denominator
			  	*(p[0]-c1[0])*(p[1]-c1[1]);
	    }
	    break;
	}
}

void CARTESIAN::getRectangleIndex(int index, int &i, int &j)
{
	i = m_cell_center[index].icoords[0];
	j = m_cell_center[index].icoords[1];
}

void CARTESIAN::getRectangleIndex(int index, int &i, int &j, int &k)
{
	i = m_cell_center[index].icoords[0];
	j = m_cell_center[index].icoords[1];
	k = m_cell_center[index].icoords[2];
}


int CARTESIAN::getRectangleComponent(int index)
{	
	return getComponent(m_cell_center[index].icoords);
}

void CARTESIAN::getRectangleCenter(
	int index, 
	double *coords)
{
	int i;
	for (i = 0; i < dim; ++i)
	    coords[i] = m_cell_center[index].m_coords[i];
}

void CARTESIAN::getRectangleCenter(
	int index0, 
	int index1, 
	double *coords)
{
	int i;
	for (i = 0; i < dim; ++i)
	{
	    coords[i] = 0.5*(m_cell_center[index0].m_coords[i] +
	    		     m_cell_center[index1].m_coords[i]);
	}
}

int CARTESIAN::getInteger(double i)
{
	return int(i+0.01);
}

bool CARTESIAN::isInteger(double i)
{
	double ret;
	ret = floor(i+0.1) - i;
	if (fabs(ret)<0.1)
	    return true;
	else
	    return false;
}


double CARTESIAN::getDistance(double *c0, double *c1)
{
	return sqrt( (c0[0]-c1[0])*(c0[0]-c1[0])
		    +(c0[1]-c1[1])*(c0[1]-c1[1]) );
}


// input : p[]
// output: q[]

void CARTESIAN::getNearestInterfacePoint(
	double *p, 
	double *q)
{
	INTERFACE *intfc = front->interf;
	double t;
	HYPER_SURF_ELEMENT *phse;
	HYPER_SURF *phs;
	nearest_interface_point(p,getComponent(p),intfc,NO_BOUNDARIES,
				NULL,q,&t,&phse,&phs);
}

int CARTESIAN::getComponent(
	double *p)
{
	return component(p,front->interf);
}

int CARTESIAN::getComponent(
	int *icoords)
{
	int index;
	switch (dim)
	{
	case 2:
	    index = d_index2d(icoords[0],icoords[1],top_gmax);
	    return top_comp[index];
	case 3:
	    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
	    return top_comp[index];
	}
}

void CARTESIAN::save(char *filename)
{
	
	RECT_GRID *rect_grid = front->rect_grid;
	INTERFACE *intfc    = front->interf;
		
	int i, j;
	int xmax = rect_grid->gmax[0];
	int ymax = rect_grid->gmax[1];
	double x, y;
	
	FILE *hfile = fopen(filename, "w");
	if(hfile==NULL)
	{
		printf("\n can't open %s in "
		       "SaveAsTecplot_rect_grid_and_interface().", filename);
		exit(0);
	}
	
	// secondly print out the interface
		
	if(exists_interface(intfc))
	{
	    CURVE		**cur;
	    CURVE		*curve;
	    BOND		*bond;
			
	    for(cur=intfc->curves; cur && *cur; cur++)
	    {
		curve = *cur;
		fprintf(hfile, "ZONE I=%d J=%d F=POINT \n", 
				curve->num_points, 1);
		bond=curve->first;
		fprintf(hfile, "%.4f %.4f \n",bond->start->_coords[0], 
				bond->start->_coords[1]);
		for(bond=curve->first; bond!=NULL; bond=bond->next)
		    fprintf(hfile, "%.4f %.4f \n",bond->end->_coords[0], 
		    		bond->end->_coords[1]);
		}					
	}		
	fclose(hfile);
}

CARTESIAN::CARTESIAN(Front &front):front(&front)
{
}

void CARTESIAN::makeGridIntfc()
{
	static bool first = YES;
	INTERFACE *grid_intfc;
	Table *T;

	FrontMakeGridIntfc(front);

	grid_intfc = front->grid_intfc;
	top_grid = &topological_grid(grid_intfc);
	lbuf = front->rect_grid->lbuf;
	ubuf = front->rect_grid->ubuf;
	top_gmax = top_grid->gmax;
	top_L = top_grid->L;
	top_U = top_grid->U;
	top_h = top_grid->h;
	dim = grid_intfc->dim;
	T = table_of_interface(grid_intfc);
	top_comp = T->components;

	switch (dim)
	{
	case 2:
	    if (first)
	    {
	    	uni_array(&array,
			(top_gmax[0]+1)*(top_gmax[1]+1),FLOAT);
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    break;
	case 3:
	    if (first)
	    {
	    	uni_array(&array,
			(top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1),FLOAT);
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    kmin = (lbuf[2] == 0) ? 1 : lbuf[2];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    kmax = (ubuf[2] == 0) ? top_gmax[2] - 1 : top_gmax[2] - ubuf[2];
	    break;
	}
}

void CARTESIAN::deleteGridIntfc()
{
	FrontFreeGridIntfc(front);
}

void CARTESIAN::scatMeshArray()
{
	FrontScatExpDualGridArray(array,front);
}

void CARTESIAN::setGlobalIndex()
{
	int i;
	int num_nodes = pp_numnodes();
	int myid = pp_mynode();
	uni_array(&n_dist,num_nodes,sizeof(int));

	gmax = front->rect_grid->gmax;
	NLblocks = 1;
	for (i = 0; i < dim; ++i)
	    NLblocks *= gmax[i];
	for (i = 0; i < num_nodes; ++i) n_dist[i] = 0;
	n_dist[myid] = NLblocks;
	pp_global_imax(n_dist,num_nodes);
	ilower = 0;
        iupper = n_dist[0];
        for (i = 1; i <= myid; i++)
        {
            ilower += n_dist[i-1];
            iupper += n_dist[i];
        }	
}


void CARTESIAN::checkStates()
{
	static double ***states;
	int i,j,k,index;
	char fname[100];
	static int step,count;
	FILE *sfile;

	if (states == NULL)
	    tri_array(&states,top_gmax[0]+1,top_gmax[1]+1,5,sizeof(double));
	if (step != front->step)
	{
	    step = front->step;
	    count = 0;
	}
	sprintf(fname,"test/state-%d-%d-%d.%d",step,count,pp_numnodes(),
			pp_mynode());
	count++;

	for (j = 0; j < top_gmax[1]; j++)
	for (i = 0; i < top_gmax[0]; i++)
	{	
	    index  = d_index2d(i,j,top_gmax);
	    states[i][j][0] = m_cell_center[index].m_state.m_U[0];
	    states[i][j][1] = m_cell_center[index].m_state.m_U[1];
	    states[i][j][2] = m_cell_center[index].m_state.m_P;
	    states[i][j][3] = m_cell_center[index].m_state.m_mu;
	    states[i][j][4] = m_cell_center[index].m_state.m_rho;
	}
	{
	    printf("fname = %s\n",fname);
	    sfile = fopen(fname,"w");
	    for (j = 0; j < top_gmax[1]; j++)
	    for (i = 0; i < top_gmax[0]; i++)
	    {
		int jj = j;
		if (pp_numnodes() == 2 && pp_mynode() == 1)
			jj += 61;
	    	fprintf(sfile,"states[%d][%d] = %f %f %f %f %f\n",i,jj,
				states[i][j][0],states[i][j][1],
				states[i][j][2],states[i][j][3],
				states[i][j][4]);
	    }
	    fclose(sfile);
	}
}

void CARTESIAN::printMeshStates(char *out_name)
{
	int i,j,k,l,index;
	char filename[100];
	FILE *outfile;

	sprintf(filename,"%s-state.ts%s",out_name,right_flush(front->step,7));
#if defined(__MPI__)
        sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */
	outfile = fopen(filename,"w");
	
	fprintf(outfile,"\nInterior states:\n");
	switch (dim)
	{
	case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
	        fprintf(outfile,"%20.16f\n",
			m_cell_center[index].m_state.m_rho);
	        fprintf(outfile,"%20.16f\n",
			m_cell_center[index].m_state.m_P);
	        fprintf(outfile,"%20.16f\n",
			m_cell_center[index].m_state.m_mu);
	    	for (l = 0; l < dim; ++l)
	            fprintf(outfile,"%20.16f\n",
			    m_cell_center[index].m_state.m_U[l]);
	    }
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
	        fprintf(outfile,"%20.16f\n",
			m_cell_center[index].m_state.m_rho);
	        fprintf(outfile,"%20.16f\n",
			m_cell_center[index].m_state.m_P);
	        fprintf(outfile,"%20.16f\n",
			m_cell_center[index].m_state.m_mu);
	    	for (l = 0; l < dim; ++l)
	            fprintf(outfile,"%20.16f\n",
			    m_cell_center[index].m_state.m_U[l]);
	    }
	}
	fclose(outfile);
}

void CARTESIAN::readMeshStates(char *restart_name)
{
	FILE *infile;
	int i,j,k,l,index;

	makeGridIntfc();

	m_rho[0] = 1.0;		m_rho[1] = 0.1;
	m_mu[0]  = PI/50.0;	m_mu[1]  = PI/50.0;
	m_smoothing_radius = top_h[0] < top_h[1] ? top_h[1] : top_h[0];
	m_smoothing_radius *= 2.0;

	infile = fopen(restart_name,"r");

	next_output_line_containing_string(infile,"Interior states:");

	switch (dim)
	{
	case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
	    	fscanf(infile,"%lf",&m_cell_center[index].m_state.m_rho);
	    	fscanf(infile,"%lf",&m_cell_center[index].m_state.m_P);
	    	fscanf(infile,"%lf",&m_cell_center[index].m_state.m_mu);
		for (l = 0; l < dim; ++l)
	    	    fscanf(infile,"%lf",&m_cell_center[index].m_state.m_U[l]);
	    }
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
	    	fscanf(infile,"%lf",&m_cell_center[index].m_state.m_rho);
	    	fscanf(infile,"%lf",&m_cell_center[index].m_state.m_P);
	    	fscanf(infile,"%lf",&m_cell_center[index].m_state.m_mu);
		for (l = 0; l < dim; ++l)
	    	    fscanf(infile,"%lf",&m_cell_center[index].m_state.m_U[l]);
	    }
	}
	fclose(infile);
	deleteGridIntfc();
}
