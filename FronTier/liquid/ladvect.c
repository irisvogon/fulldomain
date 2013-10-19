/********************************************************************************
 * 		ADVECTION.c
 ********************************************************************************/
#include "liquid.h"

// public functions

/*ADVECTION::ADVECTION()
{
}*/

ADVECTION::~ADVECTION()
{
	free_these(8,m_x,m_y,m_q,m_ux,m_uy,m_vx,m_vy,m_q1d);
	free_these(8,m_u1d,m_v1d,m_f,m_g,m_delta,m_deladd,m_fadd,m_gadd);
}

void ADVECTION::initRectGrid(
	double *L,
	double *U,
	double *h,
	int *gmax)
{
	m_mitot = gmax[0] + 1;
	m_mjtot = gmax[1] + 1;
	m_mtot = getMax(m_mitot, m_mjtot);

	uni_array(&m_x,m_mtot,sizeof(double));
	uni_array(&m_y,m_mtot,sizeof(double));
	bi_array(&m_q,m_mitot,m_mjtot,sizeof(double));
	bi_array(&m_ux,m_mitot,m_mjtot,sizeof(double));
	bi_array(&m_vx,m_mitot,m_mjtot,sizeof(double));
	bi_array(&m_uy,m_mitot,m_mjtot,sizeof(double));
	bi_array(&m_vy,m_mitot,m_mjtot,sizeof(double));

	uni_array(&m_q1d,m_mtot,sizeof(double));
	uni_array(&m_u1d,m_mtot,sizeof(double));
	uni_array(&m_v1d,m_mtot,sizeof(double));

	bi_array(&m_f,m_mitot,m_mjtot,sizeof(double));
	bi_array(&m_g,m_mitot,m_mjtot,sizeof(double));
	bi_array(&m_delta,m_mitot,m_mjtot,sizeof(double));

	uni_array(&m_deladd,m_mtot,sizeof(double));
	uni_array(&m_fadd,m_mtot,sizeof(double));
	bi_array(&m_gadd,m_mtot,2,sizeof(double));

	// define grid
	for(int i=0; i<m_mitot; i++)		
	    m_x[i] = L[0] + i*h[0];
	for(int i=0; i<m_mjtot; i++)		
	    m_y[i] = L[1] + i*h[1];
	
	m_method = METHOD_LaxWendroff;	// unstable for 160x160?

	m_limiter = LIMITER_NONE;
}

// testing
void ADVECTION::advection(void)
{
	m_T = 1;

	// set size of computational rectangle and boundary conditions
	// depending on test problem
	
	int ibc = 1;
	m_dt = 0.4 * getMin(m_dh[0], m_dh[1]);
	m_dt = m_T / int(m_T/m_dt);
	
	// expand domain by 2 points in each direction for BCs

	// initialize q and set velocities
	init(m_q,m_ux,m_vx,m_uy,m_vy,m_dh[0],m_mitot,m_mjtot);

	// extend initial data to bordering cells
	//bc(m_q,m_ux,m_vx,m_uy,m_vy,m_dh[0],m_mitot,m_mjtot,ibc,0);
	bc(m_q, ibc, 0);


	double tnph;					
	double tnp;						
	double t;
	bool bPrinted = false;

	m_method = METHOD_LaxWendroff;	// unstable for 160x160?

	m_limiter = LIMITER_VANLEER;
	//m_limiter = LIMITER_MC;

	t = m_dt;
	while(t<=m_T+0.0001)
	{
		printf(" t = %f \n", t);
		tnph = t - 0.5 * m_dt;		// time n+1/2
		tnp  = t;					// time n+1
		step(m_q, tnph);
		bc(m_q, ibc, tnp);
		t += m_dt;
	}

	char filename[200];
	sprintf(filename, "out/cartesian_simple_advection_printError.plt");
	printError(filename, m_q, m_mitot, m_mjtot);
	printMaxError(m_q, m_mitot, m_mjtot);
	printf("m_T = %f \n", m_T);
}

// private functions
double ADVECTION::getExactSolution(double coords[2])
{
	double x,y;
	x = coords[0];
	y = coords[1];
	return sin(2*PI*x) * sin(2*PI*y);
}

void ADVECTION::step(
		double **q, 
		double t)
{
	int i, j;

	// take one time step
	
	// initialize fluxes
	for(j=0; j<m_mjtot; j++)
	for(i=0; i<m_mitot; i++)
	{
		m_f[i][j] = 0;
		m_g[i][j] = 0;
		m_delta[i][j] = 0;
	}

	double tfactor;
	// perform x-sweeps
	tfactor = gt(t);
	for(j=1; j<m_mjtot-1; j++)
	{
		for(i=0; i<m_mitot; i++)
		{
			m_q1d[i] = q[i][j];	
			m_u1d[i] = m_ux[i][j] * tfactor;
			m_v1d[i] = m_vx[i][j] * tfactor;
		}
		fluxes(m_q1d,m_u1d,m_v1d,m_deladd,m_fadd,m_gadd,0,m_mitot,m_dh[0],m_dt);
		for(i=2; i<m_mitot-1; i++)
		{
			m_delta[i][j] 	+= m_deladd[i];
			m_f[i][j] 		+= m_fadd[i];
			m_g[i][j]		+= m_gadd[i][0];
			m_g[i][j+1]		+= m_gadd[i][1];
		}
	}
	// perform y-sweeps
	// note: u & v are switched and so are fadd & gadd
	for(i=1; i<m_mitot-1; i++)
	{
		for(j=0; j<m_mjtot; j++)
		{
			m_q1d[j] = m_q[i][j];
			m_u1d[j] = m_vy[i][j] * tfactor;		
			m_v1d[j] = m_uy[i][j] * tfactor;
		}
		fluxes(m_q1d,m_u1d,m_v1d,m_deladd,m_fadd,m_gadd,0,m_mjtot,m_dh[0],m_dt);
		for(j=2; j<m_mjtot-1; j++)
		{
			m_delta[i][j] 	+= m_deladd[j];
			m_g[i][j]		+= m_fadd[j];
			m_f[i][j]		+= m_gadd[j][0];
			m_f[i+1][j] 	+= m_gadd[j][1];
		}
	}
	// update q values by differencing fluxes
	for(j=2; j<m_mjtot-2; j++)
	for(i=2; i<m_mitot-2; i++)
		q[i][j] += -m_dt/m_dh[0] * (  m_f[i+1][j] - m_f[i][j] + m_g[i][j+1] - m_g[i][j]
							 + m_delta[i][j] );
}

void ADVECTION::fluxes(
		double *q, 
		double *u, 
		double *v, 
		double *deladd, 
		double *fadd, 
		double **gadd, 
		int i1, 
		int i2, 
		double dx, 
		double dt)
{
	int i, j;
	double *dq, fincr, gincr, dq2, fincr3, fract;

	uni_array(&dq,i2-i1,sizeof(double));

	// initialize fluxes
	for(i=i1; i<i2; i++)
	{
		deladd[i] = 0;
		fadd[i] = 0;
		for(int jside=0; jside<2; jside++)
			gadd[i][jside] = 0;
	}
	for(i=i1+1; i<i2; i++)
	{
		dq[i] = q[i] - q[i-1];
		if(u[i]>0)
			j = i-1;			// different from the article???
		else
			j = i;

		fadd[i] += u[i] * q[j];

		// basic first order increment for this cell based on piecewise
		// constant dq propagating normal to interface
		//deladd[j] += u[i]*dq[i];		
		if(m_method==METHOD_DonorCell) continue;	

		if(m_method==METHOD_LaxWendroff)
		{
			gincr = 0.125 * dt/dx * u[i]*v[i]*dq[i];
			gadd[i][0] 	+= - gincr;
			gadd[i][1]	+= - gincr;
			gadd[i-1][0]+= - gincr;
			gadd[i-1][1]+= - gincr;
			continue;
		}

		if(u[i]>0)
			j = i;			// different from the article???
		else
			j = i-1;

		// corrections for transverse propagation of dq
		if(v[i]>0)
			gadd[j][1] += - 0.5 * dt/dx * u[i]*v[i]*dq[i];	// upwind
		else
			gadd[j][0] += - 0.5 * dt/dx * u[i]*v[i]*dq[i];
	}


	if(m_method<=2)
		return;

	// modify F fluxes for second order uxx terms
	limiter(dq,u,i1,i2);
	for(i=i1+2; i<i2-1; i++)
	{
		fincr = 0.5 * fabs(u[i]) * (1 - fabs(u[i])*dt/dx) *dq[i];
		fadd[i] += fincr;
		if(m_method==METHOD_SecondOrderCorrectionTermsAdded)
			continue;	
		if(m_method==METHOD_LaxWendroff)
			continue;
		if(m_method>=5)			// third order correction
		{
			if(u[i]>0)
				dq2 = dq[i] - dq[i-1];
			else
				dq2 = dq[i+1] - dq[i];
			fincr3 = - u[i]/6 * (1-getSqr(u[i]*dt/dx)) * dq2;
			fadd[i] += fincr3;
			if(m_method==METHOD_TransversePropagation3)
				fincr += fincr3;
		}

		// compute fraction of dq that goesing into adjacent row
		// (transverse propagation of second order and possibly
		// third order correction waves)
		fract = fabs(v[i])*dt/dx *fincr;
		if(v[i]<0)
		{
			gadd[i][0] 	+= - fract;
			gadd[i-1][0]+= + fract;
		}
		else
		{
			gadd[i][1]	+= + fract;
			gadd[i-1][1]+= - fract;
		}
	}

	free_these(1,dq);
}
void ADVECTION::limiter(
		double *dq, 
		double *u, 
		int i1,
		int i2)
{
	if(m_limiter == LIMITER_NONE)
		return;
	double dqm1, phi;

	dqm1 = dq[i1+1];
	for(int i=i1+2; i<i2-1; i++)
	{
		if(fabs(dq[i])<0.0000001)		// changed a little bit
			continue;
		if(u[i]>0)
			phi = philim(dq[i], dqm1, m_limiter);
		else
			phi = philim(dq[i], dq[i+1], m_limiter);
		dqm1 = dq[i];
		dq[i] = dq[i] * phi;
	}
}

double ADVECTION::philim(double a, double b, LIMITER_TYPE meth)
{
	double r;
	r = b/a;
	switch(meth)
	{
	case LIMITER_MINMOD:
		return getMax(0.,getMin(1.,r));
		break;
	case LIMITER_SUPERBEE:
		return getMax(0.,getMin(1.,2*r),getMin(2.,r));
		break;
	case LIMITER_VANLEER:
		return (r+fabs(r))/(1+fabs(r));
		break;
	case LIMITER_MC:
		return getMax(0.,getMin((1+r)/2,2.,2*r));
		break;
	default:
		printf("ADVECTION::philim: shouldn't be here.\n");
	}
}

void ADVECTION::bc(
		double **q, 
		int ibc, 
		double tn)
{
	switch(ibc)
	{
	case 1:								// BC_PERIODIC:
		for(int i=0; i<m_mitot; i++)
		{
			q[i][0] = q[i][m_mjtot-4];
			if(m_mjtot==5) 
				q[i][0] = q[i][2];		// in case only one pt in y
			q[i][1] = q[i][m_mjtot-3];
			q[i][m_mjtot-2] = q[i][2];
			q[i][m_mjtot-1] = q[i][3];
		}
		for(int j=0; j<m_mjtot; j++)
		{
			q[0][j] = q[m_mitot-4][j];
			if(m_mitot==5)
				q[0][j] = q[2][j];		// in case only one pt in x
			q[1][j] = q[m_mitot-3][j];
			q[m_mitot-2][j] = q[2][j];
			q[m_mitot-1][j] = q[3][j];
		}
		break;
	case 2:								// BC_Extrapolation0:
		for(int i=0; i<m_mitot; i++)
		{
			q[i][0] = q[i][2];
			q[i][1] = q[i][2];
			q[i][m_mjtot-2] = q[i][m_mjtot-3];
			q[i][m_mjtot-1] = q[i][m_mjtot-3];
		}
		for(int j=0; j<m_mjtot; j++)
		{
			q[0][j] = q[2][j];
			q[1][j] = q[2][j];
			q[m_mitot-2][j] = q[m_mitot-3][j];
			q[m_mitot-1][j] = q[m_mitot-3][j];
		}
		break;
	case 3:								// BC_Extrapolation1:
		for(int i=0; i<m_mitot; i++)
		{
			q[i][1] = 2 * q[i][2] - q[i][3];
			q[i][0] = q[i][1];
			q[i][m_mjtot-2] = 2 * q[i][m_mjtot-3] - q[i][m_mjtot-4];
			q[i][m_mjtot-1] = q[i][m_mjtot-2];
		}
		for(int j=0; j<m_mjtot; j++)
		{
			q[1][j] = 2 * q[2][j] - q[3][j];
			q[0][j] = q[1][j];
			q[m_mitot-2][j] = 2 * q[m_mitot-3][j] - q[m_mitot-4][j];
			q[m_mitot-1][j] = q[m_mitot-2][j];
		}
		break;
	case 4:				// BC_InflowOutflow:
		break;
	}

}

void ADVECTION::init(
		double **q,
		double **ux,
		double **vx,
		double **uy,
		double **vy,
		double dx,
		int mitot,
		int mjtot)
{
	// swirling deformation flow
	double u0 = 1, xl, yb;
	for(int i=0; i<mitot; i++)
	{
		xl = m_x[i] - 0.5 * dx;
		for(int j=0; j<mjtot; j++)
		{
			yb = m_y[j] - 0.5 * dx;
			uy[i][j] = 1;
			vy[i][j] = 2;
			ux[i][j] = 1;
			vx[i][j] = 2;
		}
	}

	// set concentration profile
	double coords[2];
	for(int i=0; i<mitot; i++)
	for(int j=0; j<mjtot; j++)
	{
		coords[0] = m_x[i];
		coords[1] = m_y[j];
		q[i][j] = getExactSolution(coords);
	}
}

double ADVECTION::gt(double t)
{
	//return cos(PI*t/m_T);
	return 1;
}

double ADVECTION::getMin(double a, double b)
{
	if(a<=b)
		return a;
	else
		return b;
}
double ADVECTION::getMin(double a, double b, double c)
{
	return getMin(a, getMin(b,c));
}
int ADVECTION::getMax(int a, int b)
{
	if(a<=b)
		return b;
	else
		return a;
}
double ADVECTION::getMax(double a, double b)
{
	if(a<=b)
		return b;
	else
		return a;
}

double ADVECTION::getMax(double a, double b, double c)
{
	return getMax(a, getMax(b, c));
}

double ADVECTION::getSqr(double a)
{
	return a*a;
}

void ADVECTION::printError(char *filename, double **q, int mitot, int mjtot)
{
	printf("ADVECTION::printError: writing %s \n", filename);
	FILE *hfile = fopen(filename, "w");
	if(hfile==NULL)
	{
	    printf("\n ADVECTION::printError: "
	    		"can't open %s for writing.", filename);
	    exit(0);
	}
	fprintf(hfile, "TITLE = %s \n", filename);
	fprintf(hfile, "VARIABLES = X Y P_exact P_approx P_error \n ");
	fprintf(hfile, "ZONE I=%d J=%d F=point \n", mitot, mjtot);

	double coords[2];
	for(int i=0; i<mitot; i++)
	for(int j=0; j<mjtot; j++)
	{
		coords[0] = m_x[i];
		coords[1] = m_y[j];
		fprintf(hfile, "%f %f %f %f %f \n", coords[0], coords[1],
				getExactSolution(coords),
				q[i][j],
				getExactSolution(coords) - q[i][j]);
	}
	fclose(hfile);
}

void ADVECTION::printMaxError(double **q, int mitot, int mjtot)
{
	double coords[2];
	double max_error = 0, tmp;
	for(int i=2; i<mitot-2; i++)
	for(int j=2; j<mjtot-2; j++)
	{
		coords[0] = m_x[i];
		coords[1] = m_y[j];
		tmp = fabs(getExactSolution(coords) - q[i][j]);
		if(tmp>max_error)
			max_error = tmp;
	}
	printf("ADVECTION::printMaxError: max_error = %f \n", max_error);
}

