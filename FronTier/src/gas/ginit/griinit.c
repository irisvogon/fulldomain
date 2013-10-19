#include <ginit/ginit.h>
#include <complex.h>

//complex number 
//http://en.wikipedia.org/wiki/Complex.h

//the functions come from 
//http://www.fizyka.umk.pl/nrbook/bookcpdf.html
//numerical recipes in C

//Returns the modified Bessel function I0(x) for any real x.
float bessi0(float x)
{
	float	ax,ans;
	float	y;
	
	if((ax=fabs(x)) < 3.75) 
	{ 
	    y = x/3.75;
	    y *= y;
	    ans = 1.0 + y*(3.5156229 + y*(3.0899424 + y*(1.2067492
		+ y*(0.2659732 + y*(0.360768e-1 + y*0.45813e-2)))));
	} 
	else 
	{
	    y = 3.75/ax;
	    ans = (exp(ax)/sqrt(ax)) * (0.39894228 + y*(0.1328592e-1
		+ y*(0.225319e-2 + y*(-0.157565e-2 + y*(0.916281e-2
		+ y*(-0.2057706e-1 + y*(0.2635537e-1 + y*(-0.1647633e-1
		+ y*0.392377e-2))))))));
	}

	return ans;
}

//Returns the modified Bessel function K0(x) for positive real x.
float bessk0(float x)
{
	float	y,ans;

	if(x <= 2.0) 
	{
	    y=x*x/4.0;
	    ans=(-log(x/2.0)*bessi0(x)) + (-0.57721566 + y*(0.42278420
		+ y*(0.23069756 + y*(0.3488590e-1 + y*(0.262698e-2
		+ y*(0.10750e-3 + y*0.74e-5))))));
	} 
	else 
	{
	    y=2.0/x;
	    ans=(exp(-x)/sqrt(x)) * (1.25331414 + y*(-0.7832358e-1
		+ y*(0.2189568e-1 + y*(-0.1062446e-1 + y*(0.587872e-2
		+ y*(-0.251540e-2 + y*0.53208e-3))))));
	}
	return ans;
}

//Returns the modified Bessel function I1(x) for any real x.
float bessi1(float x)
{
	float	ax,ans;
	float	y;

	if((ax=fabs(x)) < 3.75) 
	{
	    y = x/3.75;
	    y *= y;
	    ans = ax*(0.5 + y*(0.87890594 + y*(0.51498869 + y*(0.15084934
		  + y*(0.2658733e-1 + y*(0.301532e-2 + y*0.32411e-3))))));
	} 
	else 
	{
	    y=3.75/ax;
	    ans = 0.2282967e-1 + y*(-0.2895312e-1 + y*(0.1787654e-1
		  - y*0.420059e-2));
	    ans = 0.39894228 + y*(-0.3988024e-1 + y*(-0.362018e-2
		  + y*(0.163801e-2 + y*(-0.1031555e-1 + y*ans))));
	    ans *= (exp(ax)/sqrt(ax));
	}
	
	return x < 0.0 ? -ans : ans;
}

//Returns the modified Bessel function K1(x) for positive real x.
float bessk1(float x)
{
	float	y,ans; 

	if(x <= 2.0) 
	{
	    y = x*x/4.0;
	    ans = (log(x/2.0)*bessi1(x)) + (1.0/x)*(1.0 + y*(0.15443144
		  + y*(-0.67278579 + y*(-0.18156897 + y*(-0.1919402e-1
		  + y*(-0.110404e-2 + y*(-0.4686e-4)))))));
	} 
	else 
	{
	    y = 2.0/x;
	    ans = (exp(-x)/sqrt(x))*(1.25331414 + y*(0.23498619
		  + y*(-0.3655620e-1 + y*(0.1504268e-1 + y*(-0.780353e-2
		  + y*(0.325614e-2 + y*(-0.68245e-3)))))));
	}
	
	return ans;
}

int	RI_comp_func(void*,float*);
void	RI_invicid_solution_outer(float*,float*,float*,float*,
		float,float,float,float,RI_PERT*);
void	RI_invicid_solution_inner(float*,float*,float*,float*,
		float,float,float,float,RI_PERT*);
void	RI_invicid_growth_rate(RI_PERT*,float,float);

//r > a
void	RI_invicid_solution_outer(
	float		*u,
	float		*v,
	float		*p,
	float		*d,
	float		r,
	float		z,
	float		t,
	float		dens,
	RI_PERT		*ri_param)
{
	float		a, kx, beta, K0, K1, K1a;
	float complex	expon, phi, divel;

	kx = ri_param->kx;
	beta = ri_param->beta;
	a = ri_param->a;

	divel = beta + I*kx*ri_param->U;
	expon = ri_param->amp*cexp(beta*t + I*kx*z);
	K0 = bessk0(kx*r);
	K1 = bessk1(kx*r);
	K1a = bessk1(kx*a);

	phi = -1.0*divel/kx*K0/K1a*expon;
	
	*u = creal(-1.0*divel*K0/K1a*expon);
	*v = cimag(divel*K1/K1a*expon);
	*p = cimag(-1.0*dens*divel*phi);
	*d = 0.0;
}

//r < a
void	RI_invicid_solution_inner(
	float		*u,
	float		*v,
	float		*p,
	float		*d,
	float		r,
	float		z,
	float		t,
	float		dens,
	RI_PERT		*ri_param)
{
	float		a, kx, beta, I0, I1, I1a;
	float complex	expon, phi;

	kx = ri_param->kx;
	beta = ri_param->beta;
	a = ri_param->a;

	expon = ri_param->amp*cexp(beta*t + I*kx*z);
	I0 = bessi0(kx*r);
	I1 = bessi1(kx*r);
	I1a = bessi1(kx*a);

	phi = beta/kx*I0/I1a*expon;
	
	*u = creal(beta*I0/I1a*expon);
	*v = cimag(beta*I1/I1a*expon);
	*p = cimag(-dens*beta*phi);
	*d = 0.0;
}

void	RI_invicid_growth_rate(
	RI_PERT		*ri_param,
	float		dens,
	float		sigma)
{
	float		a,kx,beta, I0,I1;

	kx = ri_param->kx;
	a = ri_param->a;

	I0 = bessi0(kx*a);
	I1 = bessi1(kx*a);
	beta = sigma*kx/(dens*a*a)*(1.0-kx*kx*a*a)*I1/I0;
	if(beta < 0)
	{
	    printf("ERROR RI_invicid_growth_rate, wave is stable.\n");
	    clean_up(ERROR);
	}

	ri_param->beta = sqrt(beta);
	
	printf("#grate = %15.8e\n", ri_param->beta);
}


//assume the domain in z direction is [-lambda/2, lambda/2]
int   RI_comp_func(
        void  *func_params,
        float *coords)
{
	RI_PERT		*params = (RI_PERT*) func_params;
	float		lambda, r, z, *cen, kx;

	kx = params->kx;
	cen = params->cen;
	lambda = 2.0*PI/kx;
	r = sqrt(sqr(coords[0]-cen[0]) + sqr(coords[1]-cen[1]));
	z = coords[2] + lambda/4.0;
	
	if(r < params->a + params->amp*sin(kx*z))
	    return 2;
	return 3;
}

void	RI_test_output(
	RI_PERT		*params,
	int		numr,
	int		numz,
	float		dens,
	float		dens1)
{
	int	i, j;
	float	lambda, R, kx;
	int	comp;
	float	z0, z, r, u, v, p, d;
	FILE	*fp;

	kx = params->kx;
	lambda = 2*PI/kx;
	R = params->R;

	fp = fopen("RI.plt", "w");

	fprintf(fp, "TITLE = \" interior state \" \n");
	fprintf(fp, "VARIABLES = ");
	fprintf(fp, "\"z\", \"r\", \"comp\", \"u\", \"v\",\"p\" \n");
	fprintf(fp,"ZONE i=%d, j=%d\n", numz, numr);

	for(i=0; i<numr; i++)
	    for(j=0; j<numz; j++)
	    {
		r = i*R/(numr-1);
		z0 = j*lambda/(numz-1) - lambda/2.0;

		z = z0 + lambda/4.0;

		if(r < params->a + params->amp*sin(kx*z))
		{
		    comp = 2;
		    RI_invicid_solution_inner(&u,&v,&p,&d,r,z,0.0,dens,params);
		}
		else
		{
		    comp = 3;
		    RI_invicid_solution_outer(&u,&v,&p,&d,r,z,0.0,dens1,params);
		}

		fprintf(fp, "%15.8e %15.8e  %d  %15.8e %15.8e %15.8e\n",
			z0, r, comp, u, v, p);
	    }

	fclose(fp);
}


int	RI_main_test()
{
	RI_PERT		ri_param;
	float		dens,dens1,sigma,a;

	a = 0.005;
	dens = 0.66;
	dens1 = 0.0066;
	sigma = 0.01*a;

	ri_param.a = a;
	ri_param.kx = 2.0*PI/(10.0*a);
	ri_param.amp = 0.05*a;
	ri_param.R = 3.0*a;
	ri_param.cen[0] = 0.0;
	ri_param.cen[1] = 0.0;
	ri_param.cen[2] = 0.0;
	ri_param.U = 0.0;
	
	RI_invicid_growth_rate(&ri_param, dens, sigma);
	RI_test_output(&ri_param, 60, 200, dens, dens1);

	return 0;
}


