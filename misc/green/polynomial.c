/*
 * Semi-Symbolic polynomial handling in c language
 * Very simple, only handles polynomials in one variable
 *
 * Linas Vepstas 16 December 2005
 */

#include <math.h>

#define MAXORDER 100

typedef double Poly[MAXORDER];

void 
poly_clear (Poly *b)
{
	int i;
	for (i=0; i<MAXORDER; i++)
	{
		(*b)[i] = 0.0;
	}
}

void 
poly_add (Poly *sum, Poly *a, Poly *b)
{
	int i;
	for (i=0; i<MAXORDER; i++)
	{
		(*sum)[i] = (*a)[i] + (*b)[i];
	}
}

void 
poly_mult (Poly *prod, Poly* a, Poly *b)
{
	int i,j;
	
	for (i=0; i<MAXORDER; i++)
	{
		double acc = 0.0;
		for (j=0; j<=i; j++)
		{
			acc += (*a)[j] * (*b)[i-j];
		}
		(*prod)[i] = acc;
	}
}

void 
poly_differntiate (Poly *deriv, Poly *in)
{
	int i;
	for (i=0; i<MAXORDER-1; i++)
	{
		(*deriv)[i] = (i+1) * (*in)[i+1];
	}
}

double 
poly_eval (Poly *a, double val)
{
	int i;
	double acc = 0.0;
	double term = 1.0;
	for (i=0; i<MAXORDER; i++)
	{
		acc += term * (*a)[i];
		term *= val;
	}

	return acc;
}

#define MAXTERMS 100

/* derivatives of the standard normal distribution 
 * with mean mu and standard deviation sigma 
 */
typedef struct 
{
	double sigma;
	double mu;
	Poly derivs[MAXTERMS];
} Gaussian;

void 
setup_gaussians (Gaussian *g, double sigma, double mu)
{
	g->sigma = sigma;
	g->mu = mu;

	int i;
	for(i=0; i<MAXTERMS; i++)
	{
		poly_clear (&(g->derivs)[i]);
	}
	g->derivs[0][0] = 1.0;
	g->derivs[1][0] = mu/(sigma*sigma);
	g->derivs[1][1] = -mu/(sigma*sigma);
	
	for(i=2; i<MAXTERMS; i++)
	{
		Poly d,p;
		poly_clear (&d);
		poly_differntiate (&d,  &g->derivs[i-1]);
		poly_clear (&p);
		poly_mult (&p, &g->derivs[i-1], &g->derivs[1]);
		poly_sum (&g->derivs[i], &d, &p);
	}
}

double 
eval_gaussian (int order, double val)
{
	double g;

}

