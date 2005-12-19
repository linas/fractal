/*
 * Semi-Symbolic polynomial handling in c language
 * Very simple, only handles polynomials in one variable
 *
 * Linas Vepstas 16 December 2005
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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
	double mean;
	Poly derivs[MAXTERMS];
} Gaussian;

Gaussian *
gaussian_new (double mu, double sigma)
{
	Gaussian *g;
	g = (Gaussian *) malloc (sizeof (Gaussian));

	g->sigma = sigma;
	g->mean = mu;

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
		poly_add (&g->derivs[i], &d, &p);
	}

	return g;
}

/* return the 'order'th derivative of a gaussian, at position val */
double 
gaussian_eval (Gaussian *g, int order, double val)
{
	double gv;
	gv = val - g->mean;
	gv *= gv;
	gv /= 2.0 * g->sigma * g->sigma;

	gv = exp (-gv);
	gv *= poly_eval (&g->derivs[order], val);

	return gv;
}


main () 
{
	Gaussian *g;

	g = gaussian_new (0.0, 1.0);

	double x;
	for (x=0.0; x< 3.0; x+=0.2)
	{
		double n0 = gaussian_eval (g, 0, x);
		double n1 = gaussian_eval (g, 1, x);
		double n2 = gaussian_eval (g, 2, x);
		double n3 = gaussian_eval (g, 3, x);
		printf ("%g	%g	%g	%g	%g\n", x, n0, n1, n2, n3);
	}

}
