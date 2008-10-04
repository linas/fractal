
/*
 * Build a takaagi-style distribution
 *
 * Linas Vepstas October 2008
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"

double triangle (double x)
{
	x -= (int) x;
	if (x < 0.5) return x;
	return 1.0-x;
}

double olde_a_k (double x, int k)
{
	int i;
	for (i=0; i<k; i++)
	{
		if (x < 0.5) x = x/(1.0-x);
		else x = (2.0*x-1)/x; 
		// else x = (1.0-x)/x; 
	}
	return triangle (x);
}

double almost_dedekind_zeta (double x)
{
	int n;
	double prod = 1.0;
	double tk = 1.0;
	for (n=0; n< 20; n++)
	{
		double mand = 1.0 + olde_a_k(x,n+1);
		prod /= 0.5 * mand*mand;
	}

	return prod;
}

double a_k (double x, int k)
{
	int i;
	for (i=0; i<k; i++)
	{
		if (x < 0.5) x = x/(1.0-x);
		else x = (2.0*x-1)/x; 
		// else x = (1.0-x)/x; 
	}
	return x;
}

double eff(double x, double lm)
{
	double f;
	if (x < 0.5) f = 1.0-x;
	else f = x;
	f = 1.0/f;
   f *= f;
	f = exp (lm * log(f));
	f *=0.5;
	return f;
}

double prod (double x, int depth, double lambda)
{
	int n;
	double prod = 1.0;
	double tk = 1.0;
	double lmp = 1.0;
	for (n=0; n< depth; n++)
	{
		// tk *= 2.0;
		// double mand = 1.0 + olde_a_k (x, n+1);
		// double mand = 1.0 + triangle(tk*x);
		// double mand =  1.25 + 0.25*sin(M_PI*tk*x);
		// double mand = 1.0/eff(a_k(x,n), lambda);
		double mand = a_k(x,n);
		prod *= eff(mand, lmp);
		lmp *= lambda;
	}

	return prod;
}

void graph(int npts, int depth, double lambda)
{
	int i;

   ContinuedFraction f;

	/* Compute the integral of the distribution */
	double acc = 0.0;
	double delta = 1.0 / (double (npts));
	for (i=0; i<npts; i++)
	{
		double x = (double) i / ((double) npts);
		double y = prod(x, depth, lambda);
		acc += y*delta;

   	f.SetRatio(i, npts);
   	double far = f.ToFarey ();

		printf ("%d	%g	%g %g	%g	%g\n", i, x, y, acc, far, acc-far);
		fflush (stdout);
	}
}


main(int argc, char *argv[])
{
	int i;

	if (argc <4)
	{
		fprintf (stderr, "Usage: %s <nbins> <depth> <lambda>\n", argv[0]);
		exit (1);
	}
	int nbins = atoi (argv[1]);
	int depth = atoi (argv[2]);
	double lambda = atof (argv[3]);

	graph (nbins, depth, lambda);
}
