
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

double eff(double x, int p)
{
	double f;
	if (x < 0.5) f = 1.0-x;
	else f = x;
	f = 1.0/f;

	double g = 1.0;
	for (int i=0; i<p; i++) g *= f;
	return g;
}

double enn(double x)
{
	return 0.5*eff(x, 2) - eff(x,3)/3.0;
}

double sum (double x, int depth, double lambda)
{
	int n;
	double acc = 0.0;
	double lmp = 1.0;
	for (n=0; n< depth; n++)
	{
		double mand = a_k(x,n);
		acc += lmp * enn(mand);
		lmp *= lambda;
	}

	return acc;
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
		double y = sum(x, depth, lambda);
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
