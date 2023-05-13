/*
 * Consider a product of takagi distributions
 *
 * XXX this is a failed experiment.
 *
 * Linas Vepstas March 2010
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// #include "Farey.h"

double triangle (double x)
{
	x -= (int) x;
	if (x < 0.5) return 2.0*x;
	return 2.0-2.0*x;
}

double takagi (double x, double lambda)
{
	int n;
	double acc = 0.0;
	double lmp = 1.0;
	double tn = 1.0;
	for (n=0; n<25; n++)
	{
		double term = triangle(tn*x);
		acc += lmp * term;
		lmp *= lambda;
		tn *= 2.0;
	}

	return acc;
}
double prod (double x, int depth, double lambda)
{
	int n;
	double acc = 1.0;
	double tlp1 = 1.0;
	double lmp = lambda;
	for (n=0; n<depth; n++)
	{
		double term = takagi(tlp1*x, lmp);
		acc *= term;
		tlp1 += 2.0;
		lmp *= lambda;
	}

	return acc;
}

void graph(int npts, int depth, double lambda)
{
	int i;

	// ContinuedFraction f;

	/* Compute the integral of the distribution */
	double acc = 0.0;
	double delta = 1.0 / (double (npts));
	for (i=0; i<npts; i++)
	{
		double x = (double) i / ((double) npts);
		double y = prod(x, depth, lambda);
		acc += y*delta;

   	// f.SetRatio(i, npts);
   	// double far = f.ToFarey ();

		// printf ("%d	%g	%g %g	%g	%g\n", i, x, y, acc, far, acc-far);
		printf ("%d	%10.6g	%10.6g	%10.6g\n", i, x, y, acc);
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
