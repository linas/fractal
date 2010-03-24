
/*
 * Attempt a linear decomposition of continuous/differentiable
 * eigenvalues into fractals.
 *
 * March 2010
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"

double triangle(double x)
{
	x -= floor(x);
	x *= 2.0;
	if (x < 1.0) return x;
	return 2.0-x;
}

double takagi(double w, double x)
{
	int i;
	double tn = 1.0;
	double wn = 1.0;
	double acc = 0.0;

	for (i=0; i< 50; i++)
	{
		acc += wn * triangle(tn*x);
		tn *= 2.0;
		wn *= w;
	}

	return acc;
}

// sigma is the reversed-bit summation 
double sigma(int l, double x)
{
	double tlp1 = 2*l+1;
	double lo = floor(tlp1*x);
	int mn1 = 4*l+1 - ((int)lo); 

	double acc = 0.0;
	double tn = 1.0;
	while (mn1)
	{
		int bit = mn1 & 0x1;
		if (bit) acc += tn;
		tn *= 0.5;
		mn1 >>= 1;
	}
	return acc;
}

// eigenvector of the dyadic sawtooth
double vect(double w, int l, double x)
{
	double cee = w + (1.0-w) * sigma(l, x);

	double tlp1 = 2*l+1;
	double eig = -0.5*cee + takagi (w, tlp1*x);

	return eig;
}


double decomp(double w, double x)
{
	int l;

	double acc = 0.0;
	for (l=0; l<2000; l++)
	{
		double al = 1.0 / ((double) 2*l+1);
		// double al = log(l) / ((double) l);
		// double al = 1.0 / (log(2*l+1) * ((double) l));
		// if (l%2 == 0) al = -al;
		acc += al * vect(w, l, x);
	}

	return acc;
}

int main(int argc, char * argv[])
{

	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s <w>\n", argv[0]);
		exit(1);
	}
	double lambda = atof(argv[1]);

	int i;
	int nsteps = 600;
	double delta = 1.0 / ((double) nsteps);
#if 1
	int  l = lambda;
	for (i=0; i<nsteps; i++)
	{
		double x = ((double) i) * delta;
		double y = sigma(l, x);
		printf ("%d	%g	%g	%g\n", i, x, y, y);
	}
exit(0);
#endif

	double w = 2.0*lambda/(lambda-1.0);

	ContinuedFraction f;

	double yprev = 0.0;
	for (i=0; i<nsteps; i++)
	{
		double x = ((double) i) * delta;

		f.SetReal(x);
		double qx = f.ToFarey();
		double y = decomp(w, qx);

		double s = (y-yprev) / delta;
		yprev = y;

		printf("%d	%g	%g	%g\n", i, x, y, s);
	}

	return 0;
}
