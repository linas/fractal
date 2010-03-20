
/*
 * Attempt a linear decomposition of continuous/differentiable
 * eigenvalues into fractals.
 *
 * March 2010
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double triangle(double x)
{
	x -= floor(x);
	x *= 2.0;
	if (x < 1.0) return x;
	return 2.0-x;
}

double tagaki(double w, double x)
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

double decomp(double w, double x)
{
	int l;

	double acc = 0.0;
	for (l=1; l<2000; l++)
	{
		double al = 1.0 / ((double) l);
		double tlp1 = 2*l+1;
		acc += al * takagi(w, tlp1 * x);
	}

	return acc;
}

int main(int agrc, char * argv[])
{

	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s <w>\n", argv[0]);
		exit(1);
	}
	double w = atof(argv[1]);

	int nsteps = 600;
	int i;
	for (i=0; i<nsteps)
	{
		double x = ((double) i) / ((double) nsteps);

		double y = decomp(x);

		printf("%d	%g	%g\n", i, x, y);
	}

	return 0;
}
