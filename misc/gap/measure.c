/*
 * measure.c
 *
 * This produces at least some of the *discontinuous* eigenfunctions
 * of the Bernoulli operator. I've probably written this code somewhere
 * before, but don't knw where.
 *
 * The aproach taken is that of working with the measure-theoretic
 * definition of the thing, using the product topology.
 *
 * Linas Oct 2008
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double sum(double lambda, unsigned int i, int dim)
{
	double lp = 1.0;
	int k=0;

	unsigned int filt = 1<<(dim-1);

	double sum = 0.0;
	for (k=0; k<dim; k++)
	{
		if (i & filt)
		{
			sum += lp;
		}
		else
		{
			sum -= lp;
		}
		filt >>= 1;
		lp *= lambda;
	}

	sum *= (1.0-lambda);
	return sum;
}

void beig(int np, double lambda)
{
	int npts = 1<<np;
	int i;

	for(i=0; i<npts; i++)
	{
		double x = ((double) i)/((double) npts);
		double y = sum(lambda, i, np);
		printf("%g	%g\n", x, y);
	}
}

int main(int argc, char * argv[])
{

	double lambda = atof(argv[1]);

	printf ("#\n# lambda=%g\n#\n", lambda);

	beig(14, lambda);
}
