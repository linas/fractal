/*
 * measure.c
 *
 * This produces at least some of the *discontinuous* eigenfunctions
 * of the Bernoulli operator. I've probably written this code somewhere
 * before.
 *
 * Linas Oct 2008
 */
#include <stdio.h>

double sum(double lambda, int i, int dim)
{
	double lp = 1.0;
	int k=0;

	double sum = 0.0;
	for (k=0; k<dim; k++)
	{
		if (i & 0x1)
		{
			sum += 1.0/lp;
		}
		else
		{
			sum -= 1.0/lp;
		}
		i >>= 1;
		lp *= lambda;
	}

	sum *= lp;
	sum /= lambda;
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

main()
{
	beig(10, 0.35);

}
