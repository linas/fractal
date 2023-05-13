/*
 * Fractal eigenfunctions of GKW.
 *
 * Used to illustrate the gkw.lyx paper.
 * Linas January 2016
 */

#include <math.h>
#include <stdio.h>

double h(double x)
{
	double y = 1.0 / x;
	return y - floor(y);
}

double rade(int n, int b, double x)
{
	for (int i=1; i<n; i++)
	{
		x = h(x);
	}
	x = 1.0 / x;
	if (b < x && x <= (b+1.0)) return 1.0;
	return 0.0;
}

double fractal(double lambda, int b, double x)
{
	double sum = 0.0;
	double ln = 1.0;
	int n = 1;
	while (1.0e-14 < ln)
	{
		sum += ln * rade(n, b, x);
		ln *= lambda;
		n++;
	}

	return sum;
}

int main(int argc, char* argv[])
{
	double x;

#ifdef RADEMACHER_GRAPH
	/* Rademacher functions for Baire space */
	/* Used to generate the graphs baire-13.eps for the gkw paper */
	for (x=0.0; x<1.0; x+=0.0001)
	{
		double r1 = rade(1,3,x);
		double r2 = rade(2,1,x);
		double r3 = rade(3,1,x);
		printf("%f	%f	%f	%f\n", x, r1, r2, r3);
	}
#endif

	double lambda = 0.5;
	double bee = 1;
	for (x=0.0; x<1.0; x+=0.00001)
	{
		double r1 = rade(1,bee,x);
		double fra = fractal(lambda, bee, x);
		double r2 = rade(1,2,x);
		double fra2 = fractal(lambda, 2, x);
		printf("%f	%f	%f	%f	%f\n", x, r1, fra, r2, fra2);
	}

	return 0;
}

