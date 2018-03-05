/*
 * beta.c
 * Sanity check beta transform vs. downshift.
 *
 * March 2018
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double downshift(double x, double K)
{
	K *= 2.0;
	if (0.5 <= x)
	{
		return K * (x - 0.5);
	}
	return K*x;
}

double beta_xform(double x, double beta)
{
	double prod = x * beta;
	return prod - floor(prod);
}

int main (int argc, char* argv[])
{
	double K = atof(argv[1]);

	double mid = 0.5;
	double par = 1.0;
	for (int i=0; i< 20; i++)
	{
		printf("%d mid=%g parry=%g diff = %g\n", i, mid, 0.5*par, mid-0.5*par);
		mid = downshift(mid, K);
		par = beta_xform(par, 2.0*K);
	}
}
