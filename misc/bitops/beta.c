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

double skippy(double beta)
{
	double K = 0.5*beta;
	double mid = K;
	double par = 1.0;
	par = beta_xform(par, beta);
	double bpoly = 0.0;
	double kpoly = 0.0;
	double bn = 1.0;
	double kn = 1.0;
	for (int i=0; i<50; i++)
	{
		double skip = 2.0*mid-par;

		bpoly += bn * skip;
		kpoly += kn * skip;

		bn /= beta;
		kn *= K;

		mid = downshift(mid, K);
		par = beta_xform(par, beta);
	}
}

int main (int argc, char* argv[])
{
#ifdef VERIFY_MIDPOINTS
	double K = atof(argv[1]);
	double beta = 2.0*K;

	double mid = K;
	double par = 1.0;
	par = beta_xform(par, beta);
	for (int i=0; i< 30; i++)
	{
		double mone = 2.0*mid - floor(2.0*mid);
		printf("%d mid=%g parry=%g diff = %g\n", i, mone, par, 2*mid-par);
		mid = downshift(mid, K);
		par = beta_xform(par, beta);
	}
#endif
}
