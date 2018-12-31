/*
 * beta.c
 * Sanity check beta transform vs. downshift.
 * Computes the exact Parry-Gel'fond expression, comparse to
 * experiments.
 *
 * March 2018
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// My downshift; beta=2K so T(x) = case bx or b(x-0.5)
double downshift(double x, double K)
{
	K *= 2.0;
	if (0.5 <= x)
	{
		return K * (x - 0.5);
	}
	return K*x;
}

// Standard beta-xform is t(x) = (b x) mod 1
double beta_xform(double x, double beta)
{
	double prod = x * beta;
	return prod - floor(prod);
}

// This is attempting to be a Witt-vector-like thing.
// This needs further work...
double witt_xform(double x, double beta, int n)
{
	double prod = pow (beta, n);

	return x - floor(x*prod)/prod;
}

#define NBITS 850

double skippy(double beta)
{
	double K = 0.5*beta;
	double mid = K;
	double par = 1.0;
	par = beta_xform(par, beta);
	double bpoly = 0.0;
	double bn = 1.0;
	for (int i=0; i<NBITS; i++)
	{
		double skip = 2.0*mid-par;

		bpoly += bn * skip;

		bn /= beta;

		mid = downshift(mid, K);
		par = beta_xform(par, beta);
	}

	return beta - bpoly;
}

double skipry(double beta)
{
	double K = 0.5*beta;
	double mid = K;
	double par = 1.0;
	par = beta_xform(par, beta);
	double kpoly = 0.0;
	double kn = 1.0;
	for (int i=0; i<NBITS; i++)
	{
		double skip = 2.0*mid-par;

		kpoly += kn * skip;

		kn *= K;

		mid = downshift(mid, K);
		par = beta_xform(par, beta);
	}
	return kpoly;
}

int main (int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s K x\n", argv[0]);
		exit(1);
	}
	double K = atof(argv[1]);
	double x = atof(argv[2]);
	double beta = 2.0*K;

	double mid = x;
	double par = x/K;
	double bn = 1.0;
	for (int i=0; i< 30; i++)
	{
		// k's are the beta-expansion
		int kn = 0;
		if (0.5 <= mid) kn = 1;

		// en's are the invariant measure
		int en = 0;
		if (par <= x) en = 1;

		double witt = witt_xform(x, beta, i);
		witt *= bn;

		double diff = witt-mid;

		printf("%d   down=%8.6f k=%d  parry=%8.6f e=%d  witt = %8.6f diff=%8.6f\n",
		       i, mid, kn, K*par, en, witt, diff);

		mid = downshift(mid, K);
		par = beta_xform(par, beta);

		bn *= beta;
	}
// #define VERIFY_MIDPOINTS
#ifdef VERIFY_MIDPOINTS
	double K = atof(argv[1]);
	double beta = 2.0*K;

	double mid = K;
	double par = 1.0;
	par = beta_xform(par, beta);
	for (int i=0; i< 30; i++)
	{
		double mone = 2.0*mid - floor(2.0*mid);
		printf("%d 2xmid mod1=%g parry=%g diff = %g\n", i, mone, par, 2*mid-par);
		mid = downshift(mid, K);
		par = beta_xform(par, beta);
	}
#endif

#ifdef POLYNOMIALS
#define NPTS 1601
	for (int i=0; i<NPTS; i++)
	{
		double x = ((double) i + 0.5) / ((double) NPTS);
		double beta = 1.0+x;
		printf("%d	%g	%g	%g\n", i, beta, skippy(beta), skipry(beta));
	}
#endif
}
