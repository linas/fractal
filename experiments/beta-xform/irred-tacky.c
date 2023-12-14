/*
 * irred-tacky.c
 *
 * Exploration of takagi-like generating functions for the beta-bitstring
 * and beta-sequence.
 *
 * December 2023
 */

#include "irred-gold.c"

// real-valued Ordinary generating function
double OGF(double (*fun)(int, double), double w, double x)
{
	double sum=0.0;
	double wn = 1.0;
	for (int i=1; i<50; i++)
	{
		sum += fun(i, x) * wn;
		wn *= w;
	}
	return sum;
}

double tak(int n, double x)
{
	long tn = 1 << n;
	double xn = x * ((double) tn);
	double xm = xn - floor(xn);
	int bn = (xm > 0.5);
	return bn;
}

int main(int argc, char* argv[])
{
	long nmax = 51;
	malloc_gold(nmax);

	double w = 0.6;

	int npts = 3*5*7*11;
	double delta = 1.0 / ((double) npts);
	for (int i=0; i<npts; i++)
	{
		double x = (i+0.5) * delta;

		double y = OGF(tak, w, x);
		printf("%d	%f	%f\n", i, x, y);
	}
}

/* --------------------------- END OF LIFE ------------------------- */
