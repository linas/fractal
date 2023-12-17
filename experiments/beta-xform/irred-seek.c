/*
 * irred-seek.c
 *
 * Exploration of sequence recursion.
 *
 * December 2023
 */

#include "irred-gold.c"

double interp(double x, int ord)
{
	long scale = 1 << ord;
	double frac = x * scale;
	long ifrac = floor (frac);

	long summit = theta_sum(scale-1);
	long isum = theta_sum(ifrac);
	double sum = isum;
	// sum /= scale;
	sum /= (double) summit;

	return sum;
}

int main(int argc, char* argv[])
{
	int ord = atoi(argv[1]);
	long nmax = 1 << ord;
	malloc_gold(nmax);

	int npts = 3*5*7*11;
	double delta = 1.0 / ((double) npts);
	for (int i=0; i<npts; i++)
	{
		double x = (i+0.5) * delta;

		double y = interp(x, ord);
		printf("%d	%f	%f\n", i, x, y);
	}
}

/* --------------------------- END OF LIFE ------------------------- */
