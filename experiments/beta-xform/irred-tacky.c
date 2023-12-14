/*
 * irred-tacky.c
 *
 * Exploration of takagi-like generating functions for the beta-bitstring
 * and beta-sequence.
 *
 * December 2023
 */

#include "irred-gold.c"

#define MAXORD 16

// real-valued Ordinary generating function
double tak(double w, double x)
{
	double sum=0.0;
	double wn = 1.0;
	long idx = 0;
	for (int i=0; i<MAXORD; i++)
	{
		int bn = (x > 0.5);
		if (bn)
		{
			idx |= 1;
		}
		int vn = is_valid_index(idx);

		if (bn && vn) sum += wn;
		wn *= w;

		idx <<= 1;
		x *= 2.0;
		x -= floor(x);
	}
	return sum;
}

int main(int argc, char* argv[])
{
	long nmax = 1 << MAXORD;
	malloc_gold(nmax);

	double w = 0.6;

	int npts = 3*5*7*11;
	double delta = 1.0 / ((double) npts);
	for (int i=0; i<npts; i++)
	{
		double x = (i+0.5) * delta;

		double y = tak(w, x);
		printf("%d	%f	%f\n", i, x, y);
	}
}

/* --------------------------- END OF LIFE ------------------------- */
