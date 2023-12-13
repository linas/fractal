/*
 * generate.c
 *
 * Exploration of generating functions for the beta-bitstring
 * and beta-sequence.
 *
 * December 2023
 */

#include "irred-gold.c"

// Ordinary generating function
double beta_ogf(double x)
{
	double sum=0.0;
	double xn = 1.0;
	for (int i=0; i<50; i++)
	{
		double beta = find_gold(i);
		if (0.5 < beta) sum += xn;
		printf("%d beta=%g sum=%g\n", i, beta, sum);
		xn *= x;
	}
	return sum;
}

int main(int argc, char* argv[])
{
	long nmax = 513;
	malloc_gold(nmax);

	printf("hello %20.16g\n", beta_ogf(0.5));
}
