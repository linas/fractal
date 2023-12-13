/*
 * generate.c
 *
 * Exploration of generating functions for the beta-bitstring
 * and beta-sequence.
 *
 * December 2023
 */

#include "irred-gold.c"

// Ordinary generating function of the mask-bits
double mask_ogf(double x)
{
	double sum=0.0;
	double xn = 1.0;
	for (int i=0; i<50; i++)
	{
		double beta = find_gold(i);
		if (0.5 < beta) sum += xn;
		// printf("%d beta=%g sum=%g\n", i, beta, sum);
		xn *= x;
	}
	return sum;
}

// Ordinary generating function of the gold values
double gold_ogf(double x)
{
	double sum=0.0;
	double xn = 1.0;
	for (int i=0; i<50; i++)
	{
		double beta = find_gold(i);
		if (0.5 < beta) sum += beta*xn;
		printf("%d beta=%g sum=%g\n", i, beta, sum);
		xn *= x;
	}
	return sum;
}

int main(int argc, char* argv[])
{
	long nmax = 513;
	malloc_gold(nmax);

	printf("Mask OGF at 1/2 = %20.16g\n", mask_ogf(0.5));
	printf("Gold OGF at 1/2 = %20.16g\n", gold_ogf(0.5));
}
