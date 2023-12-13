/*
 * irred-gf.c
 *
 * Exploration of generating functions for the beta-bitstring
 * and beta-sequence.
 *
 * December 2023
 */

#include "irred-gold.c"

// Ordinary generating function
double OGF(double (*fun)(long), double x)
{
	double sum=0.0;
	double xn = 1.0;
	for (int i=1; i<50; i++)
	{
		sum += fun(i) * xn;
		xn *= x;
	}
	return sum;
}

// The mask-bits
double mask(long n)
{
	double beta = find_gold(n);
	if (0.5 < beta) return 1.0;
	return 0.0;
}

// The golden values themselves, or zero.
double gold(long n)
{
	double beta = find_gold(n);
	if (0.5 < beta) return beta;
	return 0.0;
}

// Allowed values
double allowed(long n)
{
	long idx = 1;
	int cnt = 0;
	while (cnt < n)
	{
		double beta = find_gold(idx);
		if (0.5 < beta) cnt++;
		idx++;
	}
	idx --;
	// printf("counted %ld is %ld\n", n, idx);
	return idx;
}

int main(int argc, char* argv[])
{
	long nmax = 513;
	malloc_gold(nmax);

	printf("Mask OGF at 1/2 = %20.16g\n", OGF(mask, 0.5));
	printf("Gold OGF at 1/2 = %20.16g\n", OGF(gold, 0.5));
	printf("Allowed OGF at 1/2 = %20.16g\n", OGF(allowed, 0.5));
}
