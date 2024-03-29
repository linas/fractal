/*
 * alpha.c
 *
 * Examine the integer sequence alpha.  It is a generalized Fibonacci.
 * The linear sequence zeta is also a generalzed Fibonacci.
 * See also gamma.c for the beta-dpendent parts.
 *
 * January 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


/* Return length of bitstring. Same as ceil(log2(bitstr)). */
int bitlen(unsigned long bitstr)
{
	int len=0;
	while (bitstr) { len++; bitstr >>= 1; }
	return len;
}

int bit_k(unsigned long bitstr, int k)
{
	int nu = bitlen(bitstr);
	return (bitstr >> (nu-k-1)) & 1UL;
}

// alpha_n is an integer when nu_0=1 and turns out to be
// a generalized fibonacci.
int alpha_n(unsigned long bitstr, int n)
{
	if (n<2) return 1;

	int sum = 1;
	for (int k=1; k<n; k++)
	{
		int bk = bit_k(bitstr, k);
		if (0 == bk) continue;
		for (int m=0; m<n-k; m++)
			sum += alpha_n(bitstr, m);
	}
	return sum;
}

// Grand total number of non-zero bits.
int delta_n(unsigned long bitstr, int n)
{
	if (n<2) return 1;

	int sum = 0;
	for (int k=1; k<n; k++)
		sum += bit_k(bitstr, k);
	return sum;
}

// Generalized Fibonacci also, Duhh.
int zeta_n(unsigned long bitstr, int n)
{
	if (n<2) return 1;

	int sum = 1 + delta_n(bitstr, n);
	for (int k=1; k<n; k++)
	{
		int bk = bit_k(bitstr, k);
		if (0 == bk) continue;
		for (int m=1; m<n-k; m++)
			sum += zeta_n(bitstr, m);
	}
	return sum;
}

// ==============================================================

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s order\n", argv[0]);
		exit (1);
	}
	int nord = atoi(argv[1]);
	long nmax = 1UL << nord;

	for (long idx=1; idx<nmax; idx++)
	{
		unsigned long bitstr = 2*idx+1;
		int ord = bitlen(bitstr);
		printf("idx=%2ld ord=%d bits=", idx, ord);
		for (int n=0; n<ord; n++)
		{
			int bk = bit_k(bitstr, n);
			printf("%d", bk);
		}

		printf("  alpha=");
		for (int n=0; n<14; n++)
			printf("%d  ", alpha_n(bitstr, n));
		printf("\n");

		printf("                  ");
		for (int n=0; n<ord; n++) printf(" ");
		printf("  zeta= ");
		for (int n=0; n<14; n++)
			printf("%d  ", zeta_n(bitstr, n));

		printf("\n");
	}
}

// ==============================================================
