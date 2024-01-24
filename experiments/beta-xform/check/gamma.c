/*
 * gamma.c
 *
 * Copy of alpha.c, for the beta-dpendent parts.
 *
 * January 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "selfie.c"
#include "unutil.c"

int bit_k(unsigned long bitstr, int k)
{
	int nu = bitlen(bitstr);
	return (bitstr >> (nu-k-1)) & 1UL;
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

double gamma_n(double beta, unsigned long bitstr, int n)
{
	double sum = 0.0;
	double ben = beta;
	for (int k=1; k<n; k++)
	{
		if (0 == bit_k(bitstr, k)) continue;
		sum += t_k(beta, k) / beta;
		ben *= beta;
	}
	return sum + delta_n(bitstr, n);
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
		double beta = golden_beta(idx);

		unsigned long bitstr = 2*idx+1;
		int ord = bitlen(bitstr);
		printf("idx=%2ld ord=%d bits=", idx, ord);
		for (int n=0; n<ord; n++)
		{
			int bk = bit_k(bitstr, n);
			printf("%d", bk);
		}

		printf("  alpha=");
		for (int n=0; n<8; n++)
			printf("%f  ", gamma_n(beta, bitstr, n));
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
