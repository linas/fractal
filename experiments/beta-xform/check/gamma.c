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

double t_ke(double beta, int k)
{
	double tk = 1.0;
	for (int i=0; i<k; i++)
	{
		tk *= beta;
// EPS is needed to prevent idx=3 from becoming cyclic
#define EPS 2e-15
		if (1.0 <= tk+EPS) { tk -= 1.0; if (tk < 0.0) tk = 0.0; }
	}
	return tk;
}

double t_kb(double beta, unsigned long bitstr, int k)
{
	double tk = 0.0;
	double ben = 1.0 / beta;
	int bl = bitlen(bitstr);
	if (bl < k) k = bl;
	for (int i=0; i<k; i++)
	{
		ben *= beta;
		if (0 == bit_k(bitstr, i)) continue;
		tk += 1.0/ben;
	}
	tk *= ben;
	ben *= beta;
	tk = ben - tk;
	if (tk < 0.0) tk = 0.0;  // needed for idx=3
	return tk;
}

double gamma_n(double beta, unsigned long bitstr, int n)
{
	double sum = 0.0;
	double ben = 1.0;
	for (int k=1; k<n; k++)
	{
		ben *= beta;
		if (0 == bit_k(bitstr, k)) continue;
		sum += t_kb(beta, bitstr, k) / ben;
	}
	return delta_n(bitstr, n) - sum;
}

double gamma_na(double beta, unsigned long bitstr, int n)
{
	double sum = 0.0;
	for (int k=1; k<n; k++)
	{
		if (0 == bit_k(bitstr, k)) continue;
		double bem = 1.0;
		for (int m=0; m<k; m++)
		{
			bem *= beta;
			if (0 == bit_k(bitstr, m)) continue;
			sum += 1.0 / bem;
		}
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
		if (false == valid_gold_index(idx)) continue;

		double beta = golden_beta(idx);

		unsigned long bitstr = 2*idx+1;
		int ord = bitlen(bitstr);
		printf("idx=%2ld ord=%d bits=", idx, ord);
		for (int n=0; n<ord; n++)
		{
			int bk = bit_k(bitstr, n);
			printf("%d", bk);
		}

		for (int n=0; n<ord; n++) printf(" ");
		printf("zeta= ");
		for (int n=0; n<12; n++)
			printf("%d  ", zeta_n(bitstr, n));
		printf("\n");

		printf("    beta=%.10f", beta);
		printf("  tk=");
		for (int n=0; n<=ord; n++)
			printf("%.4f  ", t_kb(beta, bitstr, n));
		printf("\n");

      printf("                  ");
		printf("  gamma=");
		for (int n=0; n<=ord; n++)
			printf("%.4f  ", gamma_n(beta, bitstr, n));
		printf("\n");

      printf("                  ");
		printf("  bamma=");
		for (int n=0; n<=ord; n++)
			printf("%.4f  ", gamma_na(beta, bitstr, n));
		printf("\n");

		printf("---\n");
	}
}

// ==============================================================
