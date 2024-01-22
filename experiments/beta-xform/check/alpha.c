/*
 * alpha.c
 *
 * Examine the integer sequence alpha.  Lets see if its recognizable!?
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
		printf("idx=%ld  ", idx);
		for (int n=0; n<10; n++)
		{
			int an = alpha_n(bitstr, n);
			printf("%d  ", an);
		}
		printf("\n");
	}
}

// ==============================================================
