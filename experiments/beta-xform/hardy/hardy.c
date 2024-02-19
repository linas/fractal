/*
 * hardy.c
 * Hardy-space limit measure
 *
 * Feb 2024
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "selfie.c"

// This is a buggy mistake, but an interesting one.  It accidentally
// reversed the bit-sequence.  oops(z) = z^n+1 - sum_k z^k b_k
// for the usual bitsequence b_k.  This is evaluated on the circle
// |z| = 1/beta. It's minorly curious. But whatever.
complex oops(const char* bitseq, int order, double beta, double angle)
{
	complex zee = cexp(I * 2.0 * M_PI * angle) / beta;

	complex zn = 1.0;
	complex sum = 0.0;
	for (int i=0; i< order; i++)
	{
		if (bitseq[i]) sum -= zn;
		zn *= zee;
	}
	sum += zn;
	return sum;
}

// The "analytic ergodics" function, E(beta,z) = -1 + sum_n b_n zeta^n
// is ... surprisingly boring on ring zeta=1 mostly cause of
// Jentzch theorem, which is just splattering zeros on the rim.
// So I dunno; what was I thinking?
complex hardy(const char* bitseq, int order, double angle)
{
	complex zee = cexp(I * 2.0 * M_PI * angle);

	complex zn = zee;
	complex sum = -1.0;
	for (int i=0; i< order; i++)
	{
		if (bitseq[i]) sum += zn;
		zn *= zee;
	}
	return sum;
}

// -------------------------------------------------------

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s <index>\n", argv[0]);
		exit(1);
	}
	long aidx = atol(argv[1]);
	long idx = aidx;

	while (false == valid_gold_index(idx))
		idx++;

	if (idx != aidx)
		fprintf(stderr, "Warning: bump %ld to %ld\n", aidx, idx);

	unsigned long tno = 2UL * idx + 1UL;
	int order = bitlen(tno);

#define NBITS 63
	char bitseq[NBITS];
	for (int i=0; i<order; i++)
	{
		bitseq[i] = (tno & 1UL);
		tno >>= 1;
	}

	printf("#\n# Index=%ld Order=%d bitseq= ", idx, order);
	for (int i=0; i<order; i++)
		printf("%d", bitseq[i]);
	printf("\n");

	double beta = golden_beta(idx);
	printf("# Beta = %f\n", beta);

#define NPTS 1123
	for (int i=0; i< NPTS; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NPTS);
		// complex h = hardy(bitseq, order, x);
		complex h = oops(bitseq, order, beta, x);

		printf("%d	%f	%f	%f\n", i, x, creal(h), cimag(h));
	}
}
