/*
 * psibig.c
 *
 * Compute midpoints using bignum.
 * Dec 2017
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


/**
 * Compute hessenberg-matrix orthonormal basis elements wavefunction
 * midpoints using GMP bignum math.
 *
 * nbits: precision to use
 * midp: pointer to double array in which to store the midpoints
 * digs: pointer to int array in which to store bit sequence.
 * maxn: how many to compute.... use MAXN from other code
 * kind: return either floats or binary digits.
 *
 * The midpoints are NOT sorted prior to return. The are listed
 * in the sequence in which they are found.
 *
 * The first midpoint is set to zero. The second midpoint is K.
 * The third is the first "true" midpoint, set to T(K).
 *
 * Each digit records whether 0.5<=midp (true) or not (false).
 * This is the Renyi-Parry sequence. The first bit is zero;
 * the second bit is one, the third bit is 0.5 <= T(K)
 */

double
midpoint_seq(double K, int nbits, double* midp, int* digs, int maxn)
{
	mpf_set_default_prec(nbits);

	mpf_t half;
	mpf_init(half);
	mpf_set_ui(half, 1);
	mpf_div_ui(half, half, 2);

	mpf_t Kay, ex, twoK;
	mpf_init(Kay);
	mpf_init(ex);
	mpf_init(twoK);

	// This deno is about 57 bits or so
	unsigned long int deno = 121*13*17*19*23*29UL;
	deno *= 128*81*125*49UL;

	// This deno is about 29 bits or so.
	// unsigned long int deno = 64*81*125*7*11*13;

	// unsigned long int deno = 32*27*25;
	unsigned long int num = K * ((double) deno);
	while (0 == num%2 && 0 == deno%2) { num /=2; deno /=2; }
	while (0 == num%3 && 0 == deno%3) { num /=3; deno /=3; }
	while (0 == num%5 && 0 == deno%5) { num /=5; deno /=5; }
	while (0 == num%7 && 0 == deno%7) { num /=7; deno /=7; }

	mpf_set_ui(Kay, num);
	mpf_div_ui(Kay, Kay, deno);
	mpf_mul_ui(twoK, Kay, 2);

	mpf_set(ex, Kay);
	printf("# K = %lu / %lu = %20.17g = %20.17g\n",
		num, deno, K, mpf_get_d(Kay));

	if (midp) midp[0] = 0.0;
	if (digs) digs[0] = 0;

	/* OK, now start iterating the map */
	for (int i=1; i < maxn; i++)
	{
		if (midp)
		{
			midp[i] = mpf_get_d(ex);
			// printf("# midpoint %d %g\n", i, midp[i]);
		}

		if (0 <= mpf_cmp(ex, half))
		{
			mpf_sub(ex, ex, half);
			if (digs) digs[i] = 1;
		}
		else if (digs) digs[i] = 0;
		mpf_mul(ex, ex, twoK);
	}

	// Return the actual value used, which might differ
	// from the provided value.
	return mpf_get_d(Kay);
}

double big_midpoints(double K, int nbits, double* midp, int maxn)
{
	return midpoint_seq(K, nbits, midp, 0x0, maxn);
}
