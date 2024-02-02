/*
 * bigun.c
 *
 * Compute "unitary" zeros using bignum.
 * Feb 2024
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


/**
 * Compute beta sequence for a given beta value.
 *
 * digs: pointer to byte array in which to store bit sequence.
 * maxn: how many bits to compute.
 *
 * Each digit records whether 0.5<=midp (true) or not (false).
 * This is the Renyi-Parry sequence. The zeroth bit is always one.
 * the first bit is 0.5 <= (m_1 = T(beta/2))
 */

void
gen_bitseq(mpf_t beta, char* digs, int maxn)
{
	mpf_t half;
	mpf_init(half);
	mpf_set_ui(half, 1);
	mpf_div_ui(half, half, 2);

	mpf_t midp;
	mpf_init(midp);
	mpf_set(midp, beta);
	mpf_div_ui(midp, midp, 2);

	/* OK, now start iterating the map */
	for (int i=0; i < maxn; i++)
	{
		if (0 <= mpf_cmp(midp, half))
		{
			mpf_sub(midp, midp, half);
			digs[i] = 1;
		}
		else digs[i] = 0;
		mpf_mul(midp, midp, beta);
	}

	mpf_clear(half);
	mpf_clear(midp);
}

int main(int argc, char* argv[])
{
	int bprec = 500;
	mpf_set_default_prec(bprec);

	// Set beta to exactly 1.6
	mpf_t beta;
	mpf_init(beta);
	mpf_set_ui(beta, 16);
	mpf_div_ui(beta, beta, 10);

#define NBITS 400
	char bitseq[NBITS];
	gen_bitseq(beta, bitseq, NBITS);
	for (int i=0; i<70; i++)
		printf("%d", bitseq[i]);
	printf("\n");
}
