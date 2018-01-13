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


#define NBITS 4000
void big_midpoints(double K, int nbits, double* midp, int maxn)
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

	// unsigned long int deno = 128*81*125*49*121*13*17*19*23*29;
	unsigned long int deno = 64*81*125*7*11*13;
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
	printf("#\n# K = %lu / %lu = %20.17g = %20.17g\n#\n",
		num, deno, K, mpf_get_d(Kay));

	midp[0] = 0.0;

	/* OK, now start iterating the map */
	for (int i=1; i < maxn; i++)
	{
		midp[i] = mpf_get_d(ex);
printf("# duuude its %d %g\n", i, midp[i]);

		if (0 <= mpf_cmp(ex, half))
		{
			mpf_sub(ex, ex, half);
		}
		mpf_mul(ex, ex, twoK);
	}
}

