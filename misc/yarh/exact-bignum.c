/*
 * exact-bignum.c
 *
 * Implementation of exact-series. but in GMP
 *
 * Linas Vepstas December 2010
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <mp-binomial.h>
#include <mp-complex.h>

/* Integral of (x^(s-1))/(x+alpha) dx
 */
void eff (cpx_t sum, cpx_t s, mpf_t alpha, mpf_t x, int nprec)
{
	int k;
	mpf_t zero, epsi;
	mpf_init (zero);
	mpf_init (epsi);

	mpf_set_ui(epsi, 1);
	mpf_div_2exp(epsi, epsi, (int)(3.321*nprec));

	cpx_t term, ess;
	cpx_init (term);
	cpx_init (ess);

	mpf_t opxa, opxak;
	mpf_init (opxa);
	mpf_init (opxak);

	cpx_set(ess, s);
	cpx_set_ui(sum, 0, 0);

	mpf_div (opxa, x, alpha);
	mpf_add_ui (opxa, opxa, 1);
	mpf_set (opxak, opxa);

	for (k=1; k<155456123; k++)
	{
		cpx_binomial (term, ess, k);
		cpx_times_mpf (term, term, opxak);
		cpx_div_ui (term, term, k);
		if (k%2)
		{
			cpx_sub(sum, sum, term);
		}
		else
		{
			cpx_add (sum, sum, term);
		}

		/* Check to see if the last term is small enough */
		cpx_abs(zero, term);
		if (0 > mpf_cmp(zero, epsi)) break;

		mpf_mul (opxak, opxak, opxa);
	}

	cpx_clear(term);
	cpx_clear (ess);
	mpf_clear (opxa);
	mpf_clear (opxak);
	mpf_clear (epsi);
	mpf_clear (zero);
}

