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

#include <mp-complex.h>

/* Integral of (x^(s-1))/(x+alpha) dx
 */
void eff (cpx_t sum, cpx_t s, mpf_t alpha, mpf_t x)
{
	int k;

	cpx_t term, ess;
	cpx_init (term);
	cpx_init (ess);

	mpf_t opxa, opxak;
	mpf_init (opxa);
	mpf_init (opxak);

	cpx_set(ess, s);
	cpx_set_ui(sum, 0, 0);

	mpf_div (opxa, x, alpha);
	mpf_add_ui (opxa, opxa, 1)
	mpf_set(opxak, opxa);


	cpx_clear(term);
}

