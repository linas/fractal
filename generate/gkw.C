/*
 * gkw.C
 *
 * FUNCTION:
 * Display GKW operator, one matrix element per pixel.
 *
 * HISTORY:
 * December 2003
 * Pixelize Jan 2010
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"
#include "binomial.h"
#include "harmonic.h"
#include "../misc/bignum/mp-zeta.h"
#include "../misc/bignum/mp-binomial.h"

// Return the matrix element for H_mp aka the matrix element of GKW.
//
long double
ache_mp_double(int m, int p)
{
	int k;

	long double acc = 0.0L;
	long double sign = 1.0L;
	for (k=0; k<=p; k++)
	{
		long double term = zetam1 (k+m+2);
		term *= binomial (m+k+1,m);
		term *= binomial (p,k);
		term *= sign;
		acc += term;
		sign = -sign;
	}
	return acc;
}

long double
ache_mp_mp(int m, int p)
{
	int k;

	int prec = 400;
	/* Set the precision (number of binary bits) = prec*log(10)/log(2) */
	mpf_set_default_prec (3.3*prec);

	mpf_t acc, one, term, fbin;
	mpf_init (term);
	mpf_init (acc);
	mpf_init (one);
	mpf_init (fbin);

	mpz_t bin;
	mpz_init (bin);

	// long double acc = 0.0L;
	mpf_set_ui (acc, 0);
	mpf_set_ui (one, 1);

	for (k=0; k<=p; k++)
	{
		// long double term = zetam1 (k+m+2);
		fp_zeta (term, k+m+2, prec);
		mpf_sub(term, term, one);

		// term *= binomial (m+k+1,m);
		i_binomial (bin, m+k+1, m);
		mpf_set_z (fbin, bin);
		mpf_mul (term, term, fbin);

		// term *= binomial (p,k);
		i_binomial (bin, p, k);
		mpf_set_z (fbin, bin);
		mpf_mul (term, term, fbin);

		if (k%2 == 0) mpf_add (acc, acc, term);
		else mpf_sub (acc, acc, term);
	}

	return mpf_get_d (acc);
}


static double gkw_operator (double x, double y, int itermax, double param)
{
	int p = 100.0 * x + 0.5;
	int m = 100.0 * y + 0.5;
	m = 100 - m;
	double gkw = ache_mp_mp(m,p);
	// gkw = fabs(gkw);
// printf ("%d %d %f\n", m, p, gkw);
	return gkw;
}

DECL_MAKE_HEIGHT(gkw_operator);

/* --------------------------- END OF LIFE ------------------------- */
