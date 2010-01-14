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
#include "../misc/bignum/mp-binomial.h"
#include "../misc/bignum/mp-misc.h"
#include "../misc/bignum/mp-zeta.h"

// Return the matrix element for H_mp aka the matrix element of GKW.
// This implementation uses double-precision floats
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

// Return the matrix element for H_mp aka the matrix element of GKW.
// This implementation uses GMP multi-precision
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


// Return the continuous-valued version of the GKW operator.
// (the matrix elts occur at integer values)
// This implementation uses GMP multi-precision
long double
ache_smooth_mp(double m, double p)
{
	int k;

	int prec = 400;
	/* Set the precision (number of binary bits) = prec*log(10)/log(2) */
	mpf_set_default_prec (3.3*prec);

	mpf_t acc, one, term, bin;
	mpf_init (term);
	mpf_init (acc);
	mpf_init (one);
	mpf_init (bin);

	cpx_t ess, zeta;
	cpx_init(ess);
	cpx_init(zeta);

	// long double acc = 0.0L;
	mpf_set_ui (acc, 0);
	mpf_set_ui (one, 1);

	int ip = (int) floor(p);
	for (k=0; k<= ip; k++)
	{
		// long double term = zetam1 (k+m+2);
		// fp_zeta (term, k+m+2, prec);
		double km2 = k + m + 2.0L;
		cpx_set_d (ess, km2, 0.0);
		cpx_borwein_zeta(zeta, ess, prec);
		mpf_sub(term, zeta[0].re, one);
printf ("duuude k+m+2=%f\n", km2);
fp_prt("zeta -1 = ", term)

		// term *= binomial (m+k+1,m);
		double km1 = m+k+1.0L;
		fp_binomial_d (bin, km1, k+1);
		mpf_mul (term, term, bin);

		// term *= binomial (p,k);
		fp_binomial_d (bin, p, k);
		mpf_mul (term, term, bin);

		if (k%2 == 0) mpf_add (acc, acc, term);
		else mpf_sub (acc, acc, term);
	}

	return mpf_get_d (acc);
}


static double gkw_operator (double x, double y, int itermax, double param)
{
#ifdef MATRIX_ELTS
	int p = 100.0 * x + 0.5;
	int m = 100.0 * y + 0.5;
	m = 100 - m;
	double gkw = ache_mp_mp(m,p);
	// gkw = fabs(gkw);
// printf ("%d %d %g\n", m, p, gkw);
#endif

	double gkw = ache_smooth_mp(x, -y);
printf ("%f %f %g\n", x, -y, gkw);
	return gkw;
}

DECL_MAKE_HEIGHT(gkw_operator);

/* --------------------------- END OF LIFE ------------------------- */
