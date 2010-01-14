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
#include "../misc/bignum/mp-gkw.h"

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
	int prec = 400;
	/* Set the precision (number of binary bits) = prec*log(10)/log(2) */
	mpf_set_default_prec (3.3*prec);

	mpf_t matelt;
	mpf_init (matelt);

	gkw(matelt, m, p, prec);

	return mpf_get_d (matelt);
}


// Return the continuous-valued version of the GKW operator.
// (the matrix matelts occur at integer values)
// This implementation uses GMP multi-precision
long double
ache_smooth_mp(double m, double p)
{
	int prec = 400;
	/* Set the precision (number of binary bits) = prec*log(10)/log(2) */
	mpf_set_default_prec (3.3*prec);

	mpf_t acc;
	mpf_init (acc);

	gkw_smooth(acc, m, p, prec);

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
