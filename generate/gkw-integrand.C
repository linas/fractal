
/*
 * gkw-integrand.c
 *
 * FUNCTION:
 * examine the integrand of the GKW integral to find the 
 * saddle point(s).
 *
 * Linas Vepstas January 2010
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"
#include "binomial.h"
#include "harmonic.h"
#include "mp-zeta.h"


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
#define MATRIX_ELTS
#ifdef MATRIX_ELTS
	int p = 300.0 * x + 0.5;
	int m = 300.0 * y + 0.5;
	m = 300 - m;
	double gkw = ache_mp_mp(m,p);

	gkw *= exp(sqrt(m*p));
	// gkw = fabs(gkw);
// printf ("%d %d %g\n", m, p, gkw);

// 	double gkw = ache_smooth_mp(x, -y);
// printf ("%f %f %g\n", x, -y, gkw);
	return gkw;
}

DECL_MAKE_HEIGHT(gkw_operator);

/* --------------------------- END OF LIFE ------------------------- */
