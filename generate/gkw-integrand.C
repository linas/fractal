
/*
 * gkw-integrand.c
 *
 * FUNCTION:
 * examine the integrand of the GKW integral to find the 
 * saddle point(s).
 *
 * Linas Vepstas January 2010
 */


#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"
#include "binomial.h"
#include "harmonic.h"
#include "zmp.h"


// Return the integrand 
long double
grand(int m, int n, cpx_t s)
{
	mpf_t acc;
	mpf_init (acc);

	return mpf_get_d (acc);
}


static double gkw_integrand (double x, double y, int itermax, double param)
{
	static int init = 0;
	if (0 == init)
	{
		int prec = 400;
		/* Set the precision (number of binary bits) = prec*log(10)/log(2) */
		mpf_set_default_prec (3.3*prec);
		init = 1;
	}

	cpx_t s;
	cpx_init(s);
	cpx_set_d(s, x, y);

	int n = 30;
	double gkw = grand(n,n,s);

	return gkw;
}

DECL_MAKE_HEIGHT(gkw_integrand);

/* --------------------------- END OF LIFE ------------------------- */
