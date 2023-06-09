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
#include "mp-gkw.h"

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
	static int init = 0;
	if (0 == init)
	{
		/* Set the precision (number of binary bits) = prec*log(10)/log(2) */
		mpf_set_default_prec (3.3*prec);
		init = 1;
	}

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
#define MATRIX_ELTS
#ifdef MATRIX_ELTS
	int p = 100.0 * x + 0.5;
	int m = 100.0 * y + 0.5;
	m = 100 - m;
	double gkw = ache_mp_mp(m,p);

#if 0
	if (m%2 == 1) gkw = -gkw;
	if (p%2 == 1) gkw = -gkw;
#endif
	gkw *= exp(sqrt(m*p));
	// gkw = fabs(gkw);
// printf ("%d %d %g\n", m, p, gkw);

#if HYPERBOLA_FIT
	// double bola = p*m;
	// if ((bola > 2500) && (bola < 2700)) gkw = 1e30;
	// bola *= m;
	// if ((bola > 12500) && (bola < 13700)) gkw = 1e30;
	// double bola = p*m*m + p*p*m;
	// if ((bola > 12500) && (bola < 14300)) gkw = 1e30;
	// double bola = sqrt(m*m -6.0 *m*p + p*p);

	// double bola = -pow(fabs(m-p), 1.5) + pow(m+p, 1.5);
	double bola = -pow(fabs(m-p), 1.75) + pow(m+p, 1.75);
	if ((bola > 1550) && (bola < 1622)) gkw = 1e30;
	if ((bola > 15550) && (bola < 15692)) gkw = 1e30;

	// double gola = sqrt(m+ m*p + p);
	double gola = -pow(fabs(m-p), 1.7) + pow(m+p, 1.7);
	if ((gola > 1850) && (gola < 1932)) gkw = 1e30;
	if ((gola > 7850) && (gola < 7982)) gkw = 1e30;
#endif // HYPERBOLA_FIT

#if 0
   // if (p*5 == m) gkw = 1.0e30;
   // if (m*4 == p) gkw = 1.0e30;
	int mm = 1.8284 * m;
   if (mm == p) gkw = 1.0e30;
   if (mm == p+1) gkw = 1.0e30;
   if (mm+1 == p) gkw = 1.0e30;
#endif

#endif

// 	double gkw = ache_smooth_mp(x, -y);
// printf ("%f %f %g\n", x, -y, gkw);
	return gkw;
}

DECL_MAKE_HEIGHT(gkw_operator);

/* --------------------------- END OF LIFE ------------------------- */
