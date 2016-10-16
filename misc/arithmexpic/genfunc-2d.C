/*
 * Generating functions for miscellaneous arithmetic series
 * 2D phase plot.
 *
 * April 2016
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

#include <mp-zeta.h>

#include <gpf.h>
#include <modular.h>
#include <moebius.h>
#include <totient.h>

#include "brat.h"
#include "genfunc.h"



//  C and C++ is fugnuts insane in complex support.
#define complex _Complex

#define MAX_PREC 1.0e-18
int max_iter = 100000000;

/*
 * Ordinary generating function for arithmetic series
 */
double complex ordinary_genfunc(double complex x, int (*func)(int))
{
	double complex sum = 0;
	double complex xn = x;

	if (cabs(x) < 1.0e-16) return x;
	if (0.9999999 <= cabs(x)) return 0.0;

	for (int n=1; ; n++)
	{
		sum += func(n) * xn;
		xn *= x;
		if (n*cabs(xn) < MAX_PREC*cabs(sum)) break;
		if (max_iter < n) break;
	}

	return sum;
}

/*
 * Exponential generating function for arithmetic series
 */
double complex exponential_genfunc(long double complex x, int (*func)(int))
{
	long double complex sum = 0;
	long double complex xn = x;

	if (cabsl(x) < MAX_PREC) return x;

	for (int n=1; ; n++)
	{
		sum += func(n) * xn;
		xn *= x / ((long double) n+1);
		if (n*cabsl(xn) < MAX_PREC*cabsl(sum)) break;
		if (max_iter < n) break;
	}

	long double scale = expl(-cabsl(x));
	sum *= scale;
// printf("duuude %g sum=%g\n", cabs(x), cabs(sum));

	return sum;
}

double complex lambert_genfunc(double complex x, int (*func)(int))
{
	double complex sum = 0;
	double complex xn = x;

	if (cabs(x) < 1.0e-16) return x;
	if (0.9999999 <= cabs(x)) return 0.0;

	for (int n=1; ; n++)
	{
		sum += func(n) * xn / (1.0 - xn);
		xn *= x;
		if (n*cabs(xn) < MAX_PREC*cabs(sum)) break;
		if (max_iter < n) break;
	}

	return sum;
}

/* static */ double ploto(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
   double complex z = re_q + I * im_q;

      // double t = moebius_mu (i+1);
      // double t = mertens_m (i+1);
      // double t = liouville_omega (i+1);
      // double t = liouville_lambda (i+1);
      // double t = mangoldt_lambda (i+1);
      // double t = thue_morse (i+1);
      // int tm = thue_morse (i+1);

		// euler q-series aka dedekind eta,


	double complex g = ordinary_genfunc(z, totient_phi);
	// double complex g = gpf_exponential(z);
	// g *= cexp(-z);
	// g *= exp(-cabs(z)) / cabs(z);
	// double complex g = gpf_normed(z);
	// double complex g = gpf_lambert(z);

	// double rv = cabs(g);
	// double r = sqrt(re_q*re_q + im_q*im_q);
	// rv /= r;
	// rv /= r*r/log(r);
	// return rv;

	// return cabs(g);
	// return creal(g);

	return 0.5 + 0.5 * atan2(cimag(g), creal(g))/M_PI;
}

/* static */ double plot_big(double re_q, double im_q, int itermax, double param)
{
	cpx_t sum, z;
	cpx_init(sum);
	cpx_init(z);

	cpx_set_d(z, re_q, im_q);

// #define PHASE 1
#if PHASE
	// cpx_gpf_ordinary_recip(sum, z, 15);
	cpx_gpf_exponential(sum, z, 20);
	// cpx_gpf_poch_rising(sum, z, 45);
	// cpx_gpf_poch_falling(sum, z, 15);

	double rv = 0.5 + 0.5 * atan2(cpx_get_im(sum), cpx_get_re(sum))/M_PI;
	return rv;
#endif

}

DECL_MAKE_HEIGHT(ploto);
// DECL_MAKE_HEIGHT(plot_big);
