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

#include <mp-arith.h>
#include <mp-genfunc.h>
#include <mp-trig.h>

#include <isqrt.h>
#include <moebius.h>
#include <necklace.h>

#include "brat.h"



//  C and C++ is fugnuts insane in complex support.
#define complex _Complex

#define MAX_PREC 1.0e-18
int max_iter = 100000000;

/*
 * Ordinary generating function for arithmetic series
 */
double complex ordinary_genfunc(double complex x, long (*func)(long))
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
double complex exponential_genfunc(long double complex x, long (*func)(long))
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

static double mobius_exp_mag(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
// #define UNCIRCLE
#ifdef UNCIRCLE
	double tmp = re_q;
	re_q = 1.0 - im_q;
	im_q = tmp;

	// max_iter = itermax;
	max_iter = 100000;
	double theta = M_PI * im_q;
	double rr = itermax + param * re_q;
	rr = exp(rr * M_LN2);  // pow (2, itermax + param * re_q)
// printf("duuude re=%g im=%g r = %g\n", re_q, im_q, rr);
	im_q = rr*sin (theta);
	re_q = rr*cos (theta);
#endif

	long double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, moebius_mu);
	return 0.5 + 0.5 * atan2(cimag(g), creal(g))/M_PI;
	// return cabs(g);
}

static double mobius_big(double re_q, double im_q, int itermax, double param)
{
	cpx_t sum, z; cpx_init(sum); cpx_init(z);
	mpf_t val; mpf_init(val);

	cpx_set_d(z, re_q, im_q);

	cpx_exponential_genfunc(sum, z, 25, moebius_mu);
	cpx_abs(val, sum);
	double rv = mpf_get_d(val);

	cpx_clear(sum);
	cpx_clear(z);
	mpf_clear(val);
	return rv;
}

static double divisor_exp_mag(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
	long double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, divisor);
	return cabs(g);
}

static double divisor_big(double re_q, double im_q, int itermax, double param)
{
	cpx_t sum, z; cpx_init(sum); cpx_init(z);
	mpf_t val; mpf_init(val);

	cpx_set_d(z, re_q, im_q);

	cpx_exponential_genfunc(sum, z, 25, divisor);
	cpx_abs(val, sum);
	double rv = mpf_get_d(val);

	cpx_clear(sum);
	cpx_clear(z);
	mpf_clear(val);
	return rv;
}

// ========================================================

__attribute__((constructor)) void decl_things() {
	DECL_HEIGHT("mobius_exp_mag", mobius_exp_mag);
	DECL_HEIGHT("mobius_big", mobius_big);
	DECL_HEIGHT("divisor_exp_mag", divisor_exp_mag);
	DECL_HEIGHT("divisor_big", divisor_big);
}

// DECL_MAKE_HEIGHT(plot_big);
MAKE_HEIGHT;
