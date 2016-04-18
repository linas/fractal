/*
 * Generating functions for greatest prime factors.
 * 2D phase plot.
 *
 * April 2016
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>

#include <gpf.h>
#include "brat.h"
#include "gpf-gen-bignum.h"

//  C and C++ is fugnuts insane in complex support.
#define complex _Complex

#define MAX_PREC 1.0e-18
int max_iter = 100000000;

/*
 * Ordinary generating function for the greatest common factor.
 */
double complex gpf_ordinary(double complex x)
{
	double complex sum = 0;
	double complex xn = x;

	if (cabs(x) < 1.0e-16) return x;
	if (0.9999999 <= cabs(x)) return 0.0;

	for (int n=1; ; n++)
	{
		sum += gpf(n) * xn;
		xn *= x;
		if (n*cabs(xn) < MAX_PREC*cabs(sum)) break;
		if (max_iter < n) break;
	}

	return sum;
}

double complex gpf_normed(double complex x)
{
	double complex sum = 0;
	double complex xn = x;

	if (cabs(x) < 1.0e-16) return x;
	if (0.9999999 <= cabs(x)) return 0.0;

	for (int n=1; ; n++)
	{
		sum += gpf(n) * xn / ((double) n);
		xn *= x;
		if (cabs(xn) < MAX_PREC*cabs(sum)) break;
		if (max_iter < n) break;
	}

	return sum;
}

/*
 * Exponential generating function for the greatest common factor.
 */
double complex gpf_exponential(long double complex x)
{
	long double complex sum = 0;
	long double complex xn = x;

	if (cabsl(x) < MAX_PREC) return x;

	for (int n=1; ; n++)
	{
		sum += gpf(n) * xn;
		xn *= x / ((long double) n);
		if (n*cabsl(xn) < MAX_PREC*cabsl(sum)) break;
		if (max_iter < n) break;
	}

	long double scale = expl(-cabsl(x));
	sum *= scale;
// printf("duuude %g sum=%g\n", cabs(x), cabs(sum));

	return sum;
}

double complex gpf_lambert(double complex x)
{
	double complex sum = 0;
	double complex xn = x;

	if (cabs(x) < 1.0e-16) return x;
	if (0.9999999 <= cabs(x)) return 0.0;

	for (int n=1; ; n++)
	{
		sum += gpf(n) * xn / (1.0 - xn);
		xn *= x;
		if (n*cabs(xn) < MAX_PREC*cabs(sum)) break;
		if (max_iter < n) break;
	}

	return sum;
}

static double ploto(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
   double complex z = re_q + I * im_q;

	// double complex g = gpf_ordinary(z);
	double complex g = gpf_exponential(z);
	// g *= cexp(-z);
	// g *= exp(-cabs(z)) / cabs(z);
	// double complex g = gpf_normed(z);
	// double complex g = gpf_lambert(z);

	double rv = cabs(g);
	double r = sqrt(re_q*re_q + im_q*im_q);
	// rv /= r;
	rv /= r*r/log(r);

	return rv;

	// return cabs(g);
	// return creal(g);
	// return 0.5 + 0.5 * atan2(cimag(g), creal(g))/M_PI;
}

static double plot_big(double re_q, double im_q, int itermax, double param)
{
	cpx_t sum, z;
	cpx_init(sum);
	cpx_init(z);

	cpx_set_d(z, re_q, im_q);

	cpx_gpf_exponential(sum, z, 20);

	// extract
	mpf_t val;
	mpf_init(val);
	cpx_abs(val, sum);

	double rv = mpf_get_d(val);

	// Divide by z for plotting.
	double r = sqrt(re_q*re_q + im_q*im_q);
	double lr = log(r);
	// rv /= r*r / (lr*lr);
	rv /= r;

	return rv;
}

// DECL_MAKE_HEIGHT(ploto);
DECL_MAKE_HEIGHT(plot_big);
