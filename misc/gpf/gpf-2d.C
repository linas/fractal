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
	long double fact = 1.0;

	if (cabsl(x) < MAX_PREC) return x;

double scale = expl(-3.0L * cabsl(x));
// printf("duuude scale= %g\n", scale);
xn *= scale;
// fact *= expl(0.5 * cabsl(x));

	for (int n=1; ; n++)
	{
		sum += gpf(n) * (xn * fact);
		xn *= x;
		fact /= n;
		if (n*cabsl(xn*fact) < MAX_PREC*cabsl(sum)) break;
		if (max_iter < n) break;
	}
scale = expl(2.0L * cabsl(x));
scale /= cabsl(x);
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


	return cabs(g);
	// return creal(g);
	// return 0.5 + 0.5 * atan2(cimag(g), creal(g))/M_PI;
}

DECL_MAKE_HEIGHT(ploto);
