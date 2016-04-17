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

#define MAX_PREC 1.0e-8
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
double complex gpf_exponential(double complex x)
{
	double complex sum = 0;
	double complex xn = x;
	double fact = 1.0;

	if (cabs(x) < MAX_PREC) return x;

	for (int n=1; ; n++)
	{
		sum += gpf(n) * xn * fact;
		xn *= x;
		fact /= n;
		if (n*cabs(xn)*fact < MAX_PREC*cabs(sum)) break;
		if (max_iter < n) break;
	}

	return sum;
}

static double ploto(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
   double complex z = re_q + I * im_q;

	// double complex g = gpf_ordinary(z);
	// double complex g = gpf_exponential(z);
	double complex g = gpf_normed(z);

	// return cabs(g);
	// return creal(g);
	return 0.5 + 0.5 * atan2(cimag(g), creal(g))/M_PI;
}

DECL_MAKE_HEIGHT(ploto);
