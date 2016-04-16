/*
 * Generating functions for greatest prime factors.
 * 2D phase plot.
 *
 * April 2016
 */

#include <math.h>
#include <stdio.h>

#include <gpf.h>
#include "brat.h"

/*
 * Ordinary generating function for the greatest common factor.
 */
double complex gpf_ordinary(double complex x)
{
	double complex sum = 0;
	double complex xn = x;

	if (cabs(x) < 1.0e-16) return x;

	for (int n=1; ; n++)
	{
		sum += gpf(n) * xn;
		xn *= x;
		if (n*cabs(xn) < 1.0e-16*cabs(sum)) break;
	}

	return sum;
}

/*
 * Exponential generating function for the greatest common factor.
 */
double complex gpf_exponential(double x)
{
	double complex sum = 0;
	double complex xn = x;
	double complex fact = 1.0;

	if (cabs(x) < 1.0e-16) return x;

	for (int n=1; ; n++)
	{
		sum += gpf(n) * xn * fact;
		xn *= x;
		fact /= n;
		if (n*cabs(xn)*fact < 1.0e-16*cabs(sum)) break;
	}

	return sum;
}

static double ploto(double re_q, double im_q, int itermax, double param)
{
   double x = param;

   double complex z = re_q + I * im_q;

	double g = gpf_ordinary(z);

	return 0.5 + 0.5 * atan2(cimag(sm), creal(sm))/M_PI;
}

DECL_MAKE_HEIGHT(ploto);
