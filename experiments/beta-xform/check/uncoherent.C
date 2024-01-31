/*
 * uncoherent.C
 *
 * Identical to unstack.c except a complex z is introduced into the sums.
 * The theory behind this is incomplete or broken. Not working as expected.
 *
 * January 2024
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

// ==============================================================

// Arbitrary function
double nu(double x)
{
	if (x < 0.0)
	{
		if (-1e-15 < x)
			x = 0.0;
		else
			fprintf(stderr, "Error nu fail neg %g\n", x);
	}
	if (1.0 < x) fprintf(stderr, "Error nu fail pos %g\n", x);

	// return 1.0;
	// return x-0.5;
	return x - x - 0.5*0.826154195;  // appropriate for beta=1.6

	// Bernoulli poly B_2
	// The result is senstive to this being B_2.
	// Being able to integrate to exactly zero is important.
	// return x*x - x  + 1.0 / 6.0;
	// return x*x - x  + 0.16666;

	// Bernoulli poly B_3
	// return x*x*x - 1.5*x*x  + 0.5*x;

	// Bernoulli poly B_4
	// return x*x*x*x - 2.0*x*x*x  + x*x - 1.0/30.0;
}

#include "uncomplex.C"

COMPLEX zap_n(double beta, COMPLEX blam, double x, int n)
{
	COMPLEX zap = cz_n(beta, blam, x, n+1);
	zap -= cz_n(beta, blam, x, n);
	zap *= pow(blam, n+1);

	return zap;
}

COMPLEX zip_n(double beta, COMPLEX blam, double x, int n)
{
	COMPLEX zip = nuz_n(beta, blam, x, n+1);
	zip -= nuz_n(beta, blam, x, n);
	zip *= pow(blam, n+1);

	return zip;
}

COMPLEX zoop_n(double beta, COMPLEX blam, double x, int n)
{
	COMPLEX zip = nuz_n(beta, blam, x, n+1);
	zip -= (blam/beta) * nuz_n(beta, blam, x, n);
	zip *= pow(blam, n+1);

	return zip;
}

// ==============================================================

static void coherent_zero (float *array,
                             int array_size,
                             double x_center,
                             double x_width,
                             double row,
                             int itermax,
                             double beta)
{
	/* clear out the row */
	for (int j=0; j<array_size; j++) array[j] = 0.0;

	double y = row;
	for (int j=0; j<array_size; j++)
	{
		double x = ((double) j + 0.5) / ((double) array_size);
		x -= 0.5;
		x -= x_center;
		x *= x_width;

		COMPLEX z = x + I*y;
		COMPLEX blam = beta*z;

		int n = itermax;
#if 1
		double yyy = 0.33;
		// COMPLEX sum = zap_n(beta, blam, yyy, n);
		// COMPLEX sum = zip_n(beta, blam, yyy, n);
		COMPLEX sum = zoop_n(beta, blam, yyy, n);
#endif

#ifdef AVG
		COMPLEX sum = 0.0;
		for (int i=0; i<6; i++)
		{
			double yyy = 0.181*i + 0.043678;
			sum += zip_n(beta, blam, yyy, n);
		}
#endif

#if 0
		double mag = abs(sum);
		if (mag < 1.0) printf("yo at %.10f %.10f\n", x, y);
#endif

		// array[j] = abs(sum);
		array[j] = 0.5 + 0.5 * atan2(imag(sum), real(sum))/M_PI;

		double r = x*x + y*y;
		if (0.99 < r and r < 1.01) array[j] = 0.5;

		r *= beta*beta;
		if (0.99 < r and r < 1.01) array[j] = 0.75;
	}
}

DECL_MAKE_BIFUR(coherent_zero)
