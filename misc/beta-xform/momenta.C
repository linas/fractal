/*
 * momenta.C
 *
 * Feignebaum-style "invariant measure" type maps, but given a 
 * constant angular momentum.  Turns out "nothing" happens, its
 * boring, in the sense that its fully ergodic; the sum averages
 * out to zero for any non-zero momentum.  I guess this should
 * not be a surprise; its uniformly ergodic on the flat regions.
 *
 * Jan 2019
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "bitops.h"
#include "brat.h"

#define COMPLEX std::complex<double>
/*-------------------------------------------------------------------*/
/*
 */

// beta xform, except K = beta/2 and the range is different.
double downshift(double x, double K)
{
	K *= 2.0;
	if (0.5 <= x)
	{
		return K * (x - 0.5);
	}
	return K*x;
}

static void bifurcation_diagram (float *array,
                                 int array_size,
                                 double x_center,
                                 double x_width,
                                 double why,
                                 int itermax,
                                 double omega)
{
	/* clear out the row */
	for (int j=0; j<array_size; j++) array[j] = 0.0;

	COMPLEX carray[array_size];
	for (int j=0; j<array_size; j++) carray[j] = 0.0;

	int cnt=0;

	// 0 < why < 1 so why = beta - 1
	// double K = 0.5 + 0.5*why;
	double K = omega;
	COMPLEX step = cexp(why*2.0*M_PI*I);

	// printf("duuude K=%g  why=%g\n", K, why);
	for (int j=0; j<itermax; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		double x = t;  // 0 < x < 1

#ifdef RANDOMIZE_K
		t = rand();
		t /= RAND_MAX; // 0 < t < 1
		t -= 0.5;      // -0.5 < t < +0.5
		t /= 800.0; // 800 pixels tall
		t *= 0.5; // K runs from 0.5 to 1.0
		t *= x_width; // in case its zoomed.
		K = Korg + t;
#endif

		// COMPLEX z = 1.0;
		COMPLEX z = cexp(x*2.0*M_PI*I);
		/* OK, now start iterating the beta map */
		for (int iter=0; iter < 50; iter++) // number here causes banding...
		{
			x = downshift(x, K);
			z *= step;

			double en = array_size * (x-floor(x));
			int n = en;
			if (0 > n) n = 0;
			if (n >= array_size) n = array_size-1;
			carray[n] += z;
			cnt ++;
		}
	}

	for (int j=0; j<array_size; j++)
	{
		array[j] = std::abs(carray[j]);
		// array[j] = std::real(carray[j]);
		array[j] *= ((double) array_size) / ((double) cnt);
	}
}

DECL_MAKE_BIFUR(bifurcation_diagram)
