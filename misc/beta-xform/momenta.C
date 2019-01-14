/*
 * momenta.C
 *
 * Feignebaum-style "invariant measure" type maps, but given a 
 * constant angular momentum.
 *
 * Jan 2019
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "bitops.h"
#include "brat.h"

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
                                 double beta,
                                 int itermax,
                                 double omega)
{
	/* clear out the row */
	for (int j=0; j<array_size; j++) array[j] = 0.0;

	int cnt=0;

	double Korg = 0.5 + 0.5*K;
	for (int j=0; j<itermax; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		double x = t;  // 0 < x < 1

		t = rand();
		t /= RAND_MAX; // 0 < t < 1
		t -= 0.5;      // -0.5 < t < +0.5
		t /= 800.0; // 800 pixels tall
		t *= 0.5; // K runs from 0.5 to 1.0
		t *= x_width; // in case its zoomed.
		K = Korg + t;

		/* OK, now start iterating the benoulli map */
		for (int iter=0; iter < 1250; iter++)
		{
			x = downshift(x, K);

			double en = array_size * (x-floor(x));
			int n = en;
			if (0 > n) n = 0;
			if (n >= array_size) n = array_size-1;
			array[n] += 1.0;
			cnt ++;
		}
	}
}

DECL_MAKE_BIFUR(bifurcation_diagram)
