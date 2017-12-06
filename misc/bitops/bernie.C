/*
 * bernie.C
 *
 * Altered simplified algorithmic Bernoulli map
 * Dec 2017
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "bitops.h"
#include "brat.h"

/*-------------------------------------------------------------------*/
/*
 */

#if 0
static double xprod(double x, double y, int itermax, double param)
{
	// return mult_xor(x, y);
	return x*y;
}

DECL_MAKE_HEIGHT (xprod);
#endif

double bern(double x, double K)
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
                                 double K,
                                 int itermax,
                                 double omega)
{
	/* clear out the row */
	for (int j=0; j<array_size; j++) array[j] = 0.0;

	int cnt=0;

	for (int j=0; j<itermax; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		double x = t;

		/* OK, now start iterating the benoulli map */
		for (int iter=0; iter < 50; iter++) {
			x = bern(x, K);

			double en = array_size * (x-floor(x));
			int n = en;
			if (0 > n) n = 0;
			if (n >= array_size) n = array_size-1;
			array[n] += 1.0;
			cnt ++;
		}
	}
	
	for (int j=0; j<array_size; j++)
		array[j] /= cnt;
}

DECL_MAKE_BIFUR(bifurcation_diagram)

#if 0
int main ( int argc, char * argv[])
{
	double om = atof(argv[1]);
	double kb = atof(argv[2]);
}
#endif
