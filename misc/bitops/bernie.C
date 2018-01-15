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

double nocarry(double x, double K)
{
	// Notes:
	// 1) The mult_xor function is technically incorrect;
	//    its too large by a factor of two.
	//    That's why we can skip multiplying K by two.
	// 2) The mult_xor function just shuffles, so we post-multiply by
	//    two, to get a return value between zero and one.
	if (0.5 <= x)
	{
		return 2.0 * mult_xor (K, (x - 0.5));
	}
	return 2.0*mult_xor(K, x);
}

int clamp(int n)
{
	int max = 1;
	if (max < n) return max;
	return n;
}

double mangle_carry(double x, double K)
{
	//    two, to get a return value between zero and one.
	if (0.5 <= x)
	{
		x -= 0.5;
	}
	return mangle_multiply(K, x, clamp);
}

double tent(double x, double K)
{
	K *= 2.0;
	if (0.5 <= x)
	{
		return K * (1.0 - x);
	}
	return K*x;
}

double notent(double x, double K)
{
	if (0.5 <= x)
	{
		return 2.0*mult_xor (K, (1.0 - x));
	}
	return 2.0*mult_xor(K, x);
}

double feig(double x, double K)
{
	K *= 4.0;
	return K * x * (1.0 - x);
}

double nofeig(double x, double K)
{
	// Below is intersting/crazy
	// return mult_xor(K, 4.0 * x * (1.0 - x));

	// Below is same as above. Multiplying by 2 commutes with mult_xor
	// return 2.0 * mult_xor(K, 2.0 * x * (1.0 - x));
	return 4.0 * mult_xor(K, x * (1.0 - x));

	// Below is totally uniform
	// return 2.0* mult_xor(K, 2.0 * mult_xor(x, (1.0 - x)));
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
		for (int iter=0; iter < 50; iter++)
		{
			// x = bern(x, K);
			// x = nocarry(x, K);
			// x = tent(x, K);
			// x = notent(x, K);
			// x = feig(x, K);
			x = nofeig(x, K);

			double en = array_size * (x-floor(x));
			int n = en;
			if (0 > n) n = 0;
			if (n >= array_size) n = array_size-1;
			array[n] += 1.0;
			cnt ++;
		}
	}

	for (int j=0; j<array_size; j++)
		array[j] *= ((double) array_size) / ((double) cnt);
}

DECL_MAKE_BIFUR(bifurcation_diagram)

#if 0
int main ( int argc, char * argv[])
{
	double om = atof(argv[1]);
	double kb = atof(argv[2]);
}
#endif
