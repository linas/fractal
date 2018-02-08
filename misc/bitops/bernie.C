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

double island(double x, double K, double epsilon)
{
	K *= 2.0;
	if (0.5+epsilon <= x)
	{
		return K * (x - 0.5);
	}
	if (x < 0.5-epsilon)
	{
		return K*x;
	}
	return 0.25*K* (1.0 + (1.0 - 4.0*epsilon) * (0.5 - x) / epsilon);
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
	// return n+1; // mang-carry-1-more

	if (0 == n) return n;
	return n-1;  // mang-carry-1-less

	// return n%2;     // mang-carry-mod-2.png aka shift-XOR
	// return 1-n%2;   // mang-flip-mod.png

	// return n%3;

	// int max = 1;    // mang-carry-1
	// int max = 2;    // mang-carry-2
	int max = 3;    // mang-carry-3
	// int max = 6;
	if (max < n) return max;
	return n;
}

double mangle_carry(double x, double K)
{
	if (0.5 <= x)
	{
		x -= 0.5;
	}
	return 2.0* mangle_multiply(K, x, clamp);
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

	K = 0.5 + 0.5*K;
	double eps = 0.5* (1.0-K);
	eps = omega;
	for (int j=0; j<itermax; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		double x = t;

		/* OK, now start iterating the benoulli map */
		for (int iter=0; iter < 1250; iter++)
		{
			// x = bern(x, K);
			// x = nocarry(x, K);
			// x = tent(x, K);
			// x = notent(x, K);
			// x = feig(x, K);
			// x = nofeig(x, K);
			// x = mangle_carry(x, K);
			x = island(x, K, eps);

			double en = array_size * (x-floor(x));
			int n = en;
			if (0 > n) n = 0;
			if (n >= array_size) n = array_size-1;
			array[n] += 1.0;
			cnt ++;
		}
	}

	double wi = 0;
	double cut = 0.02 * cnt / ((double) array_size);
	for (int j=0; j<array_size; j++)
	{
		if (cut < array[j]) wi += 1.0;
	}

	// double norm = ((double) array_size) / ((double) cnt);
	double norm = wi / ((double) cnt);
	for (int j=0; j<array_size; j++)
		array[j] *= norm;
}

DECL_MAKE_BIFUR(bifurcation_diagram)

#if 0
int main ( int argc, char * argv[])
{
	double om = atof(argv[1]);
	double kb = atof(argv[2]);
}
#endif
