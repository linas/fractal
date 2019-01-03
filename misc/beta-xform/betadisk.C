/*
 * qdisk.C
 * Visualization of the (holomorphic) q-polynomial
 * This is exactly equal to minus alldisk.C ! Horay! It's starting to make sense!
 *
 * December 2018
 *
 */
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#include <complex.h>
#if 1
	#define COMPLEX std::complex<double>
	typedef double lodouble_t;
#else
	#define COMPLEX std::complex<__float128>
	typedef __float128 lodouble_t;
#endif

#include "brat.h"

// My downshift; beta=2K so T(x) = case bx or b(x-0.5)
lodouble_t downshift(lodouble_t x, lodouble_t K)
{
	K *= 2.0;
	if (0.5 <= x)
	{
		return K * (x - 0.5);
	}
	return K*x;
}

// Return the iterated downshift t^n(beta/2)
// The Renyi-Parry bit d_n is given by (x<tn)
lodouble_t T_n(int n, lodouble_t beta)
{
	lodouble_t tn = 0.5*beta;

	// compute T^N(b/2)
	for (int i=0; i<n; i++)
	{
		tn = downshift(tn, 0.5*beta);
	}
	return tn;
}

// Build the q-function.
COMPLEX qfunc(lodouble_t beta, COMPLEX zeta)
{
	#define SEQLEN 120
	static bool is_init = false;
	static int bit[SEQLEN];
	if (not is_init)
	{
		is_init = true;
		lodouble_t K = 0.5*beta;
		lodouble_t mid = K;
		for (int i=0; i<SEQLEN; i++)
		{
			bit[i] = 0;
			if (0.5 < mid) bit[i] = 1;
			mid = downshift(mid, K);
// bit[i] = 1 - i%2; // golden ratio n=1
// bit[i] = (0 == i%3);  // n=2 bitstring 1001001001..
// bit[i] = (0 == i%3) or (0 == (i+2)%3);  // n=3 bitstring 1101101101..
// bit[i] = (0 == i%4);  // n=4 bitstring 100010001..
// bit[i] = (0 == i%4) or (0 == (i+3)%4);  // n=6 bitstring 1100110011.
// bit[i] = (0 == i%4) or (0 == (i+3)%4) or (0 == (i+2)%4);  // n=7 bitstring 1110111011.
// printf("duuude its %d %d\n", i, bit[i]);
		}
	}

	COMPLEX zetan = zeta;

	// accumulated sum
	COMPLEX que = 1.0;

	for (int i=0; i<SEQLEN; i++)
	{
		if (bit[i]) que -= zetan;

		// compute z^n;
		zetan *= zeta;
	}

	return que;
}

// Build the exponential q-function.
COMPLEX exp_qfunc(lodouble_t beta, COMPLEX zeta)
{
	#undef SEQLEN
	#define SEQLEN 920
	static bool is_init = false;
	static int bit[SEQLEN];
	if (not is_init)
	{
		is_init = true;
		lodouble_t K = 0.5*beta;
		lodouble_t mid = K;
		for (int i=0; i<SEQLEN; i++)
		{
			bit[i] = 0;
			if (0.5 < mid) bit[i] = 1;
			mid = downshift(mid, K);
		}
	}

	COMPLEX zetan = zeta;

	// accumulated sum
	COMPLEX que = 1.0;

	for (int i=0; i<SEQLEN; i++)
	{
		if (bit[i]) que -= zetan;

		// compute z^n;
		zetan *= zeta;
		zetan /= i+1;
	}

	double r = abs(zeta);
	que *= exp (-r);
	que /= log(r+1.0);
	que *= 2.18-beta;

	return que;
}

// ================================================================

static void qpoly (float *array,
                             int array_size,
                             double x_center,
                             double x_width,
                             double row,
                             int itermax,
                             double K)
{
	/* clear out the row */
	for (int j=0; j<array_size; j++) array[j] = 0.0;

	double beta = 2.0*K;

	double y = row;
	for (int j=0; j<array_size; j++)
	{
		double x = ((double) j + 0.5) / ((double) array_size);
		x -= 0.5;
		x -= x_center;
		x *= x_width;

		// double r = x*x + y*y;
		// if (r <= 2.0)
		{
			COMPLEX zeta(x,y);
#ifdef QFUNC
			COMPLEX sum = qfunc(beta, zeta);
			// Take minus the sum, to get what alldisk.C is showing.
			sum = -sum;
			double re = real(sum);
			double im = imag(sum);
			double pha = atan2(im, re)/M_PI;
			array[j] = 0.5 + 0.5 * pha;

			if (1.0 < r and r <= 1.02) array[j] = 1;
#endif

			// A goofy game.
			COMPLEX sum = exp_qfunc(beta, zeta);
			array[j] = abs(sum);

#if 0
			double mag = abs(sum);
			if (mag < 0.15) {
				COMPLEX z = beta*zeta;
				printf("maybe zero near z=%g +i %g = %g exp(i pi %g)\n",
					real(z), imag(z), abs(z), atan2(y,x)/M_PI);
array[j] = 0.5;
for (int k=0; k<30; k++) array[j-k]=0.25 * (k%4);
			}
#endif
		}
	}
}

DECL_MAKE_BIFUR(qpoly)
