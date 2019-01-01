/*
 * qdisk.C
 * Visualization of the (holomorphic) q-polynomial
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
#define COMPLEX std::complex<double>

#include "brat.h"

// My downshift; beta=2K so T(x) = case bx or b(x-0.5)
double downshift(double x, double K)
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
double T_n(int n, double beta)
{
	double tn = 0.5*beta;

	// compute T^N(b/2)
	for (int i=0; i<n; i++)
	{
		tn = downshift(tn, 0.5*beta);
	}
	return tn;
}

// Build the q- function.
COMPLEX qfunc(double beta, COMPLEX zeta)
{
	#define SEQLEN 80
	static bool is_init = false;
	static int bit[SEQLEN];
	if (not is_init)
	{
		is_init = true;
		double K = 0.5*beta;
		double mid = K;
		for (int i=0; i<SEQLEN; i++)
		{
			bit[i] = 0;
			if (0.5 < mid) bit[i] = 1;
			mid = downshift(mid, K);
		}
	}

	COMPLEX zetan = 1.0;

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

		double r = x*x + y*y;
		if (r <= 1.0)
		{
			COMPLEX zeta = x + I*y;
			COMPLEX sum = qfunc(beta, zeta);
			array[j] = 0.5 + 0.5 * atan2(imag(sum), real(sum))/M_PI;

#if 0
			double mag = abs(sum);
			if (mag < 0.06) {
				COMPLEX z = beta*zeta;
				printf("maybe zero near z=%g +i %g = %g exp(i pi %g)\n",
					real(z), imag(z), abs(z), atan2(y,x)/M_PI);
array[j] = 0.5;
for (int k=0; k<20; k++) array[j-k]=0.25 * (k%4);
			}
#endif
		}
	}
}

DECL_MAKE_BIFUR(qpoly)
