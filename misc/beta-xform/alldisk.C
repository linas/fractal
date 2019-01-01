/*
 * alldisk.C
 * Visualization of the almost-eigen constants on the unit disk.
 * Holomorphic version of almost.c
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

// Build the "constant" function.
// This appears to be independent of the vaue of x.
// Which is surprising.
COMPLEX dcnst(double x, double beta, COMPLEX zeta)
{
	if (0.5*beta < x) return 0.0;

	double yoblo = x/beta;
	double yobhi = yoblo + 0.5;
	COMPLEX zetan = 1.0;

	double tn = 0.5*beta;

	// accumulated sum
	COMPLEX dee = 0.0;

	int cnt=0;
	while (1.0e-16 < abs(zetan))
	{
		COMPLEX term = 0.0;

		// Each dn is either zero or one... so this is easy.
		if (yoblo < tn) term += zeta;
		if (yobhi < tn) term += zeta;
		if (x < tn) term -= 1.0;

		dee += zetan * term;

		// compute T^N(b/2)
		tn = downshift(tn, 0.5*beta);

		// compute z^n;
		zetan *= zeta;

		cnt++;
		if (3000< cnt) break;
	}

	return dee;
}

// ================================================================

static void almost_zero (float *array,
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
	double undep = 0.001*itermax;

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
			COMPLEX sum = dcnst(undep, beta, zeta);
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

DECL_MAKE_BIFUR(almost_zero)
