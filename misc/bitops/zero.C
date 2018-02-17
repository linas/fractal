/*
 * zero.C
 * Visualization of complex zeros of golden polynomials
 *
 * Ferbruary 2018
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

static void golden_zero (float *array,
                             int array_size,
                             double x_center,
                             double x_width,
                             double row,
                             int itermax,
                             double K)
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
		COMPLEX zn = 1.0;

		// Use itermax as the encoding for the bit-string.
		COMPLEX sum = 0.0;
		int bitstr = 2*itermax+1;
   	while (bitstr)
		{
			if (bitstr & 0x1) sum += zn;
			bitstr >>= 1;
			zn *= z;
		}

		sum = zn - sum;

		// And now, switch to the asuymptotic series.
		if (K < 0.0) sum /= zn; 

 		array[j] = abs(sum);

		double r = x*x + y*y;
		if (0.99 < r and r < 1.01) array[j] = 0.5;
	}
}

DECL_MAKE_BIFUR(golden_zero)
