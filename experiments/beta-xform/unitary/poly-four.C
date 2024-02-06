/*
 * poly-four.C
 *
 * Fourier at different betas.
 * Dec 2017
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "selfie.c"
#include "brat.h"

#define COMPLEX std::complex<double>

/*-------------------------------------------------------------------*/
void beta_fourier (float *array,
                   int array_size,
                   double x_center,
                   double x_width,
                   double y,
                   int itermax,
                   double param)
{
	/* clear out the row */
	for (int j=0; j<array_size; j++) array[j] = 0.0;

	COMPLEX phi = cexp(I*2.0*M_PI*y);
	for (int idx=0; idx<itermax; idx++)
	{
		if (false == valid_gold_index(idx)) continue;

		unsigned long tno = 2UL * idx + 1UL;
		int len = bitlen(tno);

		COMPLEX phin = phi;
		COMPLEX sum = -1.0;
		for (int i=0; i< len; i++)
		{
			if (tno & 1UL) sum += phin;
			tno >>= 1;
			phin *= phi;
		}

		double re = real(sum);
		double im = imag(sum);

		double beta = golden_beta(idx);
		int j = array_size* (beta - 1.0);
		array[j] = sqrt(re*re + im*im);
		// array[j] = 0.5 + atan2(im, re) / (2.0*M_PI);
	}
}

DECL_MAKE_BIFUR(beta_fourier)
