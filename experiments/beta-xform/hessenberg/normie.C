/*
 * normie.C
 * Derived from betadisk.C/zetadisk.C
 *
 * October 2020
 */
#include <gmp.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#include "psibig.c"

#include <complex.h>
#if 1
	#define COMPLEX std::complex<double>
	typedef double lodouble_t;
#else
	#define COMPLEX std::complex<__float128>
	typedef __float128 lodouble_t;
#endif

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

// Build the normalization integral function.
//  sum_n=0 zeta^n T^n(beta/2)
COMPLEX normie(lodouble_t beta, COMPLEX zeta)
{
	int nterms = -log(1.0e-7) / log(beta);

	// This computes the sum T^n(beta/2)
	COMPLEX sum = 0;
	COMPLEX zetan = 1;

	lodouble_t K = 0.5*beta;
	lodouble_t mid = K;
	for (int i=0; i<nterms; i++)
	{
		sum += mid * zetan;
		mid = downshift(mid, K);

		// compute z^n;
		zetan *= zeta;
	}

	return sum;
}

// Build the normalization integral function.
//  sum_n=0 zeta^n T^n(beta/2)
// same as above, high-precision
COMPLEX big_normie(lodouble_t beta, COMPLEX zeta)
{
	int nterms = -log(1.0e-7) / log(beta);

	double midpoints[nterms];

	big_midpoints(0.5*beta, nterms+20, midpoints, nterms);

	// This computes the sum T^n(beta/2)
	COMPLEX sum = 0;
	COMPLEX zetan = 1;

	for (int i=0; i<nterms; i++)
	{
		sum += zetan * midpoints[i+1];

		// compute z^n;
		zetan *= zeta;
	}

	return sum;
}

// ================================================================

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s mod(z) arg(z)\n", argv[0]);
		exit(1);
	}

	double absz = atof(argv[1]);
	double theta = atof(argv[2]);
	COMPLEX phase = cexp(I*theta*M_PI);
	COMPLEX zee = absz * phase;

	printf("#\n# mod(z)=%g arg(z)=%g\n#\n", absz, theta);

#define NPTS 2503
	for (int j=0; j<NPTS; j++)
	{
		double x = ((double) j + 0.5) / ((double) NPTS);
		double beta = x + 1.0;

		COMPLEX zeta = zee / beta;

		// COMPLEX sum = normie(beta, zeta);
		COMPLEX sum = big_normie(beta, zeta);
		double mag = abs(sum);
		// double mag = sqrt(re*re + im*im);
		double arg = atan2(imag(sum), real(sum));

		printf("%d	%g	%g	%g\n", j, beta, mag, arg);
		fflush(stdout);
	}
}
