/*
 * zetadisk.C
 * Utter fail.
 * This just explores **obvious** properties of polynomials,
 * and for some reason I failed to realize this for hours. Stupid me.
 * The code works, but the meaning of the code is valueless.
 *
 * Radial slice through the (holomorphic) q-polynomial
 * This is a copy of betadisk.C but one-dimensional.
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

// Build the q-function.
// This is minus E(beta; zeta) where
// E(beta; zeta) = -1 + zeta * sum_n=0 zeta^n d_n(1/2)
COMPLEX qfunc(lodouble_t beta, COMPLEX zeta, int label)
{
	// #define SEQLEN 820
	#define SEQLEN 40
	static bool is_init = false;
	static int bit[SEQLEN+1];
	if (not is_init)
	{
		// is_init = true;

#define REGULAR_FLOAT_POINT
#ifdef REGULAR_FLOAT_POINT
		// This computes the bit sequence d_n(1/2)
		// AKA the bit-sequence Theta(T^n(beta/2) - 1/2)
		// This ... works ... but bignum does a sanity check.
		lodouble_t K = 0.5*beta;
		lodouble_t mid = K;
		for (int i=0; i<SEQLEN; i++)
		{
			bit[i] = 0;
			if (0.5 < mid) bit[i] = 1;
			mid = downshift(mid, K);
		}
#endif

// #define POLYNOMIAL_BITSTRINGS
#ifdef POLYNOMIAL_BITSTRINGS
		// Some hand-built infinite bit-strings, as an example:
		// bit[i] = 1 - i%2; // golden ratio n=1
		// bit[i] = (0 == i%3);  // n=2 bitstring 1001001001..
		// bit[i] = (0 == i%3) or (0 == (i+2)%3);  // n=3 bitstring 1101101101..
		// bit[i] = (0 == i%4);  // n=4 bitstring 100010001..
		// bit[i] = (0 == i%4) or (0 == (i+3)%4);  // n=6 bitstring 1100110011.
		// bit[i] = (0 == i%4) or (0 == (i+3)%4) or (0 == (i+2)%4);  // n=7 bitstring 1110111011.
		printf("Polynomial label n=%d\n", label);
		int pattern[24];
		int patlen = 0;
		for (int j=0; j<24; j++)
		{
			pattern[j] = label%2;
			label >>= 1;
			if (pattern[j]) patlen = j+1;
		}
		printf("bitstring length=%d\n", patlen+1);

#ifdef FINITE_STRINGS
		for (int i=0; i<SEQLEN; i++) bit[i] = 0;
		for (int i=0; i<patlen; i++)
		{
			bit[i] = pattern[patlen-i-1];
			printf("bit[%d] is %d\n", i, bit[i]);
		}
		bit[patlen] = 1;
		printf("bit[%d] is %d\n", patlen, bit[patlen]);
#endif

#define INFINITE_STRINGS
#ifdef INFINITE_STRINGS
		for (int i=0; i<SEQLEN; i++)
		{
			bit[i] = pattern[patlen - (i+1)%(patlen+1)];
			if (i < 40) printf("bit[%d] is %d\n", i, bit[i]);
		}
#endif // INFINITE_STRINGS

#endif // POLYNOMIAL_BITSTRINGS


// #define BIGNUM_MIDPOINTS
#ifdef BIGNUM_MIDPOINTS
		double K = 0.5*beta;
		midpoint_seq(K, SEQLEN+50, 0x0, bit, SEQLEN+1);

		// Above generates a sequence starting with zero;
		// below uses sequence starting with one...
		for (int i=0; i<SEQLEN; i++)
		{
			bit[i] = bit[i+1];
		}
#endif
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

// ================================================================

int main(int argc, char* argv[])
{
// #define RADIAL_SLICE
#ifdef RADIAL_SLICE
	// The radial slize looks at E(beta;z) along a radial line of zeta.
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s K theta\n", argv[0]);
		exit(1);
	}

	double K = atof(argv[1]);
	double beta = 2.0*K;

	double theta = atof(argv[2]);
	COMPLEX phase = cexp(I*theta*M_PI);

	printf("#\n# K=%g beta=%g theta=%g\n#\n", K, beta, theta);
#endif

	int label = 0;

#define NPTS 800
	for (int j=0; j<NPTS; j++)
	{
		double x = ((double) j + 0.5) / ((double) NPTS);

#ifdef RADIAL_SLICE
		double r = 2.0*x + 1.0;

		COMPLEX zeta = r* phase;

		COMPLEX sum = qfunc(beta, zeta, label);
		double mag = abs(sum);
		// double mag = sqrt(re*re + im*im);

		printf("%d	%g	%g\n", j, r, mag);
#endif

// #define BETA_SLICE
#ifdef BETA_SLICE
		// The beta slice fixes r=3 and looks at the magnitude
		// of E(beta;z) as beta is varied. The magnitude is huge.
		// it is always an integer power of r. (!)
		double beta = x + 1.0;
		double r = 3.0;
		double logr = log(r);

		printf("%d	%g", j, beta);

		for (double the=0.0; the< 1.0; the+=0.1)
		{
			COMPLEX phase = cexp(I*the*M_PI);
			COMPLEX zeta = r* phase;
			COMPLEX sum = qfunc(beta, zeta, label);
			double mag = abs(sum);
			double ponent = log(mag) / logr;
			printf("	%g", ponent);
		}
		printf("\n");
#endif

#define ZERO_COUNT
#ifdef ZERO_COUNT
		// The zero count fixes r=3 and looks for phase continuities
		// to count the number of zeros.
		double beta = x + 1.0;
		double r = 3.0;

		int num_zeros = 0;
		int num_up = 0;
		int num_down = 0;
		COMPLEX zeta = r;
		COMPLEX sum = qfunc(beta, zeta, label);
		double ere = real(sum);
		double eim = imag(sum);
		// ph ranges from -pi to +pi.
		double phprev = atan2(eim, ere);

#define NSTEPS 5000
		for (int a=0; a<NSTEPS; a++)
		{
			double the = ((double) a + 1) / ((double) NSTEPS);
			COMPLEX phase = cexp(I*the*2.0*M_PI);
			COMPLEX zeta = r* phase;
			COMPLEX sum = qfunc(beta, zeta, label);
			double ere = real(sum);
			double eim = imag(sum);

			// ph ranges from -pi to +pi.
			double ph = atan2(eim, ere);

			// This double-counts, but that's OK.
			if (ph*phprev < 0.0)
			{
				num_zeros ++;
				if (0.0 < ph) num_up ++;
				else num_down ++;
// printf("duude a=%d the=%g up=%d dn=%d phprev=%g  ph=%g\n", a, the, num_up, num_down, ph, phprev);
			}
			phprev = ph;
		}
		printf("%d	%g	%d	%d	%d\n", j, beta, num_zeros, num_up, num_down);
		if (0 != num_zeros % 2)
		{
			fprintf(stderr, "FAIL!!!!\n");
			exit(1);
		}
#endif
	}
}
