/*
 * irred-disk.c
 *
 * Exploration of Lambert-series-inspired generating functions
 * I.E. generation with the the beta polynomials.
 * The Lambert series placed a pole at each of the roots of unity.
 * This does likewise, but for the roots of beta.
 *
 * December 2023
 */

#include "brat.h"
#include "irred-gold.c"

#include <iostream>
#include <complex.h>
#define COMPLEX std::complex<double>
// #define COMPLEX double complex

// Complex polynomial
/* Return the beta value corresponding to the n'th golden polynomial.
 * It is be constructed from the bit string of (2n+1). Construction
 * is the mid-point construction: repeated iteration of the midpoint
 * 1/2 with this beta will (re-)generate the same bitstring, until
 * returning to the midpoint. The bits are just whether the orbit went
 * left or right of midpoint. The length of the orbit will be log_2(2n+1).
 */
COMPLEX cpx_golden_poly(long n, COMPLEX x)
{
	COMPLEX acc = 0.0;
	COMPLEX xn = 1.0;
	unsigned long bitstr = 2*n+1;
	while (bitstr)
	{
		if (bitstr%2 == 1) acc += xn;
		xn *= x;
		bitstr >>= 1;
	}
// printf("duuude n=%d x=%20.16g beta=\n", n, x, xn-acc);
	return xn - acc;
}

COMPLEX golden_recip(long n, COMPLEX x)
{
	return 1.0 / cpx_golden_poly(n, x);
}

#define MAXSUM 500

// complex-valued Ordinary generating function
COMPLEX COGF(COMPLEX (*fun)(long, COMPLEX), COMPLEX x)
{
	COMPLEX sum=0.0;
	COMPLEX xn = 1.0;
	for (int i=1; i<MAXSUM; i++)
	{
		sum += fun(i,x) * xn;
		xn *= x;
		if (abs(xn) < 1.0e-14) break;
	}
	return sum;
}

// coomplex-valued Exponential generating function
COMPLEX CEGF(COMPLEX (*fun)(long, COMPLEX), COMPLEX x)
{
	COMPLEX sum=0.0;
	COMPLEX xn = 1.0;
	for (int i=1; i<MAXSUM; i++)
	{
		sum += fun(i,x) * xn;
		xn *= x / ((double)(i+1));
		if (abs(xn) < 1.0e-14) break;
	}
	// std::cout << "zee=" << x << " sum=" << sum << std::endl;
	return sum;
}

static double beta_disk(double re_q, double im_q, int itermax, double param)
{
	static bool init = false;
	if (not init)
	{
		long nmax = 1UL << 20;
		malloc_gold(nmax);
		init = true;
	}

	COMPLEX zee = re_q + I * im_q;
	// COMPLEX og = COGF(cpx_golden_poly, zee);
	// COMPLEX og = COGF(golden_recip, zee);
	// COMPLEX og = CEGF(cpx_golden_poly, zee);
	COMPLEX og = CEGF(golden_recip, zee);

#if 1
	double faby = abs(og);
	double abz = abs(zee);

	// This norm is totally wrong.
	double norm = abz * exp(-0.6666*abz);
	// printf("u %g %g %g %g \n", re_q, im_q, faby, norm);
	faby *= norm;
	return faby;
#endif

#if 0
	double frea = real(og);
	double fima = imag(og);
	double phase = atan2 (fima, frea);
   phase += M_PI;
   phase /= 2.0*M_PI;
   return phase;
#endif
}

DECL_MAKE_HEIGHT(beta_disk);

/* --------------------------- END OF LIFE ------------------------- */
