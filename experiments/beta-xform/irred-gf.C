/*
 * irred-gf.c
 *
 * Exploration of generating functions for the beta-bitstring
 * and beta-sequence.
 *
 * December 2023
 */

#include "brat.h"
#include "irred-gold.c"

#include <iostream>
#include <complex.h>
#define COMPLEX std::complex<double>
// #define COMPLEX double complex

// real-valued Ordinary generating function
double OGF(double (*fun)(long), double x)
{
	double sum=0.0;
	double xn = 1.0;
	for (int i=1; i<50; i++)
	{
		sum += fun(i) * xn;
		xn *= x;
	}
	return sum;
}

// coomplex-valued Ordinary generating function
COMPLEX COGF(double (*fun)(long), COMPLEX x)
{
	COMPLEX sum=0.0;
	COMPLEX xn = 1.0;
	for (int i=1; i<500; i++)
	{
		sum += fun(i) * xn;
		xn *= x;
		if (abs(xn) < 1.0e-14) break;
	}
	return sum;
}

// coomplex-valued Exponential generating function
COMPLEX CEGF(double (*fun)(long), COMPLEX x)
{
	COMPLEX sum=0.0;
	COMPLEX xn = 1.0;
	for (int i=1; i<500; i++)
	{
		sum += fun(i) * xn;
		xn *= x / ((double)(i+1));
		if (abs(xn) < 1.0e-14) break;
	}
	// std::cout << "zee=" << x << " sum=" << sum << std::endl;
	return sum;
}

// The mask-bits
double mask(long n)
{
	double beta = find_gold(n);
	if (0.5 < beta) return 1.0;
	return 0.0;
}

// The golden values themselves, or zero.
double gold(long n)
{
	double beta = find_gold(n);
	if (0.5 < beta) return beta;
	return 0.0;
}

// Allowed values
double allowed(long n)
{
	long idx = 1;
	int cnt = 0;
	while (cnt < n)
	{
		double beta = find_gold(idx);
		if (0.5 < beta) cnt++;
		idx++;
	}
	idx --;
	// printf("counted %ld is %ld\n", n, idx);
	return idx;
}

#if 0
int main(int argc, char* argv[])
{
	long nmax = 513;
	malloc_gold(nmax);

	printf("Mask OGF at 1/2 = %20.16g\n", OGF(mask, 0.5));
	printf("Gold OGF at 1/2 = %20.16g\n", OGF(gold, 0.5));
	printf("Allowed OGF at 1/2 = %20.16g\n", OGF(allowed, 0.5));
}
#endif

static double beta_disk(double re_q, double im_q, int itermax, double param)
{
	static bool init = false;
	if (not init)
	{
		long nmax = 513;
		malloc_gold(nmax);
		init = true;
	}

	COMPLEX zee = re_q + I * im_q;
	// COMPLEX og = COGF(mask, zee);
	COMPLEX og = CEGF(mask, zee);

	double faby = abs(og);
	double abz = abs(zee);
	double norm = abz * abz * exp(-abz);
	// printf("u %g %g %g %g \n", re_q, im_q, faby, norm);
	faby *= norm;
	return faby;

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
