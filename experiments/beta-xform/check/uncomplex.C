/*
 * uncomplex.C
 *
 * Identical to unstack.c except a complex z is introduced into the sums.
 *
 * January 2024
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "unutil.c"

#define COMPLEX std::complex<double>

// ==============================================================

COMPLEX ez_nk(double beta, COMPLEX blam, double x, int n, int k);

// Due to use in intermediate values, it can (commonly) happen that 1<x
COMPLEX cz_n(double beta, COMPLEX blam, double x, int n)
{
	if (x < 0.0) fprintf(stderr, "Error fail neg %d %g\n", n, x);
	if (0 == n) return 0.0;

	COMPLEX sum = 0.0;
	for (int k=0; k<n; k++)
	{
		if (0 == b_k(beta, k)) continue;
		sum += ez_nk(beta, blam, x, n, k);
	}
	return sum;
}

COMPLEX ez_nk(double beta, COMPLEX blam, double x, int n, int k)
{
	if (n <= k) fprintf(stderr, "Error enk fail index %d <= %d\n", n, k);
	double tk = t_k(beta, k);
	double ben = pow(beta, n);
	double bek = pow(beta, k);
	double arg = 1.0 + x/ben - tk/bek;
	COMPLEX sum = nu(arg);

	COMPLEX bln = pow(blam, n);
	if (0 == k) return sum / bln;

	double cnst = beta -1.0 - (beta*tk -1.0)/bek;
	double xen = x / ben;

	double bem = beta;
	COMPLEX blem = blam;
	for (int m=1; m < n-k; m++)
	{
		double arg = cnst + bem * xen;
		sum += blem * cz_n(beta, blam, arg, m);
		bem *= beta;
		blem *= blam;
	}
	return sum / bln;
}

COMPLEX hz_nk(double beta, COMPLEX blam, double x, int n, int k)
{
	if (n < k) fprintf(stderr, "Error hnk fail index %d <= %d\n", n, k);

	if (n == 0 && k == 0) return nu(x);
	if (k == 0) return 0.0;

	double tk = t_k(beta, k);
	double bek = pow(beta, k);
	COMPLEX blk = pow(blam, k);
	if (n == k)
	{
		double arg =  1.0 + (x-tk)/bek;
		return nu(arg) / blk;
	}

	double arg = beta - 1.0 + (x +1.0 - beta*tk)/bek;
	return cz_n(beta, blam, arg, n-k) / blk;
}

COMPLEX nuz_n(double beta, COMPLEX blam, double x, int n)
{
	COMPLEX sum = cz_n(beta, blam, x, n);
	for (int k=0; k<= n; k++)
	{
		double tk = t_k(beta, k);
		if (tk <= x) continue;
		sum += hz_nk(beta, blam, x, n, k);
	}
	return sum;
}

// ==============================================================
