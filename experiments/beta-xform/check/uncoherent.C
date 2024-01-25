/*
 * uncoherent.C
 *
 * Identical to unstack.c except a complex z is introduced into the sums.
 *
 * January 2024
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "unref.c"
#include "unutil.c"

#define COMPLEX std::complex<double>

// ==============================================================

// Arbitrary function
double nu(double x)
{
	if (x < 0.0) fprintf(stderr, "Error nu fail neg %g\n", x);
	if (1.0 < x) fprintf(stderr, "Error nu fail pos %g\n", x);

	// return 1.0;
	return x-0.5;
	// return x - 0.5 + 0.08684;  // appropriate for beta=1.6

	// Bernoulli poly B_2
	// The result is senstive to this being B_2.
	// Being able to integrate to exactly zero is important.
	// return x*x - x  + 1.0 / 6.0;
	// return x*x - x  + 0.16666;

	// Bernoulli poly B_3
	// return x*x*x - 1.5*x*x  + 0.5*x;

	// Bernoulli poly B_4
	// return x*x*x*x - 2.0*x*x*x  + x*x - 1.0/30.0;
}

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

	if (0 == k) return sum / ben;

	double cnst = beta -1.0 - (beta*tk -1.0)/bek;
	double xen = x / ben;

	double bem = beta;
	COMPLEX blem = blam;
	for (int m=1; m < n-k; m++)
	{
		double arg = cnst + bem * xen;
		sum += blem * cz_n(beta, blam, arg, m);
		beta *= beta;
		blem *= blam;
	}
	COMPLEX blan = pow(blam, n);
	return sum / blan;
}

COMPLEX hz_nk(double beta, COMPLEX blam, double x, int n, int k)
{
	if (n < k) fprintf(stderr, "Error hnk fail index %d <= %d\n", n, k);

	if (n == 0 && k == 0) return nu(x);
	if (k == 0) return 0.0;

	double tk = t_k(beta, k);
	double bek = pow(beta, k);
	if (n == k)
	{
		double arg =  1.0 + (x-tk)/bek;
		return nu(arg) / bek;
	}

	double arg = beta - 1.0 + (x +1.0 - beta*tk)/bek;
	COMPLEX blek = pow(blam, k);
	return cz_n(beta, blam, arg, n-k) / blek;
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
// #include "un.c"  // for unit testing only. Copy of unwrap.c w/o main()

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		fprintf(stderr, "Usage: %s beta n\n", argv[0]);
		exit (1);
	}
	double beta = atof(argv[1]);
	int n = atoi(argv[2]);

}
