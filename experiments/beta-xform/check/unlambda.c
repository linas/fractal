/*
 * lambda.c
 *
 * Version of unstack.c with lambda inserted into it.
 *
 * January 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "unutil.c"

// ==============================================================

double el_nk(double beta, double blam, double x, int n, int k);

// Due to use in intermediate values, it can (commonly) happen that 1<x
double cl_n(double beta, double blam, double x, int n)
{
	if (x < 0.0) fprintf(stderr, "Error fail neg %d %g\n", n, x);
	if (0 == n) return 0.0;

	double sum = 0.0;
	for (int k=0; k<n; k++)
	{
		if (0 == b_k(beta, k)) continue;
		sum += el_nk(beta, blam, x, n, k);
	}
	return sum;
}

double el_nk(double beta, double blam, double x, int n, int k)
{
	if (n <= k) fprintf(stderr, "Error enk fail index %d <= %d\n", n, k);
	double tk = t_k(beta, k);
	double ben = pow(beta, n);
	double bek = pow(beta, k);
	double arg = 1.0 + x/ben - tk/bek;
	double sum = nu(arg);

	double bln = pow(blam, n);
	if (0 == k) return sum / bln;

	double cnst = beta -1.0 - (beta*tk -1.0)/bek;
	double xen = x / ben;

	double bem = beta;
	double blm = blam;
	for (int m=1; m < n-k; m++)
	{
		double arg = cnst + bem * xen;
		sum += blm * cl_n(beta, blam, arg, m);
		bem *= beta;
		blm *= blam;
	}
	return sum / bln;
}

double hl_nk(double beta, double blam, double x, int n, int k)
{
	if (n < k) fprintf(stderr, "Error hnk fail index %d <= %d\n", n, k);

	if (n == 0 && k == 0) return nu(x);
	if (k == 0) return 0.0;

	double tk = t_k(beta, k);
	double bek = pow(beta, k);
	double blk = pow(blam, k);
	if (n == k)
	{
		double arg =  1.0 + (x-tk)/bek;
		return nu(arg) / blk;
	}

	double arg = beta - 1.0 + (x +1.0 - beta*tk)/bek;
	return cl_n(beta, blam, arg, n-k) / blk;
}

double nul_n(double beta, double blam, double x, int n)
{
	double sum = cl_n(beta, blam, x, n);
	for (int k=0; k<= n; k++)
	{
		double tk = t_k(beta, k);
		if (tk <= x) continue;
		sum += hl_nk(beta, blam, x, n, k);
	}
	return sum;
}

// ==============================================================
