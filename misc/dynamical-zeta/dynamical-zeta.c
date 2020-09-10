/*
 * dynamical-zeta.c
 *
 * Quick and cheesy dynamical zeta defined as
 *
 * zeta(X;s) = sum_{n=1} n^{-s} exp 2i\pi x_n
 *
 * for some sequence X={x_n}
 *
 * We're gonna start with the Bernoulli sequence.
 *
 * Linas Vepstas September 2020
 */

#include <math.h>

double dyn_zeta_bern(double x, double s)
{
	double sum = 0.0;
	for (int n=1; n<56; n++)
	{
		double term = pow (n, s);
	}
	return sum;
}
