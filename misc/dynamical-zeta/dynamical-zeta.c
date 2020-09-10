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

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

complex double dyn_zeta_bern(double x, double s)
{
	complex double sum = 0.0;
	for (int n=1; n<56; n++)
	{
		double term = pow (n, -s);
		complex double clock = cexp( 2.0* M_PI * I * x);
		sum += term * clock;
		x = 2*x;
		if (1.0 < x) x -= 1.0;	
	}
	return sum;
}

int main(int argc, char * argv[])
{
	if (argc < 1)
	{
		fprintf(stderr, "Usage: %s ess\n", argv[0]);
		exit(1);
	}

	double ess = atof(argv[1]);
	int WIDTH = 1000;
	printf("#\n# s=%g\n#\n", ess);
	for (int i=0; i<WIDTH; i++)
	{
		double x = ((double) i) / (double) WIDTH;
		complex double y = dyn_zeta_bern(x, ess);
		printf("%d	%g	%g	%g\n", i, x, creal(y), cimag(y));
	}
}
