/*
 * zeta-fractal.c
 *
 * Experimental sum. Not tractable, and creted due to 
 * a mistake, in the first place. Dead end.
 * Linas Vepstas november 2010
 */
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

complex phi(unsigned int p, int k, int r, complex lambda, double x)
{
	complex sum = 0.0;
	complex lamn = lambda;
	unsigned int pn = p;

	while(1)
	{
		double arg = 2.0*M_PI*x*pn*(p*k+r);
		sum += lamn * cexp(I*arg);
		lamn *= lambda;
		pn *= p;
		if (cabs(lamn) < 1.0e-16) break;
	}

	return sum;
}

#define NTERMS 2020

complex remi (complex ess, double x)
{
	int p;
	complex sum = 0.0;

	for (p=2; p<NTERMS; p++)
	{
		complex lambda = cexp(-ess * log(p));
		complex ph = phi(p, 0, 1, lambda, x);
		sum += ph;
	}

	return sum;
}

int main (int argc, char *argv[])
{
	double x, sigma, tau;
	complex ess = 0.5 + I*14.13472514173;
	ess = 0.5 + I*21.022039638771;
	ess = 0.5 + I*25.0108575801;

	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s <sigma> <tau>\n", argv[0]);
		exit(1);
	}

	sigma = atof(argv[1]);
	tau = atof(argv[2]);
	printf("#\n# s = %f+i %f\n#\n", sigma, tau);

	ess = sigma + I*tau;

	for (x=0; x<=1.0; x+= 0.0005)
	{
		complex ph = remi(ess, x);

		double re = creal(ph);
		double im = cimag(ph);
		printf ("%f\t%g\t%g\n", x, re, im);
	}

	return 0;
}
