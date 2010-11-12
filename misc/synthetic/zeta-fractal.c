
#include <complex.h>
#include <math.h>
#include <stdio.h>

complex phi(unsigned int p, int k, int r, complex lambda, double x)
{
	complex sum = 0.0;
	complex ln = 1.0;
	unsigned int pn = 1;

	while(1)
	{
		double arg = 2.0*M_PI*x*pn*(p*k+r);
		sum += ln * cexp(I*arg);
		ln *= lambda;
		pn *= p;
		if (cabs(ln) < 1.0e-12) break;
	}

	return sum;
}

complex remi (complex ess, double x)
{
	int p;
	complex sum = 0.0;

	for (p=2; p<20; p++)
	{
		complex lambda = cexp(-ess * log(p));
		complex ph = phi(p, 0, 1, lambda, x);
		sum += ph;
	}

	return sum;
}

int main (int argc, char *argv[])
{
	double x = 0.3;
	complex ess = 0.5 + I*4.3;

	for (x=0; x<=1.0; x+= 0.005)
	{
		complex ph = remi(ess, x);

		double re = creal(ph);
		double im = cimag(ph);
		printf ("%f\t%g\t%g\n", x, re, im);
	}

	return 0;
}
