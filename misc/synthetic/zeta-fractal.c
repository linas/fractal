
#include <complex.h>
#include <stdio.h>

complex phi(int p, int k, int r, complex lambda, double x)
{
	complex sum = 0.0;

	return sum;
}

int main (int argc, char *argv[])
{
	double x = 0.3;
	complex la = 0.5 + I*4.3;
	complex ph = phi(2,0,1,la, x);

	double re = creal(ph);
	double im = cimag(ph);
	printf ("its %g %g\n", re, im);

	return 0;
}
