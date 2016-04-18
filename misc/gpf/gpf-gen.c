/*
 * Generating functions for greatest prime factors.
 *
 * April 2016
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gpf.h>
#include "gpf-gen-bignum.h"

/*
 * Ordinary generating function for the greatest common factor.
 */
double gpf_ordinary(double x)
{
	double sum = 0;
	double xn = x;

	if (x < 1.0e-16) return x;

	for (int n=1; ; n++)
	{
		sum += gpf(n) * xn;
		xn *= x;
		if (n*xn < 1.0e-16*sum) break;
	}

	return sum;
}

/*
 * Exponential generating function for the greatest common factor.
 */
double gpf_exponential(double x)
{
	double sum = 0;
	double xn = x;

	if (x < 1.0e-16) return x;

	for (int n=1; ; n++)
	{
		sum += gpf(n) * xn;
		xn *= x / ((double) n);

		if (n*xn < 1.0e-16*sum) break;
	}

	return sum;
}

double gpf_bignum_exponential(double x, double theta)
{
	cpx_t sum, z;
	cpx_init(sum);
	cpx_init(z);

	theta *= 2.0 * M_PI;
	cpx_set_d(z, x*cos(theta), x*sin(theta));

	cpx_gpf_exponential(sum, z, 20);

	mpf_t val;
	mpf_init(val);
	cpx_abs(val, sum);
	double rv = mpf_get_d(val);
	return rv;
}

int main(int argc, char* argv[])
{
#ifdef ORD
	for (double x=0.0; x< 1.0; x+= 0.002)
	{
		double y = gpf_ordinary(x);
		double z = gpf_exponential(x);
		printf("%g\t%g\t%g\n", x, y, z);
	}
#endif
#ifdef EXPO
	for (double x=0.0; x< 675.0; x+= 0.5)
	{
		double y = gpf_exponential(x);
		double z = y * exp(-x);
		printf("%g\t%g\t%g\n", x, y, z);
		fflush(stdout);
	}
#endif
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s <r>\n", argv[0]);
		exit(1);
	}
	double dom = atof(argv[1]);
	printf("#\n# Max = %g\n#\n", dom);
	for (double x=0.0; x< dom; x+= 0.001*dom)
	{
		double w0 = gpf_bignum_exponential(x, 0.0);
		double w1_2 = gpf_bignum_exponential(x, 1.0/2.0);
		double w1_3 = gpf_bignum_exponential(x, 1.0/3.0);
		double w1_4 = gpf_bignum_exponential(x, 1.0/4.0);
		double w1_5 = gpf_bignum_exponential(x, 1.0/5.0);
		double w1_6 = gpf_bignum_exponential(x, 1.0/6.0);
		printf("%g\t%g\t%g\t%g\t%g\t%g\t%g\n", x, w0, w1_2, w1_3, w1_4, w1_5, w1_6);
		fflush(stdout);
	}
}
