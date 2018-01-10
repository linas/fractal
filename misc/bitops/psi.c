/*
 * psi.c
 *
 * Attempt at orthonormal fns.
 * January 2018
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double psi_0(double x, double K)
{
	if (x<K) return 1.0/sqrt(K);
	return 0.0;
}

double psi_1(double x, double K)
{
	double step = 2.0*K*K-K;
	if (x<step) return 1.0/sqrt(2.0*step);
	if (x<K) return -1.0/sqrt(2.0*(K-step));
	return 0.0;
}

double xfer(double x, double K, double (*fun)(double, double))
{
	if (K<x) return 0.0;
	double res = x / (2.0*K);
	double elf = fun(res, K) + fun(0.5+res, K); 
	elf /= 2.0*K;
	return elf;
}

double part(double x, double K, double (*fun)(double, double), int which)
{
	if (K<x) return 0.0;
	double res = x / (2.0*K);
	double elf = fun(0.5*which+res, K); 
	elf /= 2.0*K;
	return elf;
}

int main(int argc, char* argv[])
{
	double K = atof(argv[1]);

#define NPTS 201
	double s = 0.0;
	for (int i=0; i< NPTS; i++)
	{
		double x = ((double) i + 0.5) / ((double) NPTS);
		double (*f)(double, double) = psi_1;
		double ef = f(x, K);
		double p0 = part(x, K, f, 0);
		double p1 = part(x, K, f, 1);
		double y = xfer(x, K, f);
		s += y / ((double) NPTS);
		printf("%d	%g	%g %g	%g	%g	%g\n", i, x, ef, p0, p1, y, s);
	}
}
