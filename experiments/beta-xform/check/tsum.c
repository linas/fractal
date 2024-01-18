/*
 * tsum.c
 *
 * Validate that tsum and csum are right.
 * Yes, it does check out correctly. Yay!
 * January 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double tsum(double beta, int n)
{
	double sum=0.0;
	double tn = 1.0;
	double bn = 1.0;
	for (int i=0; i<=n; i++)
	{
		sum += tn / bn;

		tn *= beta;
		if (1.0 <= tn) tn -= 1.0;

		bn *= beta;
	}
	return sum;
}

double c_sub_n(double beta, int n)
{
	double sum=0.0;
	double tn = 1.0;
	double bn = beta;
	for (int i=0; i<=n; i++)
	{
		if (beta*tn >= 1.0) sum += 1.0 / bn;

		tn *= beta;
		if (1.0 <= tn) tn -= 1.0;

		bn *= beta;
	}
	return 1.0 - sum;
}

double t_sub_n(double beta, int n)
{
	double tn = 1.0;
	for (int i=0; i<n; i++)
	{
		tn *= beta;
		if (1.0 <= tn) tn -= 1.0;
	}
	return tn;
}

double csum(double beta, int n)
{
	double sum=0.0;
	for (int i=0; i<n; i++)
		sum += c_sub_n(beta, i);

	return sum;
}

double tnksum(double beta, int n)
{
printf("----\n");
	double sum=0.0;
	double tn = 1.0;
	double bn = 1.0;
	double dkm1 = 0.0;
	for (int i=0; i<=n; i++)
	{
		double term = tn + dkm1 * (n-i+1);
		double frac = term / bn;
		sum += frac;
printf("k=%d tk=%g  dk=%g  term=%g  frac=%g  sum=%g\n",
i, tn, dkm1, term, frac, sum);
		if (beta*tn >= 1.0)
			dkm1 = 1.0;
		else
			dkm1 = 0.0;

		tn *= beta;
		if (1.0 <= tn) tn -= 1.0;

		bn *= beta;
	}
	sum -= n;
printf("Final n=%d sum = %g\n", n, sum);
	return sum;
}

int main (int argc, char* argv[])
{
	double beta = atof(argv[1]);

#if 0
	for (int n=0; n<20; n++)
	{
		double ts = tsum(beta, n);
		double cs = csum(beta, n);
		double cn = c_sub_n(beta, n);
		double norm = ts - cs;
		printf("its %d cn=%g	%g %g 	%g\n", n, cn, ts, cs, norm-1.0);
	}
#endif

// #define WTF
#ifdef WTF
	for (int n=0; n<20; n++)
	{
		double ts = tsum(beta, n);
		double cs = csum(beta, n);
		double cn = c_sub_n(beta, n);
		double tnk = tnksum(beta, n);
		printf("its %d cn=%g	tn=%g  %g 	%g\n", n, cn, ts, cs, tnk);
	}
#endif

	for (int n=1; n<20; n++)
	{
		double cn = c_sub_n(beta, n-1);
		double tn = t_sub_n(beta, n);
		double cs = pow(beta, n) * cn;
		double diff = tn - cs;
		printf("its %d cn=%g	tn=%g  %g  %g\n", n, cn, tn, cs, diff);
	}
}
