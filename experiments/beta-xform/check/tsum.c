/*
 * tsum.c
 *
 * Validate that tsum and csum are right.
 * Yes, it does check out correctly. Yay!
 * January 2024
 */

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

double csum(double beta, int n)
{
	double sum=0.0;
	for (int i=0; i<n; i++)
		sum += c_sub_n(beta, i);

	return sum;
}

int main (int argc, char* argv[])
{
	double beta = atof(argv[1]);
	for (int n=0; n<20; n++)
	{
		double ts = tsum(beta, n);
		double cs = csum(beta, n);
		double cn = c_sub_n(beta, n);
		double norm = ts - cs;
		printf("its %d cn=%g	%g %g 	%g\n", n, cn, ts, cs, norm-1.0);
	}
}
