/*
 * tsum.c
 *
 * Validate that tsum and csum are right.
 * January 2024
 */

#include <stdio.h>

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

double csum(double beta, int n)
{
	double sum=0.0;
	double tn = 1.0;
	double bn = beta;
	for (int i=0; i<n; i++)
	{
		if (beta*tn >= 1.0) sum += 1.0 / bn;

		tn *= beta;
		if (1.0 <= tn) tn -= 1.0;

		bn *= beta;
	}
	return 1.0 - sum;
}

int main (int argc, char* argv[])
{

	double beta = 1.6;
	for (int n=0; n<20; n++)
	{
		double ts = tsum(beta, n);
		double cs = csum(beta, n);
		printf("its %d  %g %g\n", n, ts, cs);
	}
}
