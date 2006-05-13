
/*
 * plouffe.c
 *
 * Explore some of the sums given by Simon Plouffe
 *
 * Linas May 2006
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "bernoulli.h"
#include "binomial.h"
#include "harmonic.h"

static double 
ess_tee_k (int k, double x, int cyn)
{
	int n;

	double sum = 0.0;
	double ex = exp (x);
	double exn = ex;
	for (n=1; n<300; n++)
	{
		double term = pow (n, -k);
		term /= (exn + cyn);
		sum += term;

		if (term < 1.0e-16*sum) break;
		exn *= ex;
	}

	return sum;
}

double ess_k (int k, double x)
{
	return ess_tee_k (k,x,-1);
}
 
double tee_k (int k, double x)
{
	return ess_tee_k (k,x,1);
}

double eye_k (int k, double x)
{
	int j;

	x /= 2.0*M_PI;

	double cyn = 1.0;
	if ((k+1)%2) cyn = -1.0;

	double zt = 1.0 + cyn * pow (x, k-1);
	zt *= zetam1 (k) + 1.0;

	double sum = 0.0;
	double xn = 1.0;
	for (j=0; j<=(k+1)/2; j++)
	{
		double term = bernoulli (2*j) * bernoulli (k+1-2*j);
		term /= factorial (2*j) * factorial (k+1-2*j);
		term *= xn;

		if (j%2) 
		{
			sum -= term;
		}
		else
		{
			sum += term;
		}

		xn *= x*x;
	}

	sum *= cyn * pow (2.0*M_PI, k);

	sum += zt;
	sum *= -0.5;

	return sum;
}

int test_tee_ident (void)
{
	int errcnt = 0;
	int k;
	for (k=3; k<50; k+=2)
	{
		double x;
		for (x=0.05; x<20; x+=0.12345)
		{
			double tee = tee_k (k,x);
			double sum = tee - ess_k(k,x) + 2.0*ess_k(k, 2.0*x);
			sum = fabs (sum/tee);
			if (sum > 2.0e-14) 
			{
				printf ("Error: tee ident failed, k=%d x=%g sum=%g\n", k, x, sum);
				errcnt ++;
			}
		}
	}

	return errcnt;
}
 
main (int argc, char * argv[])
{
	if (argc < 3)
	{
		fprintf (stderr, "Usage: %s <k> <x>\n", argv[0]);
		exit (1);
	}
	int k = atoi(argv[1]);
	double x = atof(argv[2]);

	test_tee_ident ();
}
