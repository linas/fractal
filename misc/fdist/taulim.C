/*
 * taulim.C
 *
 * Explore limiting behaviour of the GKW eigenfuncs at zero.
 *
 * March 2010
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"

double triangle(double x)
{
	x -= floor(x);
	if (x < 0.5) return 2.0*x;
	return 2.0-2.0*x;
}

#if MONOSUM
double sum(int k, int n, int lmax)
{
	int l;

	double tk = pow(2.0, k-n);
	double acc = 0.0;
#define LMAX 10123
	for (l=0; l<lmax; l++)
	{
		double tlp1 = 2*l+1;

		double alpha = 1.0;
		acc += alpha * triangle (tlp1*tk) / tlp1;
	}

	acc *= n+1;
	return acc;
}
#endif

#if ALLSUM
double delta (int l)
{
	// return pow(0.5, l) / ((double) (2*l+1));
	double d = 1.0 / ((double) (2*l+1));
	d = pow (d, 1.1);
	return d;
}

double delta_sum(int lambda, int p)
{
	if (p>30) return delta(lambda);

	int tp = 1<<p;
	int m = 0;
	double acc =0.0;
	while(1)
	{
		double term = delta(lambda + m * tp);
		acc += term;
		m ++;
		if (fabs(term/acc) < 1.0e-10) break;
	}
	return acc;
}

double lambda_sum(int p)
{
	int tp = 1<<p;
	double otp = 1.0 / ((double) tp);

	double acc = 0.0;
	int lambda;
	int smax = p/2 -1;
	for (lambda = 0; lambda < smax; lambda++)
	{
		double tlp1 = 2*lambda+1;
		acc += triangle(tlp1*otp) * delta_sum(lambda, p);
	}
	return acc;
}

double w_sum(int n, double w)
{
	int k;
	double acc = 0.0;
	double wk = 1.0;
	for (k=0; k<n-1; k++)
	{
		acc += wk * lambda_sum(n-k);
		wk *= 0.5*w;
	}

	acc *= n+1;
	acc *= w-1.0;
	acc *= 0.25;

	return acc;
}
#endif

double tw(double w, int l, int n)
{
	int k;
	double tlp1 = 2*l+1;
	int tpn = 1<<n;
	double tn = 1.0 / ((double) tpn);
	double wh = 1.0;
	double acc = 0.0;
	for (k=0; k<n; k++)
	{
		acc += wh * triangle(tlp1 * tn);

		tn *= 2.0;
		wh *= 0.5 * w;
	}

	return acc;
}

double takagi_w(double w, double x)
{
	int k;
	double tn = 1.0;
	double wh = 1.0;
	double acc = 0.0;
	for (k=0; k<30; k++)
	{
		acc += wh * triangle(tn * x);

		tn *= 2.0;
		wh *= w;
	}

	return acc;
}

double dual(double w, int l, int n)
{

	double v = tw(w,l,n);
	v /= 2*l+1;
	v *= n+1;
	// v *= 0.25 * (w-1.0);
	return v;
}

double dual_sum(double w, double x)
{
	int l, n;

	double acc = 0.0;
	for (n=2; n<20; n++)
	{
		l = 3 * (1<<(n-1));
		double tlp1 = 2*l+1;

		double alpha = 1.0 / dual(w,l,n);
		acc += alpha * takagi_w (w, tlp1*x) / tlp1;
	}

	return acc;
}

main(int argc, char * argv[])
{
	int i, n, l;
	int npts;
	int k = 2;
	double delta;

	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s blah\n", argv[0]);
		exit(1);
	}

	// k = atoi(argv[1]);
	// int lmax = atoi(argv[2]);
	double w = atof(argv[1]);

#if COEFFS_GRAPH
	// for (n=1; n<30; n++)
	for (l=0; l<600; l++)
	{
		// double y = sum(k,n,lmax);
		// double y = w_sum(n,w);

		printf ("%d", l);
		for (n=2; n<16; n++)
		{
			double y = dual(w,l, n);
			printf ("\t%g", y);
		}
		printf("\n");
		fflush(stdout);
		// printf ("its %d sum=%g\n", n, y);
	}
#endif

	npts = 600;

	ContinuedFraction f;

	delta = 1.0 / ((double) npts);
	double yprev = 0.0;
	for (i=0; i<npts; i++)
	{
		double x = i*delta;
		f.SetReal(x);
		double qx = f.ToFarey();
		double y = dual_sum(w, qx);
		double dy = (y-yprev) /delta;
		printf("%d	%g	%g	%g\n", i, x, y, dy);
		yprev = y;
	}
}
