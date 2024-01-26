/*
 * lambda.c
 *
 * Version of unstack.c with lambda inserted into it.
 *
 * January 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "unref.c"
#include "unutil.c"

// ==============================================================

// Arbitrary function
double nu(double x)
{
	if (x < 0.0) fprintf(stderr, "Error nu fail neg %g\n", x);
	if (1.0 < x) fprintf(stderr, "Error nu fail pos %g\n", x);

	// return 1.0;
	return x-0.5;
	// return x - 0.5 + 0.08684;  // appropriate for beta=1.6

	// Bernoulli poly B_2
	// The result is senstive to this being B_2.
	// Being able to integrate to exactly zero is important.
	// return x*x - x  + 1.0 / 6.0;
	// return x*x - x  + 0.16666;

	// Bernoulli poly B_3
	// return x*x*x - 1.5*x*x  + 0.5*x;

	// Bernoulli poly B_4
	// return x*x*x*x - 2.0*x*x*x  + x*x - 1.0/30.0;
}

// ==============================================================

double e_nk(double beta, double blam, double x, int n, int k);

// Due to use in intermediate values, it can (commonly) happen that 1<x
double c_n(double beta, double blam, double x, int n)
{
	if (x < 0.0) fprintf(stderr, "Error fail neg %d %g\n", n, x);
	if (0 == n) return 0.0;

	double sum = 0.0;
	for (int k=0; k<n; k++)
	{
		if (0 == b_k(beta, k)) continue;
		sum += e_nk(beta, blam, x, n, k);
	}
	return sum;
}

double e_nk(double beta, double blam, double x, int n, int k)
{
	if (n <= k) fprintf(stderr, "Error enk fail index %d <= %d\n", n, k);
	double tk = t_k(beta, k);
	double ben = pow(beta, n);
	double bek = pow(beta, k);
	double arg = 1.0 + x/ben - tk/bek;
	double sum = nu(arg);

	double bln = pow(blam, n);
	if (0 == k) return sum / bln;

	double cnst = beta -1.0 - (beta*tk -1.0)/bek;
	double xen = x / ben;

	double bem = beta;
	double blm = blam;
	for (int m=1; m < n-k; m++)
	{
		double arg = cnst + bem * xen;
		sum += blm * c_n(beta, blam, arg, m);
		bem *= beta;
		blm *= blam;
	}
	return sum / bln;
}

double h_nk(double beta, double blam, double x, int n, int k)
{
	if (n < k) fprintf(stderr, "Error hnk fail index %d <= %d\n", n, k);

	if (n == 0 && k == 0) return nu(x);
	if (k == 0) return 0.0;

	double tk = t_k(beta, k);
	double bek = pow(beta, k);
	double blk = pow(blam, k);
	if (n == k)
	{
		double arg =  1.0 + (x-tk)/bek;
		return nu(arg) / blk;
	}

	double arg = beta - 1.0 + (x +1.0 - beta*tk)/bek;
	return c_n(beta, blam, arg, n-k) / blk;
}

double nu_n(double beta, double blam, double x, int n)
{
	double sum = c_n(beta, blam, x, n);
	for (int k=0; k<= n; k++)
	{
		double tk = t_k(beta, k);
		if (tk <= x) continue;
		sum += h_nk(beta, blam, x, n, k);
	}
	return sum;
}

// ==============================================================
// #include "un.c"  // for unit testing only. Copy of unwrap.c w/o main()

int main(int argc, char* argv[])
{
	if (argc != 4)
	{
		fprintf(stderr, "Usage: %s beta lambda n\n", argv[0]);
		exit (1);
	}
	double beta = atof(argv[1]);
	double lambda = atof(argv[2]);
	int n = atoi(argv[3]);

	double blam = beta * lambda;

#define PRINT_NU
#ifdef PRINT_NU

	double scale = lambda * beta;
	scale = lambda;
	double scan = pow(scale, n);
	int imax = 814;
	for (int i=0; i< imax; i++)
	{
		double x = (((double) i) + 0.5) / ((double) imax);

		double y = gp_invar(beta, x);
		printf("%d	%f	%f", i, x, y);

#define NIT 6
		double plm = scan;
		for (int j=0; j<NIT; j++)
		{
			double y = nu_n(beta, blam, x, n+j);
			y *= plm;
			plm *= scale;
			printf("	%f", y);
		}
		printf("\n");
		fflush(stdout);
	}
#endif
}
