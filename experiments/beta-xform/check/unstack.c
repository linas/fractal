/*
 * unstack.c
 *
 * Verify to unwrapped recursion relations for generalized
 * stretch-cut-stack map. The are the ones after noticing that
 * the unwrap simplifies even more.
 *
 * January 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// ==============================================================

// Arbitrary function
double nu(double x)
{
	if (x < 0.0) fprintf(stderr, "Error nu fail neg %g\n", x);
	if (1.0 < x) fprintf(stderr, "Error nu fail pos %g\n", x);

	return 1.0;
	// return x-0.5;

	// Bernoulli poly B_2
	// return x*x - x  + 1.0 / 6.0;

	// Bernoulli poly B_3
	// return x*x*x - 1.5*x*x  + 0.5*x;

	// Bernoulli poly B_4
	// return x*x*x*x - 2.0*x*x*x  + x*x - 1.0/30.0;
}

// ==============================================================
// Return endpoint iterate.
double t_k(double beta, int k)
{
	double tk = 1.0;
	for (int i=0; i<k; i++)
	{
		tk *= beta;
		if (1.0 <= tk) tk -= 1.0;
	}
	return tk;
}

// Return midpoint iterate digit b_k = d_k(1/2)=theta(beta t_k-1)
int b_k(double beta, int k)
{
	double tk = t_k(beta, k);
	if (beta*tk >= 1.0) return 1;
	return 0;
}

// Return the q_k constant
double q_k(double beta, int k)
{
	double sum = 0.0;
	double bej = 1.0;
	for (int j=0; j<=k; j++)
	{
		sum += b_k(beta, j) / bej;
		bej *= beta;
	}
	return sum * bej/beta;
}

// ==============================================================

double e_nk(double beta, double x, int n, int k);

double c_n(double beta, double x, int n)
{
	if (n < 1) fprintf(stderr, "Error fail index %d\n", n);
	if (x < 0.0) fprintf(stderr, "Error fail neg %d %g\n", n, x);
	if (1.0 < x) fprintf(stderr, "Error c_n fail pos n=%d x= %g\n", n, x);
	double sum = 0.0;
	for (int k=0; k<n; k++)
	{
		if (0 == b_k(beta, k)) continue;
		sum += e_nk(beta, x, n, k);
	}
	return sum;
}

double e_nk(double beta, double x, int n, int k)
{
	if (n <= k) fprintf(stderr, "Error enk fail index %d <= %d\n", n, k);
	double tk = t_k(beta, k);
	if (tk <= x) return 0.0;

	double ben = pow(beta, n);
	double bek = pow(beta, k);
	double arg = 1.0 + x/ben - tk/bek;
	double sum = nu(arg);

	double cnst = beta -1.0 - (beta*tk -1.0)/bek;
	double xen = x / ben;

	double bem = beta;
	for (int m=1; m < n-k; m++)
	{
		double arg = cnst + bem * xen;
		sum += bem * c_n(beta, arg, m);
		bem *= beta;
	}
	return sum/ben;
}

// ==============================================================

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		fprintf(stderr, "Usage: %s beta n\n", argv[0]);
		exit (1);
	}
	double beta = atof(argv[1]);
	int n = atoi(argv[2]);

#define PRINT_CEE
#ifdef PRINT_CEE

	int imax = 14;
	// double delta = 1.0 / ((double) imax);
	for (int i=0; i< imax; i++)
	{
		double x = (((double) i) + 0.5) / ((double) imax);
		double y = c_n(beta, x, n);

		printf("%d	%f	%g\n", i, x, y);

		fflush(stdout);
	}
#endif

// #define QCHECK
#ifdef QCHECK
	for (int k=0; k<= n; k++)
	{
		int bek = b_k(beta, k);
		double qk = q_k(beta, k);
		double bkp1 = pow(beta, k+1);
		double tk = t_k(beta, k+1);
		double rek = bkp1-qk;
		printf("%d %d qk=%f	b^{k+1}-qk=%f	t{k+1}=%f	%g\n",
			k, bek, qk, rek, tk, rek-tk);
	}
#endif

#if 0
	for (int k=0; k<= n; k++)
	{
		int bek = b_k(beta, k);
		double tk = t_k(beta, k);
		double rek = beta - 1.0 + (1.0 -tk - bek) * pow(beta, -k);
		printf("%d	%d	%f	%f\n", k, bek, tk, rek);
	}
#endif

#if 0
	for (int k=0; k<= n; k++)
	{
		int bek = b_k(beta, k);
		double tk = t_k(beta, k);
		double rek = (bek+tk)/beta;
		printf("%d	%d	%f	%f\n", k, bek, tk, rek);
	}
#endif
}
