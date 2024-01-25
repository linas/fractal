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

// Due to use in intermediate values, it can (commonly) happen that 1<x
double c_n(double beta, double x, int n)
{
	if (x < 0.0) fprintf(stderr, "Error fail neg %d %g\n", n, x);
	if (0 == n) return 0.0;

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
	double ben = pow(beta, n);
	double bek = pow(beta, k);
	double arg = 1.0 + x/ben - tk/bek;
	double sum = nu(arg);

	if (0 == k) return sum / ben;

	double cnst = beta -1.0 - (beta*tk -1.0)/bek;
	double xen = x / ben;

	double bem = beta;
	for (int m=1; m < n-k; m++)
	{
		double arg = cnst + bem * xen;
		sum += bem * c_n(beta, arg, m);
		bem *= beta;
	}
	return sum / ben;
}

double h_nk(double beta, double x, int n, int k)
{
	if (n < k) fprintf(stderr, "Error hnk fail index %d <= %d\n", n, k);

	if (n == 0 && k == 0) return nu(x);
	if (k == 0) return 0.0;

	double tk = t_k(beta, k);
	double bek = pow(beta, k);
	if (n == k)
	{
		double arg =  1.0 + (x-tk)/bek;
		return nu(arg) / bek;
	}

	double arg = beta - 1.0 + (x +1.0 - beta*tk)/bek;
	return c_n(beta, arg, n-k) / bek;
}

double nu_n(double beta, double x, int n)
{
	double sum = c_n(beta, x, n);
	for (int k=0; k<= n; k++)
	{
		double tk = t_k(beta, k);
		if (tk <= x) continue;
		sum += h_nk(beta, x, n, k);
	}
	return sum;
}

// ==============================================================
// #include "un.c"  // for unit testing only. Copy of unwrap.c w/o main()

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

#define NIT 6
	double sum[NIT];
	for (int j=0; j<NIT; j++) sum[j] = 0.0;

	// double lambda = 1.0;
	double lambda = 1.0 / beta;
	// double lambda = 1.0 / (beta*beta);
	// double lambda = 1.0 / (beta*beta*beta);
	double lamn = pow(lambda, n);

	int imax = 814;
	double delta = 1.0 / ((double) imax);
	for (int i=0; i< imax; i++)
	{
		double x = (((double) i) + 0.5) / ((double) imax);

		double y = gp_invar(beta, x);
		// double y = gp_n1(beta, x);
		// double y = gp_n2(beta, x);
		// double y = gp_n3(beta, x);
		// double y = gp_quad_n1(beta, x);
		printf("%d	%f	%f", i, x, y);

		double lscale = lamn;
		for (int j=0; j<NIT; j++)
		{
			// double y = c_n(beta, x, n+j);
			double y = nu_n(beta, x, n+j);
			y /= lscale;
			lscale *= lambda;

			sum[j] += y * delta;
			printf("	%f", y);
		}
		printf("\n");
		fflush(stdout);
	}

	printf("#\n# ");
	for (int j=0; j<NIT; j++)
		printf(" %g", sum[j]);
	printf("\n#\n");
#endif

// #define UNIT_TEST
#ifdef UNIT_TEST

	// Compare new code to old code from unwrap.c. Everything passes.
	int imax = 14;
	for (int i=0; i< imax; i++)
	{
		double x = (((double) i) + 0.5) / ((double) imax);
		for (int n=0; n<6; n++)
		{
			// Test passes; although nu overflows.
			// Overflow is normal cause this is artificial.
			for (int k=1; k<n; k++)
			{
				double oy = esum_n_k(beta, n, k, x);
				double ny = e_nk(beta, x, n, k);
				printf("%d %f %d %d  %f %f   %g\n", i, x, n, k, oy, ny, oy-ny);
				if (1.0e-12 < fabs(oy-ny)) printf("------------FAIL\n");
			}

			// Test passes; although nu overflows.
			// Overflow is normal cause this is artificial.
			for (int k=1; k<n; k++)
			{
				double oy = hsum_n_k(beta, n, k, x);
				double ny = h_nk(beta, x, n, k);
				printf("%d %f %d %d  %f %f   %g\n", i, x, n, k, oy, ny, oy-ny);
				if (1.0e-12 < fabs(oy-ny)) printf("------------FAIL\n");
			}

			// Test passes, no overflows, everything is OK.
			double oy = cee_n(beta, n, x);
			double ny = c_n(beta, x, n);
			printf("%d %f %d  %f %f   %g\n", i, x, n, oy, ny, oy-ny);
			if (1.0e-12 < fabs(oy-ny)) printf("------------FAIL\n");

			// Test passes, no overflows, everything is OK.
			double oy = unu_n(beta, n, x);
			double ny = nu_n(beta, x, n);
			printf("%d %f %d  %f %f   %g\n", i, x, n, oy, ny, oy-ny);
			if (1.0e-12 < fabs(oy-ny)) printf("------------FAIL\n");

			// Test basic h_n1 -- this passes; although nu overflows.
			// Overflow is normal cause this is artificial.
			if (n<2) continue;
			double oy = hsum_n_1(beta, n, x);
			// double oy = hsum_n_k(beta, n, 1, x);
			// double ny = c_n(beta, (1.0+x)/beta, n-1)/beta;
			double ny = h_nk(beta, x, n, 1);
			printf("%d %f %d  %f %f   %g\n", i, x, n, oy, ny, oy-ny);
			if (1.0e-12 < fabs(oy-ny)) printf("------------FAIL\n");
		}
		printf("\n");
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
