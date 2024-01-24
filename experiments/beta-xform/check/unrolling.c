/*
 * unrolling.c
 *
 * Verify to unrolled recursion relations for generalized
 * stretch-cut-stack map. Result: they are correct, and do
 * reproduce older results. But also obsolete; the unwrap.c
 * version is newer, better.
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
// Recursive forms

// Forward decl
double g_n_1(double beta, int n, double x);
double gsum_n_1(double beta, int n, double x);

// Return the g_n_k constant from the "generalized stretch-cut-stack"
// section of paper. This is computed using the recursing formula.
double g_n_k(double beta, int n, int k, double x)
{
	if (n < k) return 0.0;
	if (n < 0) { fprintf(stderr, "Error: negative n\n"); return 0.0; }
	if (0 == k) return nu(x);

	if (1 == k) return g_n_1(beta, n, x);
	// if (1 == k) return gsum_n_1(beta, n, x);

	// Recurse
	double bkm1 = b_k(beta, k-1);
	double arg = (x + bkm1)/ beta;
	double guh = g_n_k(beta, n-1, k-1, arg);
	return guh / beta;
}

// Return the f_n_k constant from the "generalized stretch-cut-stack"
// section of paper. This is computed using the recursing formula.
double f_n_k(double beta, int n, int k, double x)
{
	if (n <= k) return 0.0;
	if (k < 0) fprintf(stderr, "Error can't negative k\n");
	double arg = x / beta;
	double sum = f_n_k(beta, n-1, k, arg) + g_n_k(beta, n-1, k, arg);
	return sum / beta;
}

// Return the g_n_1 constant from the "generalized stretch-cut-stack"
// section of paper. This is computed with the recursive formula.
double g_n_1(double beta, int n, double x)
{
	double arg = (x + 1.0) / beta;
	double sum = nu(arg);
	for (int k=0; k<n-1; k++)
	{
		if (b_k(beta, k))
			sum += f_n_k(beta, n-1, k, arg);
	}
	return sum / beta;
}

// ==============================================================
// Series summations

double gsum_n_k(double beta, int n, int k, double x);

// Return the g_n_1 constant from the "generalized stretch-cut-stack"
// section of paper. This is computed with the series sum.
double gsum_n_1(double beta, int n, double x)
{
	double bej = 1.0 / beta;
	double arg = (x + 1.0) / beta;

	double sum = bej * nu(arg);
	for (int j=0; j<n-1; j++)
	{
		bej /= beta;
		arg /= beta;

		double bitso = 0.0;
		for (int k=0; k< n-j-1; k++)
		{
			if (0 == b_k(beta, k)) continue;
			bitso += gsum_n_k(beta, n-j-2, k, arg);
		}
		sum += bej * bitso;
	}

	return sum;
}

// Return the g_n_k constant from the "generalized stretch-cut-stack"
// section of paper. This is computed using the summation formula.
double gsum_n_k(double beta, int n, int k, double x)
{
	if (n < k) return 0.0;
	if (0 == k) return nu(x);

	// Loop.
	double arg = 0.0;
	double bei = 1.0;
	for (int i=1; i<k; i++)
	{
		bei /= beta;
		arg += b_k(beta, i) * bei;
	}
	double bek = bei;
	arg += x * bek;
	return gsum_n_1(beta, n-k+1, arg) * bek;
}

// Return the f_n_k constant from the "generalized stretch-cut-stack"
// section of paper. This is computed using the summation formula.
double fsum_n_k(double beta, int n, int k, double x)
{
	if (n <= k) return 0.0;
	if (k < 0) fprintf(stderr, "Error can't negative k\n");
	double sum = 0.0;
	double bej = 1.0;
	for (int j=1; j <= n-k; j++)
	{
		bej /= beta;
		sum += bej * gsum_n_k(beta, n-j, k, x*bej);
	}
	return sum;
}

// ==============================================================
// Density expressions

// Iterated density
double rho_n(double beta, int n, double x)
{
	double sum = 0.0;
	for (int k=0; k <=n; k++)
	{
		double tk = t_k(beta, k);
		if (tk > x)
			sum += gsum_n_k(beta, n, k, x);
	}
	return sum;
}

// Forward decl
double nu_n(double beta, int n, double x);

// Iterated density
double dee_n(double beta, int n, double x)
{
	double sum = 0.0;
	for (int k=0; k <n; k++)
	{
		sum += nu_n(beta, k, x);
		if (b_k(beta, k))
			sum -= fsum_n_k(beta, n, k, x);
	}
	return sum;
}

// Iterated density
double nu_n(double beta, int n, double x)
{
	if (0 == n) return nu(x);

	return rho_n(beta, n, x) - dee_n(beta, n, x);
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

// #define UNIT_TEST
#ifdef UNIT_TEST
	printf("#\n# beta=%g\n#\n", beta);

	int imax = 20;
	for (int i=0; i< imax; i++)
	{
		double x = (((double) i) + 0.5) / ((double) imax);
		// double y = nu_n(beta, n, x);

#if 0
		// N=1 test cases pass
		double ef10 = nu (x/beta)/beta;
		double f10 = f_n_k(beta, 1, 0, x);
		printf("%d	%g f10=%g ef10=%g\n", i, x, f10, ef10);

		double eg11 = nu ((x+1.0)/beta)/beta;
		double g11 = g_n_k(beta, 1, 1, x);
		printf("%d	%g g11=%g eg11=%g\n", i, x, g11, eg11);

		// N=2 efff test cases pass
		double ef21 = g_n_k(beta, 1, 1, x/beta)/beta;
		double f21 = f_n_k(beta, 2, 1, x);
		printf("%d	%g f21=%g ef21=%g\n", i, x, f21, ef21);

		double ef20 = nu(x/beta)/beta + f_n_k(beta, 1, 0, x/beta)/beta;
		double f20 = f_n_k(beta, 2, 0, x);
		printf("%d	%g f20=%g ef20=%g\n", i, x, f20, ef20);

		double eg21 = nu((x+1.0)/beta)/beta + f_n_k(beta, 1, 0, (x+1.0)/beta)/beta;
		double g21 = g_n_k(beta, 2, 1, x);
		printf("%d	%g g21=%g eg21=%g\n", i, x, g21, eg21);

		double arg = (x+b_k(beta, 1)) / beta;
		double eg22 = g_n_k(beta, 1, 1, arg) /beta;
		double g22 = g_n_k(beta, 2, 2, x);
		printf("%d	%g g22=%g eg22=%g\n", i, x, g22, eg22);

		for (int k=0; k < 6; k++)
		{
			for (int n=k+1; n<10; n++)
			{
				double efn1 = f_n_k(beta, n, k, x);
				double fn1 = fsum_n_k(beta, n, k, x);
				printf("%d	%f   %d %d  fn1=%f efn1=%f  diff=%g\n",
					i, x, n, k, fn1, efn1, fn1-efn1);
			}
			printf("---------\n");
		}
#endif

		for (int n=0; n<10; n++)
		{
			double egn1 = g_n_1(beta, n, x);
			double gn1 = gsum_n_1(beta, n, x);
			printf("%d	%g   %d  gn1=%g egn1=%g  diff=%g\n",
				i, x, n, gn1, egn1, gn1-egn1);
		}
		printf("---------\n");
	}
#endif

#define PRINT_DEN
#ifdef PRINT_DEN
#define NIT 3
	double sum[NIT];
	for (int j=0; j<NIT; j++) sum[j] = 0.0;

	double lambda = 1.0;
	// double lambda = 1.0 / beta;
	// double lambda = 1.0 / (beta*beta);
	// double lambda = 1.0 / (beta*beta*beta);
	double lamn = pow(lambda, n);

	int imax = 301;
	double delta = 1.0 / ((double) imax);
	for (int i=0; i< imax; i++)
	{
		double x = (((double) i) + 0.5) / ((double) imax);
		double y = gp_invar(beta, x);
		// double y = gp_n2(beta, x);
		// double y = gp_n3(beta, x);
		// double y = gp_quad_n1(beta, x);
		printf("%d	%g	%g", i, x, y);

		double lscale = lamn;
		for (int j=0; j<NIT; j++)
		{
			double y = nu_n(beta, n+j, x);
			y /= lscale;
			lscale *= lambda;

			sum[j] += y * delta;
			printf("	%g", y);
		}
		printf("\n");
		fflush(stdout);
	}

	printf("#\n# ");
	for (int j=0; j<NIT; j++)
		printf(" %g", sum[j]);
	printf("\n#\n");
#endif
}
