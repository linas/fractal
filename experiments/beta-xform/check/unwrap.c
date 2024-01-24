/*
 * unwrap.c
 *
 * Verify to unwrapped recursion relations for generalized
 * stretch-cut-stack map. These are the ones after noticing that
 * the unroll could be made prettier. Works, but is obsoleted by
 * code in unstack.c which is even more unrolled.
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
	if (x < 0.0 || 1.0 < x) exit(1);

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
double h_n_1(double beta, int n, double x);

// Return the h_n_k constant from the "generalized stretch-cut-stack"
// section of paper. This is computed using the recursing formula.
double h_n_k(double beta, int n, int k, double x)
{
	if (n < k) return 0.0;
	if (n < 0) { fprintf(stderr, "Error: negative n\n"); return 0.0; }
	if (0 == n && 0 == k) return nu(x);
	if (0 == k) return 0.0;

	if (1 == k) return h_n_1(beta, n, x);

	// Recurse
	double bkm1 = b_k(beta, k-1);
	double arg = (x + bkm1) / beta;
	double guh = h_n_k(beta, n-1, k-1, arg);
	return guh / beta;
}

// Return the e_n_k constant from the "generalized stretch-cut-stack"
// section of paper. This is computed using the recursing formula.
double e_n_k(double beta, int n, int k, double x)
{
	if (n <= k) return 0.0;
	if (k < 0) { fprintf(stderr, "Error: can't negative k\n"); return 0.0; }
	double arg = x / beta;
	double sum = h_n_k(beta, n-1, k, arg);
	if (k < n) sum += e_n_k(beta, n-1, k, arg);
	return sum / beta;
}

// Return the h_n_1 constant from the "generalized stretch-cut-stack"
// section of paper. This is computed with the recursive formula.
double h_n_1(double beta, int n, double x)
{
	if (n < 1) { fprintf(stderr, "Error: badness\n"); return 0.0; }

	double arg = (x + 1.0) / beta;
	if (1 == n) return nu(arg) / beta;

	double sum = 0.0;
	for (int k=0; k<n-1; k++)
	{
		if (b_k(beta, k))
			sum += e_n_k(beta, n-1, k, arg);
	}
	return sum / beta;
}

// ==============================================================
// Series summations

double hsum_n_k(double beta, int n, int k, double x);

// Return the h_n_1 constant from the "generalized stretch-cut-stack"
// section of paper. This is computed with the series sum.
// Does pointless inline; which does nothing but increase complexity.
double hsumin_n_1(double beta, int n, double x)
{
	double arg = (x + 1.0) / beta;
	if (1 == n) return nu(arg) / beta;

	double sum = 0.0;
	for (int k=0; k<= n-2; k++)
	{
		if (0 == b_k(beta, k)) continue;

		double bej = 1.0 / (beta * beta);
		double arg = (x + 1.0) / (beta * beta);

		for (int j=2; j <= n-k; j++)
		{
			sum += bej * hsum_n_k(beta, n-j, k, arg);
			bej /= beta;
			arg /= beta;
		}
	}

	return sum;
}

// Foward decl
double hsum_n_1(double beta, int n, double x);

// Return the h_n_k constant from the "generalized stretch-cut-stack"
// section of paper. This is computed using the summation formula.
double hsum_n_k(double beta, int n, int k, double x)
{
	if (n < k) return 0.0;
	if (n < 0) { fprintf(stderr, "Error: negative n\n"); return 0.0; }
	if (0 == n && 0 == k) return nu(x);
	if (0 == k) return 0.0;

	if (1 == k) return hsum_n_1(beta, n, x);

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
	return hsum_n_1(beta, n-k+1, arg) * bek;
}

// Return the e_n_k constant from the "generalized stretch-cut-stack"
// section of paper. This is computed using the summation formula.
double esum_n_k(double beta, int n, int k, double x)
{
	if (n <= k) return 0.0;
	if (k < 0) fprintf(stderr, "Error can't negative k\n");
	double sum = 0.0;
	double bej = 1.0;
	for (int j=1; j <= n-k; j++)
	{
		bej /= beta;
		sum += bej * hsum_n_k(beta, n-j, k, x*bej);
	}
	return sum;
}

// Return the h_n_1 constant from the "generalized stretch-cut-stack"
// section of paper. This is computed with the recursive formula.
double hsum_n_1(double beta, int n, double x)
{
	if (n < 1) { fprintf(stderr, "Error: badness\n"); return 0.0; }

	double arg = (x + 1.0) / beta;
	if (1 == n) return nu(arg) / beta;

	double sum = 0.0;
	for (int k=0; k<n-1; k++)
	{
		if (b_k(beta, k))
			sum += esum_n_k(beta, n-1, k, arg);
	}
	return sum / beta;
}

// ==============================================================
// Density expressions

// Iterated density
double sigma_n(double beta, int n, double x)
{
	double sum = 0.0;
	for (int k=0; k <=n; k++)
	{
		double tk = t_k(beta, k);
		if (tk > x)
			sum += hsum_n_k(beta, n, k, x);
	}
	return sum;
}

// Iterated density
double cee_n(double beta, int n, double x)
{
	double sum = 0.0;
	for (int k=0; k <n; k++)
	{
		if (b_k(beta, k))
			sum += esum_n_k(beta, n, k, x);
	}
	return sum;
}

// Iterated density
double nu_n(double beta, int n, double x)
{
	if (0 == n) return nu(x);
	return sigma_n(beta, n, x) + cee_n(beta, n, x);
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

	int imax = 314;
	double delta = 1.0 / ((double) imax);
	for (int i=0; i< imax; i++)
	{
		double x = (((double) i) + 0.5) / ((double) imax);

		double y = gp_invar(beta, x);
		// double y = gp_n1(beta, x);
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
