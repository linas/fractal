/*
 * unwrap.c
 *
 * Verify to unwrapped recursion relations for generalized
 * stretch-cut-stack map. The are the ones after noticing that
 * the unroll could be made prettier.
 *
 * January 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Invariant measure for beta shift (not beta transform)
double invar(double beta, double x)
{
	double midpnt = 0.5*beta;
	double obn = 1.0;
	double sum = 0.0;
	double norm = 0.0;
	for (int i=0; i<1000; i++)
	{
		if (x < midpnt) sum += obn;
		norm += midpnt*obn;

		if (0.5 < midpnt) midpnt -= 0.5;
		midpnt *= beta;
		obn /= beta;
		if (obn < 1e-15) break;
	}
	return sum / norm;
}

double gp_invar(double beta, double x)
{
	return 0.5*beta*invar(beta, 0.5*beta*x);
}

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
// Density expressions

// Iterated density
double sigma_n(double beta, int n, double x)
{
	double sum = 0.0;
	for (int k=0; k <=n; k++)
	{
		double tk = t_k(beta, k);
		if (tk > x)
			sum += h_n_k(beta, n, k, x);
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
			sum += e_n_k(beta, n, k, x);
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
	double lamn = pow(lambda, n);

	int imax = 314;
	double delta = 1.0 / ((double) imax);
	for (int i=0; i< imax; i++)
	{
		double x = (((double) i) + 0.5) / ((double) imax);
		double y = gp_invar(beta, x);
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
