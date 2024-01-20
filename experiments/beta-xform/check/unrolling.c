/*
 * unrolling.c
 * Verify to unrolled recursion relations for generalized
 * stretch-cut-stack map.
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

// Arbitrary function
double nu(double x)
{
	// return 1.0;
	// return x-0.5;

	// Bernoulli poly B_2
	// return x*x - x  + 1.0 / 6.0;

	// Bernoulli poly B_3
	return x*x*x - 1.5*x*x  + 0.5*x;

	// Bernoulli poly B_4
	return x*x*x*x - 2.0*x*x*x  + x*x - 1.0/30.0;
}

// Forward decl
double g_n_1(double beta, int n, double x);

// Return the g_n_k constant from the "generalized stretch-cut-stack"
// section of paper. This is computed using the recursing formula.
double g_n_k(double beta, int n, int k, double x)
{
	if (n < k) return 0.0;
	if (0 == k) return nu(x);

	if (1 == k) return g_n_1(beta, n, x);

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
// section of paper.
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

// Iterated density
double rho_n(double beta, int n, double x)
{
	double sum = 0.0;
	for (int k=0; k <=n; k++)
	{
		double tk = t_k(beta, k);
		if (tk > x)
			sum += g_n_k(beta, n, k, x);
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
			sum -= f_n_k(beta, n, k, x);
	}
	return sum;
}

// Iterated density
double nu_n(double beta, int n, double x)
{
	if (0 == n) return nu(x);

	return rho_n(beta, n, x) - dee_n(beta, n, x);
}

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		fprintf(stderr, "Usage: %s beta n\n", argv[0]);
		exit (1);
	}
	double beta = atof(argv[1]);
	int n = atoi(argv[2]);

#if UNIT_TEST
	printf("#\n# beta=%g\n#\n", beta);

	int imax = 20;
	for (int i=0; i< imax; i++)
	{
		double x = (((double) i) + 0.5) / ((double) imax);
		// double y = nu_n(beta, n, x);

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
	}
#endif

#define PRINT_DEN
#ifdef PRINT_DEN
#define NIT 3
	double sum[NIT];
	for (int j=0; j<NIT; j++) sum[j] = 0.0;

	// double lambda = 1.0 / beta;
	double lambda = 1.0 / (beta*beta);
	// double lambda = 1.0 / (beta*beta*beta);
	double lamn = pow(lambda, n);

	int imax = 301;
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
		printf(" %g", sum[j] - 1.0);
	printf("\n#\n");
#endif
}
