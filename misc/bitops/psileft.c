/*
 * psileft.c
 * Compute the Bergman (left-hand) polynomial coefficients.
 * The bergman.C file does similar, but only uses Z.
 *
 * February 2018
 */
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define NOMAIN
#include "psi.c"
#include "psibig.c"


/** Evaluate the entries in the Bergman vector at some fixed z */
void bergman_vect(double Kay, complex z, complex* poly, int maxn)
{
	poly[0] = 1.0;
	for (int n=1; n<maxn; n++)
	{
		complex acc = z * poly[n-1];
		for (int k=0; k<n; k++)
		{
			acc -= hess(Kay, k, n-1) * poly[k];
		}
		acc /= hess(Kay, n, n-1);
		poly[n] = acc;
	}
}

/**
 * Just like above, but without the z values.
 * Returns the matrix entry p_{n,j} for the bergman polynomial matrix.
 * Computed recusrively, in closed form.
 */
double bergman_oper(double Kay, int n, int j)
{
	if (0 == n && 0 == j) return 1.0;
	if (n < j) return 0.0;
	if (j < 0) return 0.0;

	double acc = bergman_oper(Kay, n-1, j-1);
	for (int k=0; k<n; k++)
	{
		acc -= hess(Kay, k, n-1) * bergman_oper(Kay, k, j);
	}
	acc /= hess(Kay, n, n-1);
	
	return acc;
}

/** Evaluate polynomial with coefficeints given in poly,
 * up to degree */
complex eval_poly(complex z, complex* poly, int degree)
{
	complex zn = 1;
	complex acc = 0.0;
	for (int n=0; n<=degree; n++)
	{
		acc += zn * poly[n];
		zn *= z;
	}
	return acc;
}

/** Validate the two polynomial formulas one against the other.
 * They are bothe the same formula, we are only checking for
 * programming bugs.
 */
void cross_check_code(double Kay, int maxn, complex z)
{
	complex poly[maxn];
	bergman_vect(Kay, z, poly, maxn);

	for (int n=0; n<maxn; n++)
	{
		complex coef[n+1];
		for (int k=0; k<=n; k++)
		{
			coef[k] = bergman_oper(Kay, n, k);
		}
		complex pn = eval_poly(z, coef, n+1);

		complex diff = poly[n] - pn;
		double adiff = cabs(diff);
		if (1.0e-12 < adiff)
		{
			printf("Error at n=%d, z=%g + I %g diff=%g\n",
			        n, creal(z), cimag(z), adiff);
			printf("pzn = %g	%g\n", creal(pn), cimag(pn));
			printf("berg = %g	%g\n", creal(poly[n]), cimag(poly[n]));
		}
	}
}

/**
 * Sample a range of z's to make sure the computations
 * are all OK. As of right now, this test is passing.
 */
void cross_check_poly(double Kay, int maxn)
{
	complex z;
	int ixmax = 20;
	for (int ix = 0; ix< ixmax; ix++)
	{
		double x = ((double) ix + 0.5) / ((double) ixmax);
		for (int iy = 0; iy < ixmax; iy++)
		{
			double y = ((double) iy + 0.5) / ((double) ixmax);
			z = x + I * y;
			z = x * cexp (I * 2.0 * M_PI * y);
			cross_check_code(Kay, maxn, z);
		}
	}
}

/**
 * Return the product hess-transpose times bergman
 */
double hess_trans_berg(double Kay, int n, int m)
{
	double acc = 0;

	// berg is zero if k < m
	// hess iz zero if n+1 < k
	for (int k=m; k<=n+1; k++)
	{
		acc += hess(Kay, k, n) * bergman_oper(Kay, k, m);
	}
	return acc;
}

/**
 * Unit test the shift product
 */
void verify_shift(double Kay, int maxn)
{
	for (int m=0; m< maxn; m++)
	{
		double sum = 0.0;
		for (int n=0; n<=maxn; n++)
		{
			double htp = hess_trans_berg(Kay, m, n);
			double poly = bergman_oper(Kay, m, n-1);
			sum += htp;
			double diff = fabs(htp - poly);
			if (1.0e-12 < diff)
			{
				printf("Error: htp[%d, %d] = %g  berg=%g diff=%g\n",
				        m, n, htp, poly, diff);
			}
		}
		sum = fabs(sum);
		if (0 < m && 1.0e-12 < sum)
		{
			printf("Error: col %d sum=%g should be zero\n\n", m, sum);
		}
	}
}

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s K maxn\n", argv[0]);
		exit(1);
	}
	double K = atof(argv[1]);
	int maxn = atoi(argv[2]);
	printf("#\n# K=%g\n#\n", K);

	// find_midpoints(K);
	big_midpoints(K, 400, midpoints, MAXN);
	sequence_midpoints(K, MAXN);

// #define UNIT_TEST_POLY
#ifdef UNIT_TEST_POLY
	// The unit test is currently passing.
	cross_check_poly(K, maxn);
#endif

	verify_shift(K, maxn);
}
