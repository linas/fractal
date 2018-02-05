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
double bergman_oper(double Kay, int n, int j);
double bergman_oper_raw(double Kay, int n, int j)
{
	if (0 == n && 0 == j) return 1.0;
	if (n < j) return 0.0; // lower triangular
	if (j < 0) return 0.0;

	double acc = bergman_oper(Kay, n-1, j-1);
	for (int k=0; k<n; k++)
	{
		acc -= hess(Kay, k, n-1) * bergman_oper(Kay, k, j);
	}
	acc /= hess(Kay, n, n-1);
	
	return acc;
}

double bergman_oper(double Kay, int n, int j)
{
	static bool init = false;
	static double bcache[MAXN][MAXN];
	if (!init)
	{
		init = true;
		for (int p=0; p<MAXN; p++)
		{
			for (int q=0; q<MAXN; q++)
			{
				bcache[p][q] = -1e40;
			}
		}
	}
	double val = bcache[n][j];
	if (-1e30 < val) return val;

	val = bergman_oper_raw(Kay, n, j);
	bcache[n][j] = val;
	return val;
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
	// hess is zero if n+1 < k
	for (int k=m; k<=n+1; k++)
	{
		acc += hess(Kay, k, n) * bergman_oper(Kay, k, m);
	}
	return acc;
}

/**
 * Unit test the shift product. That is, H^T P = PK should hold.
 * This test is currently passing.
 */
void verify_shift(double Kay, int maxn)
{
	for (int m=0; m< maxn; m++)
	{
		double sum = 0.0;
		double absum = 0.0;
		for (int n=0; n<=maxn; n++)
		{
			double htp = hess_trans_berg(Kay, m, n);
			double poly = bergman_oper(Kay, m, n-1);
			sum += htp;
			absum += fabs(htp);
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
			printf("Error: col %d sum=%g should be zero fabs=%g\n\n", m, sum, absum);
		}
	}
}

/**
 * Return the right inverse of the bergman polynomial operator.
 * That is, return R such that PR = I the identity matrix,
 * Note that both R and P are lower-triangular.
 */
double rinverse(double Kay, int m, int n);
double rinverse_raw(double Kay, int m, int n)
{
	if (m<0 || n<0) return 0.0;
	if (m < n) return 0.0; // lower triangular
	if (m==n) return 1.0 / bergman_oper(Kay, n, n);

	double acc = 0.0;

	for (int k=n; k<m; k++)
	{
		acc -= bergman_oper(Kay, m, k) * rinverse(Kay, k, n);
	}
	acc /= bergman_oper(Kay, m, m);
	return acc;
}

double rinverse(double Kay, int n, int j)
{
	static bool init = false;
	static double bcache[MAXN][MAXN];
	if (!init)
	{
		init = true;
		for (int p=0; p<MAXN; p++)
		{
			for (int q=0; q<MAXN; q++)
			{
				bcache[p][q] = -1e40;
			}
		}
	}
	double val = bcache[n][j];
	if (-1e30 < val) return val;

	val = rinverse_raw(Kay, n, j);
	bcache[n][j] = val;
	return val;
}

/**
 * Unit test the inverse matrix.  When multiplied, it should
 * give the identity matrix.  Both P P^-1 = I = P^-1 P are tested.
 * Currently, this unit test is passing.
 */
void verify_inverse(double Kay, int nmax)
{
	for (int n=0; n<nmax; n++)
	{
		for (int m=0; m<nmax; m++)
		{
			double rprod = 0.0;
			double lprod = 0.0;
			for (int k=0; k<=nmax; k++)
			{
				rprod += bergman_oper(Kay, n, k) * rinverse(Kay, k, m);
				lprod += rinverse(Kay, n, k) * bergman_oper(Kay, k, m);
			}

			// Verify that we got the identity matrix
			if (n != m && 1.0e-12 < fabs(rprod))
			{
				printf("Error: right prod expecting zero at %d %d got %g\n", n, m, rprod);
			}
			if (n == m && 1.0e-12 < fabs(1.0 - rprod))
			{
				printf("Error: right prod expecting one at %d %d got %g\n", n, m, rprod);
			}

			// And the other one.
			if (n != m && 1.0e-12 < fabs(lprod))
			{
				printf("Error: left prod expecting zero at %d %d got %g\n", n, m, lprod);
			}
			if (n == m && 1.0e-12 < fabs(1.0 - lprod))
			{
				printf("Error: left prod expecting one at %d %d got %g\n", n, m, lprod);
			}
		}
	}
}

/**
 * Return the product bergman^-1 times hess-transpose times bergman
 * viz P^-1 H^T P .. this should just be the right-shift, nothing more.
 */
double binv_hess_trans_berg(double Kay, int n, int m)
{
	double acc = 0;

	// rinv is zero if n < k
	// H^TP is zero if m+1 < k
	for (int k=0; k<=n; k++)
	{
		acc += rinverse(Kay, n, k) * hess_trans_berg(Kay, k, m);
	}
	return acc;
}

/**
 * Verify that the above really does give the shift.
 * Yes, it really does.  This unit test passes.
 * So everything is sound.
 */
void verify_koopman(double Kay, int nmax)
{
	for (int n=0; n<nmax; n++)
	{
		for (int m=0; m<nmax; m++)
		{
			double koop = binv_hess_trans_berg(Kay, n, m);

			// Verify that we got the right-shift matrix
			if (n+1 != m && 1.0e-12 < fabs(koop))
			{
				printf("Error: koop expecting zero at %d %d got %g\n", n, m, koop);
			}
			if (n+1 == m && 1.0e-12 < fabs(1.0 - koop))
			{
				printf("Error: koop expecting one at %d %d got %g\n", n, m, koop);
			}
		}
	}
}

/**
 * If the domain of orthogonality was the unit disk,
 * (or is it was the disk r< 1/2K) then the below would
 * return the identity i.e. delta_mn . But it doesn't so
 * neither of these are the domain.  So wtf? What's the domain?
 */
double ortho(double Kay, int m, int n)
{
	double acc = 0.0;
	int nm = n<m ? m : n;
	// int nm = n<m ? n : m;
	for (int k=0; k<nm; k++)
	{
		double prod = bergman_oper(Kay, m, k) * bergman_oper(Kay, n, k);
		// double prod = rinverse(Kay, m, k) * rinverse(Kay, n, k);
		prod /= ((double) (k+1));
		// prod *= pow(2.0*Kay, -2*(k+1));
		acc += prod;
	}
	acc *= M_PI;
	return acc;
}

/**
 * scratch
 */
double scratch(double Kay, int m, int n)
{
	double acc = 0.0;
	int nm = n<m ? m : n;
	for (int k=0; k<nm; k++)
	{
		double prod = hess(Kay, k, m) * hess(Kay, k, n);
		acc += prod;
	}
	return acc;
}

double herm(double Kay, int m, int n)
{
	double acc = 0.0;
	int nm = n<m ? m : n;
	for (int k=0; k<nm; k++)
	{
		double prod = rinverse(Kay, m, k) * rinverse(Kay, n, k);
		acc += prod;
	}
	return acc;
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

K = 0.5 + K*0.001 * 0.32;
	printf("#\n# K=%g\n#\n", K);

	// find_midpoints(K);
	big_midpoints(K, 400, midpoints, MAXN);
	sequence_midpoints(K, MAXN);

// #define UNIT_TEST_POLY
#ifdef UNIT_TEST_POLY
	// The polynomial unit test is currently passing.
	cross_check_poly(K, maxn);

	// The shift unit test is currently passing.
	verify_shift(K, maxn);

	// The inverse test is currently passing.
	verify_inverse(K, maxn);

	// The koopman test is currently passing.
	verify_koopman(K, maxn);
#endif


#if PRINT_MATRIX
	for (int n=0; n<maxn; n++)
	{
		for (int m=0; m<maxn; m++)
		{
			// double poly = bergman_oper(K, n, m);
			// double inv = rinverse(K, n, m);
			// double ort = ortho(K, n, m);
			double sym = herm(K, n, m);
			// double koop = binv_hess_trans_berg(K, n, m);
			double try = sym;
			printf("[%d	%d] =	%g\n", n, m, try);
		}
		printf(" ---------------\n");
	}
#endif

#ifdef COLUMN_RATIOS
	// Well, its symmetric, so it doesn't matter: rows or cols.
	int n = 1;
	double prev = herm(K, n, 0);
	for (int m=1; m<maxn; m++)
	{
		double sym = herm(K, n, m);
		double rat = sym/prev;
		printf("%d	%d %g %g\n", n, m, sym, rat);
		fflush(stdout);
		prev = sym;
	}
#endif

#define HI_ACC_RATIOS
#ifdef HI_ACC_RATIOS
	int n = maxn;
	int m = maxn;
	// m = maxn*0.723;
	// m = maxn*0.411;
	// m = maxn*1.711;
	double prev = herm(K, n, m-1);
	double sym = herm(K, n, m);
	double rat = sym/prev;
	printf("%g	%d	%d %g %g\n", K, n, m, sym, rat);
#endif

#if DIAGONAL_ELEMENTS
	double ob = 1.0 / (2.0*K);
	double obn = ob;
	double prev = 1.0;
	for (int n=1; n<maxn; n++)
	{
		double poly = bergman_oper(K, n, n);
		double ren = poly * obn;
		double rat = poly / prev; 
		double inv = 1.0 / rat;
		double sub = hess(K, n, n-1);
		printf("%d	%g	%g	%g	%g	%g\n", n-1, poly, ren, rat, inv, sub);
		obn *= ob;
		prev = poly;
	}
#endif

#if BIZARRO
	// Don't know if this means anything, but for K=0.8
	// it appears that the product below H^T H has only one
	// non-zero entry, at [6,6] having value=0.625=1/1.6
	// All other entries in the 6th row are zero.
	// Other rows are not like this.
	for (int n=0; n<maxn; n++)
	{
		double acc = 0.0;
		int m = 6;
		for (int k=0; k<maxn; k++)
		{
			acc += hess(K, k, m) * hess(K, k, n);
		}
		printf("%d	%g\n", n, acc);
	}
#endif
}
