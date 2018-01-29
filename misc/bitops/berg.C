/*
 * berg.C
 *
 * Bergman polynomials that correspond to the beta-transform.
 *
 * The main result is that the holomorphic function that corresponds
 * to the invariant measure is zero inside the unit circle, and infinite
 * outside of it.
 *
 * January 2018
 */

#include <math.h>
#include <complex.h>

#include "brat.h"

#define NOMAIN
#include "psi.c"
#include "psibig.c"
#include "psifp.c"

#define WTF_COMPLEX std::complex<double>

// Bergman polynomials for a Hessenberg operator, as per
// Edward B. Saff and Nikos Stylianopoulos, Asymptotics for Hessenberg
// Recursive implmentation gets real slow for large n.
WTF_COMPLEX bergman_recursive(double Kay, int n, WTF_COMPLEX z)
{
	if (0 == n) return 1.0;
	WTF_COMPLEX acc = 0;
	for (int j=0; j< n; j++)
	{
		acc += hess(Kay, j, n-1) * bergman_recursive(Kay, j, z);
	}
	acc = z * bergman_recursive(Kay, n-1, z) - acc;
	acc /= hess(Kay, n,n-1);
	return acc;
}

// Non-recursive, much faster version, assumes values
// for k<n available in the cache array.
WTF_COMPLEX bergman_cache(WTF_COMPLEX cache[], double Kay, int n, WTF_COMPLEX z)
{
	if (0 == n) return 1.0;

	WTF_COMPLEX acc = 0;
	for (int j=0; j< n; j++)
	{
		acc += hess(Kay, j, n-1) * cache[j];
	}
	acc = z * cache[n-1] - acc;
	acc /= hess(Kay, n, n-1);
	cache[n] = acc;
	return acc;
}

// Wrapper around the cached version.
WTF_COMPLEX bergman(double Kay, int n, WTF_COMPLEX z)
{
	if (0 == n) return 1.0;
	WTF_COMPLEX cache[n];
	cache[0] = 1.0;
	for (int j=1; j<=n; j++)
	{
		bergman_cache(cache, Kay, j, z);
	}
	return cache[n];
}

// Ruelle-Frobenius-Perron Bergman eigenfunction
WTF_COMPLEX rfp_bergman(double* vec, double Kay, WTF_COMPLEX z, int smax)
{
	int n=smax;

	// Precompute the cache of values.
	WTF_COMPLEX cache[n];
	cache[0] = 1.0;
	for (int j=1; j<=n; j++)
	{
		bergman_cache(cache, Kay, j, z);
	}

	// Compute the vector.
	WTF_COMPLEX acc = 0.0;
	for (int j=0; j<=n; j++)
	{
		acc += vec[j] * cache[j];
	}

// printf(" its %g %g  so %g %g\n", real(z), imag(z), real(acc), imag(acc));
	return acc;
}

double invariant_domain(double re_q, double im_q, int itermax, double Kay)
{
	static double* fpvec = NULL;
	static bool init = false;
	if (not init)
	{
		init = true;
		// find_midpoints(Kay);
		big_midpoints(Kay, 400, midpoints, MAXN);
		sequence_midpoints(Kay);

		fpvec = (double*) malloc(itermax * sizeof(double));
		get_fp_eigenvector(Kay, fpvec, itermax);
	}

	WTF_COMPLEX z = re_q + I * im_q;
	// WTF_COMPLEX pz = bergman(Kay, itermax, z);
	WTF_COMPLEX pz = rfp_bergman(fpvec, Kay, z, itermax);

#define MAG
#ifdef MAG
	double rv = abs(pz);
#else
	double rv = arg(pz);
	rv += M_PI;
	rv /= 2.0 * M_PI;
#endif

double zabs = abs(z);
if(abs(zabs - 1.0/(2.0*Kay)) < 0.003) return 0.35;
if(abs(zabs - 1.0) < 0.003) return 0.5;

	// printf("duuude uh %g for %g %d\n", rv, Kay, itermax);
	return rv;
}

DECL_MAKE_HEIGHT(invariant_domain)

#if 0
int main (int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s K dim\n", argv[0]);
		exit(1);
	}
	double K = atof(argv[1]);
	int dim = atoi(argv[2]);

	// find_midpoints(K);
	big_midpoints(K, 400, midpoints, MAXN);
	sequence_midpoints(K);

	double lam = 1.0;
	for (int i=0; i<dim; i++)
	{
		double h = hess(K, i+1, i);
		double d = hess(K, i, i);
		double d1 = hess(K, i-1, i);
		double d2 = hess(K, i-2, i);
		lam /= h;
		printf("%d	%g	%g	%g	%g	%g\n", i, h, d, d1, d2, lam);
	}
}

#endif
