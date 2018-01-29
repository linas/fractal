/* 
 * berg.C
 *
 * Crazy Bergman polynomials
 *
 * January 2018
 */

#include <math.h>
#include <complex.h>

#include "brat.h"

#define NOMAIN
#include "psi.c"
#include "psibig.c"

#define WTF_COMPLEX std::complex<double>

// Bergman polynomials for a Hessenberg operator, as per 
// Edward B. Saff and Nikos Stylianopoulos, Asymptotics for Hessenberg 
WTF_COMPLEX bergman(double Kay, int n, WTF_COMPLEX z)
{
	if (0 == n) return 1.0;
	WTF_COMPLEX acc = 0;
	for (int j=0; j< n; j++)
	{
		acc += hess(Kay, j, n-1) * bergman(Kay, j, z);
	}
	acc = z * bergman(Kay, n-1, z) - acc;
	acc /= hess(Kay, n,n-1);
	return acc;
}

double invariant_domain(double re_q, double im_q, int itermax, double Kay)
{
	static bool init = false;
	if (not init)
	{
		init = true;
		// find_midpoints(Kay);
		big_midpoints(Kay, 400, midpoints, MAXN);
		sequence_midpoints(Kay);
	}

	WTF_COMPLEX z = re_q + I * im_q;
	WTF_COMPLEX pz = bergman(Kay, itermax, z);
	//	double rv = abs(pz);
	rv = arg(pz);
	rv += M_PI;
	rv /= 2.0 * M_PI;
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
