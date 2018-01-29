/* 
 * berg.c
 *
 * Crazy bergman polynomials
 *
 * January 2018
 */

#include <math.h>
#include <complex.h>

#define NOMAIN
#include "psi.c"
#include "psibig.c"

double complex bergman(double Kay, int n, double complex z)
{
	if (0 == n) return 1.0;
	double complex acc = 0;
	for (int j=0; j< n; j++)
	{
		acc += hess(Kay, j, n-1) * bergman(Kay, j, z);
	}
	acc = z * bergman(Kay, n-1, z) - acc;
	acc /= hess(Kay, n,n-1);
	return acc;
}

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

 
