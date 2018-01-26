/*
 * psiunitary.c
 *
 * Try to see if the downshift operator really is unitary (or not)
 * I mean given the eigenvalues...
 *
 * Linas Vepstas Nauary 2018
 */

#define NOMAIN
#include "psi.c"
#include "psibig.c"

/*
 * Compute L times L^T  (L-transpose)
 */
double right_prod (double K, int m, int n)
{
	double acc = 0.0;
	for (int i=0; i< 4000; i++)
	{
		acc += hess(K, m,i) * hess (K, n,i);
	}
}

void chk(double K)
{
	for (int m = 0; m<10; m++)
	{
		for (int n=0; n<10; n++)
		{
			double p = right_prod(K, m, n);
			printf("LL^T[%d, %d] = %g\n", m, n, p);
		}
	}
}

int main (int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s K\n", argv[0]);
		exit(1);
	}
	double K = atof(argv[1]);

	// find_midpoints(K);
	big_midpoints(K, 400, midpoints, MAXN);
	sequence_midpoints(K);
	chk(K);
}
