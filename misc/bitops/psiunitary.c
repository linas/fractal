/*
 * psiunitary.c
 *
 * Try to see if the downshift operator really is unitary (or not)
 * I mean given the eigenvalues...
 *
 * Linas Vepstas January 2018
 */

#define NOMAIN
#include "psi.c"
#include "psibig.c"
#include "psifp.c"

/* Return the matrix elements of the downshift operator,
 * but with the primary frobenius-perron eigenvector subtracted.
 * Then rescale by 2K, so that the resuulting matrix has eigenvalues
 * on the unit cicle, only. (And at least one at zero).
 */
double decay(double K, int m, int n, int dim)
{
	static double* fpvec = NULL;
	if (NULL == fpvec)
	{
		fpvec = malloc(dim * sizeof(double));
		get_fp_eigenvector(K, fpvec, dim);
	}

	return 2.0 * K * (hess(K, m, n) - fpvec[m]*fpvec[n]);
}

/*
 * Compute L times L^T  (L-transpose)
 */
double right_prod (double K, int m, int n, int dim)
{
	double acc = 0.0;
	for (int i=0; i< dim; i++)
	{
		acc += decay(K, m, i, dim) * decay (K, n, i, dim);
		// acc += decay(K, i, m, dim) * decay (K, i, n, dim);
	}
	return acc;
}

void chk(double K, int dim)
{
	for (int m = 0; m<10; m++)
	{
		for (int n=0; n<10; n++)
		{
			double p = right_prod(K, m, n, dim);
			if (fabs(p) < 1.0e-14) p = 0.0;
			else
				printf("LL^T[%d, %d] = %g\n", m, n, p);
		}
	}
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
	chk(K, dim);
}
