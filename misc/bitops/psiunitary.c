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


#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

/*
 * Get the Frobenius-Perron eigenvector (the one with eigenvalue 1.0).
 * Write it into the array "vec" which should be of at least dimension
 * dim.
 *
 * This uses the GSL rouines to compute it. The correctness of the beast
 * was previously verified.
 */
void get_fp_eigenvector(double K, double* vec, int dim)
{
	int mxi = MAXN-1;
mxi = dim;
	double* matrix = malloc(mxi*mxi*sizeof(double));

	// Do we need this, or the transpose ??? Does it matter?
	int idx = 0;
	for (int i=0; i< mxi; i++)
	{
		for (int j=0; j< mxi; j++)
		{
			double h = hess(K, i, j);
			matrix[idx] = h;
			idx++;
		}
	}

	// Magic incantation to diagonalize the matrix.
	gsl_matrix_view m = gsl_matrix_view_array (matrix, mxi, mxi);
	gsl_vector_complex *eval = gsl_vector_complex_alloc (mxi);
	gsl_matrix_complex *evec = gsl_matrix_complex_alloc (mxi, mxi);
	gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (mxi);
	gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);
	gsl_eigen_nonsymmv_free (w);


	/* Get the eigenvector for the largest eigennvalue, which
	 * should be 1.0.  Graph it to make sure it matches what
	 * we already know it should be.
	 */
	gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);

	gsl_complex eval_0 = gsl_vector_complex_get (eval, 0);

	double x = GSL_REAL(eval_0);
	double y = GSL_IMAG(eval_0);
	double mag = sqrt(x*x+y*y);
	printf ("# eigenvalue = %20.18g + i %20.18g Magnitude=%g\n#\n", x, y, mag);

	/* Get the eigenvector */
	gsl_vector_complex_view evec_0 = gsl_matrix_complex_column (evec, 0);

	gsl_complex z = gsl_vector_complex_get(&evec_0.vector, 0);
	double sgn = GSL_REAL(z);
	if (0 < sgn) sgn = 1.0; else sgn = -1.0;

	mag = 0.0;
	for (int j=0; j< mxi; j++)
	{
		gsl_complex z = gsl_vector_complex_get(&evec_0.vector, j);
		// printf("# %d	%g	%g\n", j, GSL_REAL(z), GSL_IMAG(z));
		double re = sgn * GSL_REAL(z);
		mag += re * re;
		vec[j] = re;
	}
	mag = sqrt(mag);
	printf("# Vector length = %20.18g\n#\n", mag);
}

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
