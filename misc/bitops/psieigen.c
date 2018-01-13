/*
 * psieigen.c
 *
 * Eigenvalues of downshift from Hessenberg matrix functions.
 * Construct wave functions which put the downshift operator
 * into Hessenberg form.
 *
 * January 2018
 */
#define NOMAIN 1
#include "psi.c"
#include "psibig.c"

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

/*
 * Eigenvalue and eigenvector problem.
 * GSL has a simpler API than Lapack, so try that first.
 * XXX Except its not stable...
 */
void eigen(double K, int dim)
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

	// Magic incantation
	gsl_matrix_view m = gsl_matrix_view_array (matrix, mxi, mxi);
	gsl_vector_complex *eval = gsl_vector_complex_alloc (mxi);
	gsl_matrix_complex *evec = gsl_matrix_complex_alloc (mxi, mxi);
	gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (mxi);
	gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);
	gsl_eigen_nonsymmv_free (w);
	gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);

	
	for (int i=0; i< mxi; i++)
	{
		gsl_complex eval_i = gsl_vector_complex_get (eval, i);

		double x = GSL_REAL(eval_i);
		double y = GSL_IMAG(eval_i);
		double mag = sqrt(x*x+y*y);
		printf ("%d	%g	%g	%g\n", i, x, y, mag);
#if 0
		printf ("eigenvector = \n");
		for (j=0; j< mxi; j++)
		{
			gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
			printf("%g + %gi\n", GSL_REAL(z), GSL_IMAG(z));
		}
#endif
	}

	gsl_vector_complex_free(eval);
	gsl_matrix_complex_free(evec);
	free(matrix);
}


int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s K dim\n", argv[0]);
		exit(1);
	}
	double K = atof(argv[1]);
	int dim = atoi(argv[2]);
	printf("#\n# GSL variant K=%g\n#\n", K);

	// find_midpoints(K);
	big_midpoints(K, 4000, midpoints, MAXN);
	sequence_midpoints(K);
	verify_ortho();
	eigen(K, dim);
}
