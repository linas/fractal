/*
 * psieign.c
 *
 * Eigenvalues of downshift from Hessenberg matrix functions.
 * Construct wave functions which put the downshift operator
 * into Hessenberg form.
 *
 * January 2018
 */
#define NOMAIN 1
#include "psi.c"

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

/* 
 * Eigenvalue and eigenvector problem.
 * Lapack, Eigenvalues only, use: HSEQR
 * Eigenvalues an eigenvectors: HSEQR TREVC
 *
 * GSL has a simpler API, so try that first.
 */
void eigen(double K, int dim)
{
	int mxi = MAXN-1;
mxi = dim;
	double* matrix = malloc(mxi*mxi*sizeof(double));

	// Do we need this, or the transpose ???
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
		printf ("%d eigenvalue = %g + %gi\n", i, GSL_REAL(eval_i), GSL_IMAG(eval_i));
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
	printf("#\n# K=%g\n#\n", K);

#if 0
#define NPTS 201
	double s = 0.0;
	for (int i=0; i< NPTS; i++)
	{
		double x = ((double) i + 0.5) / ((double) NPTS);
		double (*f)(double, double) = psi_1;
		double ef = f(x, K);
		double p0 = part(x, K, f, 0);
		double p1 = part(x, K, f, 1);
		double y = xfer(x, K, f);
		s += y / ((double) NPTS);
		printf("%d	%g	%g %g	%g	%g	%g\n", i, x, ef, p0, p1, y, s);
	}
#endif

	find_midpoints(K);
	verify_ortho();
	eigen(K, dim);
}
