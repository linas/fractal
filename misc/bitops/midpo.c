/*
 * midpo.c
 *
 * Stupid wild guess
 *
 * January 2018
 */

#include <math.h>
#include <complex.h>

#define NOMAIN
#include "psi.c"
#include "psibig.c"

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

// Stupid wild-ass guess -- stick midpoints into a shift matrix
// what happens? anything at all? doe this have any meaning?
void swag(double K, int mxi)
{
	double* matrix = (double*) malloc(mxi*mxi*sizeof(double));
	int idx = 0;
	for (int i=0; i< mxi; i++)
	{
		for (int j=0; j< mxi; j++)
		{
			double h = midpoints[i+j];
			// double h = midpoints[mid_sequence[i+j]];
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

	gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);

	for (int i=0; i< mxi; i++)
	{
		gsl_complex lamb_i = gsl_vector_complex_get (eval, i);
		double x = GSL_REAL(lamb_i);
		double y = GSL_IMAG(lamb_i);
		double mag = sqrt(x*x+y*y);
		printf("%d	%g	%g	%g\n", i, mag, x, y);
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

	swag(K,dim);
}
