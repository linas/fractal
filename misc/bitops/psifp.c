/*
 * psifp.c
 *
 * Return the Ruelle-Frobenius-Perron eigenvector
 * i.e. the one with eigenvalue=1
 *
 * Linas Vepstas January 2018
 */

#ifndef NOMAIN
#define NOMAIN
#include "psi.c"
#include "psibig.c"
#undef NOMAIN
#endif

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
	double* matrix = (double*) malloc(mxi*mxi*sizeof(double));

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


	/* Get the eigenvector for the largest eigenvalue, which
	 * should be 1.0.  Graph it to make sure it matches what
	 * we already know it should be.
	 */
	gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);

	gsl_complex eval_0 = gsl_vector_complex_get (eval, 0);

	double x = GSL_REAL(eval_0);
	double y = GSL_IMAG(eval_0);
	double mag = sqrt(x*x+y*y);
	printf ("# eigenvalue = %20.18g + i %20.18g Magnitude=%g\n#\n", x, y, mag);

	if (1.0e-10 < fabs(1.0 - mag))
		printf ("# ERROR !!!!!!!!! NOT THE UNIT EIGENVALUE!\n#\n");

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

#ifndef NOMAIN
int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s K\n", argv[0]);
		exit(1);
	}
	double K = atof(argv[1]);
	printf("#\n# K=%g\n#\n", K);
	printf("#\n# eigenvector componentns\n#\n");

	// find_midpoints(K, MAXN);
	big_midpoints(K, 400, midpoints, MAXN);
	sequence_midpoints(K, MAXN);

#define NPTS 201
	double* fpvec = (double*) malloc(NPTS * sizeof(double));
	get_fp_eigenvector(K, fpvec, NPTS);

	for (int i=0; i< NPTS; i++)
	{
		if (fpvec[i] < 0.0 && fabs(fpvec[i])<1.0e15) fpvec[i] = 0.0;
		printf("%d	%g\n", i, fpvec[i]);
	}
}
#endif
