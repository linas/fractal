
/*
 * sho.c:
 * compute hyperblic simple harmonic oscillator eigenvalues
 * Linas Vepstas November 2006
 */

#include <math.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>

void find_eigen(int nsz)
{
	gsl_matrix_complex *matrix;

	matrix = gsl_matrix_complex_alloc (nsz+1, nsz+1);
	gsl_matrix_complex_set_zero(matrix);

	int n;
	for (n=0; n<nsz; n++)
	{
		gsl_complex z;
		double y = sqrt ((2*n+2)*(2*n+1));
		// double y = sqrt ((2*n+2)*(2*n+3));
		// double y = 1;
		GSL_SET_COMPLEX(&z,0.0,y);
		gsl_matrix_complex_set(matrix, n,n+1, z);
		GSL_SET_COMPLEX(&z,0.0,-y);
		gsl_matrix_complex_set(matrix, n+1,n, z);
	}
	
	gsl_vector *eigenvals;
	eigenvals = gsl_vector_alloc (nsz+1);
	
	gsl_eigen_herm_workspace * work;
	work = gsl_eigen_herm_alloc (nsz+1);
	
	gsl_eigen_herm(matrix, eigenvals, work);
	gsl_sort_vector(eigenvals);

	for (n=0; n<nsz+1; n++)
	{
		double e = gsl_vector_get(eigenvals, n);
		double x = ((double) n)/ ((double) (nsz+1));
		x -= 0.5;
		e /= ((double) (nsz+1));
		printf ("%d	%g	%g\n", n, x, e);
	}
	
	gsl_eigen_herm_free(work);
	gsl_matrix_complex_free (matrix);
}

main(int argc, char * argv[])
{
	int n = atoi (argv[1]);
	find_eigen (n);
}
