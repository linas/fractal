
/* gsl-diag.c
 *
 * Attempt to diagonalize the continued-fraction map 
 * using the gnu gsl library eigenvector routines.
 * 
 * The effing thing appears to yeild divergent 
 * eigenvalues, so ugh, this numerical approach
 * seems not to work (again). Arghhh this is frustrating!
 *
 * Linas Dec 2003
 */

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include "zetafn.h"

#define MS 100

// set up the G-matrix, regulated with value tee
void fill_matrix (double *data, int dim, double tee)
{
	int i,j;

	for (i=0; i<dim; i++)
	{
		for (j=0; j<dim; j++)
		{
			double sign = 1.0;
			if ((i+j)%2 ==1) sign = -1.0;
			double z = 1.0 + zetam1(i+j+2);
			z *= factorial (i+j+1);
			z /= factorial (i) * factorial (j);
			z *= sign;
			z *= sqrt ((double) (i+1)*(j+1));
			double reg = i+j;
			reg = exp (-tee*tee*reg*reg);
			data [i*dim+j] = z *reg;
		}
	}
}

void solve (double * data, int dim)
{
	gsl_matrix_view m = gsl_matrix_view_array (data, dim, dim);

	gsl_vector *eigval = gsl_vector_alloc (dim);
	gsl_matrix *eigvec = gsl_matrix_alloc (dim, dim);

	gsl_eigen_symmv_workspace * gsw;
	gsw = gsl_eigen_symmv_alloc (dim);

	gsl_eigen_symmv (&m.matrix, eigval, eigvec, gsw);
   gsl_eigen_symmv_free (gsw);
	
	// we don't want to sort the eignevalues !!
	// gsl_eigen_symmv_sort (eigval, eigvec, 
   //                      GSL_EIGEN_SORT_ABS_DESC);
    
	int i;

	for (i = 0; i < dim; i++)
	{
		double eval_i = gsl_vector_get (eigval, i);
		gsl_vector_view evec_i = gsl_matrix_column (eigvec, i);

		printf ("eigenvalue %d = %g\n", i, eval_i);
		// printf ("eigenvector = \n");
		// gsl_vector_fprintf (stdout, &evec_i.vector, "%g");
	}
}

void doit (int dim)
{
	double *data;
	data = malloc (MS*MS*sizeof(double));

	double tee;
	for (tee=0.8; tee>1.0e-5; tee *= 0.8)
	{
		printf ("================== dim=%d tee=%f below\n", dim, tee);
		fill_matrix (data, dim, tee);
		solve (data, dim);
	}
}

main(int argc, char *argv[])
{
	int dim = 10;
	if (argc ==2)
	{
		dim = atoi (argv[1]);
	}
	doit(dim);
}

