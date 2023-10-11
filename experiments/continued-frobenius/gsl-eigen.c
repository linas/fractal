
/* gsl-eigen.c
 *
 * bi-diagonalize the H matrix
 */


#include <stdio.h>
#include <gsl/gsl_linalg.h>

#include "ache.h"

double * 
h_trunc (int dim)
{
	static double * h = NULL;

	h = (double *) realloc (h,  dim*dim*sizeof (double));

	int i,j;
	for (i=0; i<dim; i++)
	{
		for (j=0; j<dim; j++)
		{
			h[i*dim+j] = ache_mp(i,j);
		}
	}
	return h;
}

int
main (void)
{
	int dim;
	for (dim=2; dim<25; dim++)
	{
		double *ht_data = h_trunc (dim);
		gsl_matrix_view ht = gsl_matrix_view_array (ht_data, dim, dim);
		gsl_vector *tau_U = gsl_vector_alloc (dim);
		gsl_vector *tau_V = gsl_vector_alloc (dim-1);
		gsl_linalg_bidiag_decomp (&ht.matrix, tau_U, tau_V);

		gsl_vector *diag = gsl_vector_alloc (dim);
		gsl_vector *superdiag = gsl_vector_alloc (dim-1);
		gsl_linalg_bidiag_unpack_B (&ht.matrix, diag, superdiag);

		printf ("duude dim=%d diag = \n", dim);
		gsl_vector_fprintf (stdout, diag, "%g");
		printf ("\nduude super diag = \n");
		gsl_vector_fprintf (stdout, superdiag, "%g");

		printf ("-------------\n\n");
	}

	return 0;
}

