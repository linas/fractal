/**
 *
 * Extrapolate eigenvalue and eigenvector finding
 * Call lapack from C code.
 *
 * Get Schur factorization,
 * G = OTO^ where O is orthogonal (O^ == O transpose)
 * T is triangular.
 *
 * Use LAPACK xGEEV routine
 * specifically, complex double precision:
 * ZGEEV
 *
 * Apply to find eigenvalues of the Fourier-of-Minkowski-question
 * operator
 *
 * Linas Vepstas September 2004
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "lapack.h"
#include "gral-simple.h"


/* dry run -- get the working dimension */
int 
cgetworkdim (int dim, complex double *matrix, 
        complex double *eigenvalues, 
		  complex double *left_eigen, complex double *right_eigen,
		  complex double *workspace, double *w2)
{

	char jobvl = 'V';
	char jobvr = 'V';
	int workdim;
	int info;

	workdim = -1;

	zgeev_ (&jobvl, &jobvr, &dim, matrix, &dim, 
	       eigenvalues, 
			 left_eigen, &dim, right_eigen, &dim, 
			 workspace, &workdim, w2, &info);

	return (int) (workspace[0]+0.5);
}

int 
cgeteigen (int dim, complex double *matrix, 
        complex double *eigenvalues, 
		  complex double *left_eigen, complex double *right_eigen,
		  int workdim, complex double *workspace, double *w2)
{
	char jobvl = 'V';
	char jobvr = 'V';
	int info;

	zgeev_ (&jobvl, &jobvr, &dim, matrix, &dim, 
	       eigenvalues, 
			 left_eigen, &dim, right_eigen, &dim, 
			 workspace, &workdim, w2, &info);

}


main (int argc, char * argv[]) 
{
	complex double *mat;
	complex double *regmat;
	complex double *eval;
	complex double *lev;
	complex double *rev;
	complex double *work;
	double *w2;
	int dim;
	int workdim;
	int i,j, k;
	
	dim = 28;
	set_npts(923123);

	if (argc < 2)
	{
		fprintf (stderr, "Usage: %s <dim>\n", argv[0]);
		exit (-1);
	}
	dim = atoi (argv[1]);

	printf ("#\n#\n");
	printf ("# Eigenvectors of fourier of question operator\n");
	printf ("# or maybe something else, depending on the version of this code\n");
	printf ("# Numerically solved to rank=%d\n", dim);
	printf ("#\n#\n");

	mat = (complex double *) malloc (dim*dim*sizeof (complex double));
	regmat = (complex double *) malloc (dim*dim*sizeof (complex double));
	eval = (complex double *) malloc (dim*sizeof (complex double));
	lev = (complex double *) malloc (dim*dim*sizeof (complex double));
	rev = (complex double *) malloc (dim*dim*sizeof (complex double));
	workdim = 4*dim*dim;
	work = (complex double *) malloc (workdim*sizeof (complex double));
	w2 = (double *) malloc (2*dim*sizeof (double));

	/* Insert values for the operator */
	for (i=0; i<dim; i++)
	{
		for (j=0; j<dim; j++)
		{
			/* Note transposed matrix'ing for FORTRAN */
			// mat[i+j*dim] = ache_mp(i,j);
			// mat[i+j*dim] = mtm_svd(i,j);

			long double re;
			long double im;
			int m = i-dim/2;
			int n = j - dim/2;
			make_elt (m,n, &re, &im);

			complex double z = re + I * im;
			z *= exp(-0.1*(m*m+n*n)/(dim*dim));

			mat[i+j*dim] = z;
			// printf ("mat(%d, %d) = %g\n", i,j,mat[i+j*dim]);
		}
		// printf("\n");
	}

	
	int wd = cgetworkdim (dim, mat, eval, lev, rev, work, w2);
	printf ("# recommended dim=%d actual dim=%d\n#\n", wd, workdim);
	// workdim = wd;
	
	work = (complex double *) realloc (work, workdim*sizeof (complex double));

	// Now, start regulating
	double t;
	for (t=32.0; t>1.0e-4; t /= sqrt(sqrt(sqrt(2))))
	{
		for (i=0; i<dim; i++)
		{
			for (j=0; j<dim; j++)
			{
				int m = i-dim/2;
				int n = j - dim/2;
				double reg = exp(-t*(m*m+n*n)/(dim*dim));
				regmat[i+j*dim] = reg * mat[i+j*dim];
			}
		}

		// examine the edges of the matrix, see how bad things are.
		double edgemax = 0.0;
		for (i=0; i<dim; i++)
		{
			if (edgemax < cabs(regmat[i])) edgemax = cabs(regmat[i]);
			if (edgemax < cabs(regmat[i+(dim-1)*dim])) edgemax = cabs(regmat[i+(dim-1)*dim]);
			if (edgemax < cabs(regmat[i*dim])) edgemax = cabs(regmat[i*dim]);
			if (edgemax < cabs(regmat[dim-1+i*dim])) edgemax = cabs(regmat[dim-1+i*dim]);
		}

		cgeteigen (dim, regmat, eval, lev, rev, workdim, work, w2);

		/* print the eigenvalues */
		printf("%g\t%g", t, edgemax);

		for (i=0; i<25; i++)
		{
			printf ("\t%g\t%g\t%g", cabs(eval[i]), creal(eval[i]), cimag(eval[i]));
		}
		printf ("\n");
	}
}
