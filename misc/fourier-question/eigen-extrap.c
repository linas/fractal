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
 * specifically, real double precision:
 * DGEEV
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
getworkdim (int dim, double *matrix,
        double *eigenvalues_re, double *eigenvalues_im, 
		  double *left_eigen, double *right_eigen,
		  double *workspace)
{

	char jobvl = 'V';
	char jobvr = 'V';
	int workdim;
	int info;

	workdim = -1;

	dgeev_ (&jobvl, &jobvr, &dim, matrix, &dim, 
	       eigenvalues_re, eigenvalues_im, 
			 left_eigen, &dim, right_eigen, &dim, 
			 workspace, &workdim, &info);

	return (int) (workspace[0]+0.5);
}

int 
geteigen (int dim, double *matrix, 
        double *eigenvalues_re, double *eigenvalues_im, 
		  double *left_eigen, double *right_eigen,
		  int workdim, double *workspace)
{

	char jobvl = 'V';
	char jobvr = 'V';
	int info;

	dgeev_ (&jobvl, &jobvr, &dim, matrix, &dim, 
	       eigenvalues_re, eigenvalues_im, 
			 left_eigen, &dim, right_eigen, &dim, 
			 workspace, &workdim, &info);

}


main (int argc, char * argv[]) 
{
	double *mat;
	double *regmat;
	double *ere;
	double *eim;
	double *lev;
	double *rev;
	double *work;
	int dim;
	int workdim;
	int i,j, k;
	
	dim = 28;
	set_npts(4123123);

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

	mat = (double *) malloc (dim*dim*sizeof (double));
	regmat = (double *) malloc (dim*dim*sizeof (double));
	ere = (double *) malloc (dim*sizeof (double));
	eim = (double *) malloc (dim*sizeof (double));
	lev = (double *) malloc (dim*dim*sizeof (double));
	rev = (double *) malloc (dim*dim*sizeof (double));
	workdim = 4*dim*dim;
	work = (double *) malloc (workdim*sizeof (double));

	/* Compute values for the operator */
	for (i=0; i<dim; i++)
	{
		for (j=0; j<dim; j++)
		{
			/* Note transposed matrix'ing for FORTRAN */
			// mat[i+j*dim] = ache_mp(i,j);
			// mat[i+j*dim] = sst(i,j);
			// mat[i+j*dim] = binomial(i,j) * exp (-(i+j)*0.2/dim);
			// mat[i+j*dim] = mtm_svd(i,j);
			// mat[i+j*dim] = mmt_svd(i,j);

			long double re;
			long double im;
			int m = i-dim/2;
			int n = j - dim/2;
			make_elt (m,n, &re, &im);

			mat[i+j*dim] = re;
			// printf ("mat(%d, %d) = %g\n", i,j,mat[i+j*dim]);
		}
		// printf("\n");
	}

	int wd = getworkdim (dim, mat, ere, eim, lev, rev, work);
	printf ("# recommended dim=%d actual dim=%d\n#\n", wd, workdim);
	// workdim = wd;
	
	work = (double *) realloc (work, workdim*sizeof (double));

	// Now, start regulating
	for (t=0.5; t>0.0001; t /= sqrt(sqrt(sqrt(2))))
	{
		for (i=0; i<dim; i++)
		{
			for (j=0; j<dim; j++)
			{
				int m = i-dim/2;
				int n = j - dim/2;
				double reg = exp(-t*(2*m-n)*(2*m-n)/dim*dim);
				regmat[i+j*dim] = reg * mat[i+j*dim];
			}
		}

		// examine the edges of the matrix, see how bad things are.
		double edgemax = 0.0;
		for (i=0; i<dim; i++)
		{
			if (edgemax < fabs(regmat[i])) edgemax = fabs(regmat[i]);
			if (edgemax < fabs(regmat[i+(dim-1)*dim])) edgemax = fabs(regmat[i+(dim-1)*dim]);
			if (edgemax < fabs(regmat[i*dim])) edgemax = fabs(regmat[i*dim]);
			if (edgemax < fabs(regmat[dim-1+i*dim])) edgemax = fabs(regmat[dim-1+i*dim]);
		}

		geteigen (dim, regmat, ere, eim, lev, rev, workdim, work);

		/* print the eigenvalues */
		printf("%g\t%g", t, edgemax);

		for (i=0; i<11; i++)
		{
			printf ("\t%20.15g\t%g", ere[i], eim[i]);
		}
		printf ("\n");
	}
	
}
