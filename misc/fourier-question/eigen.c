/**
 *
 * Eigenvalue and eigenvector finding
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

#include <complex.h>
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
	complex double *eval;
	complex double *lev;
	complex double *rev;
	complex double *work;
	double *w2;
	int dim;
	int workdim;
	int i,j, k;
	
	dim = 28;
	set_npts(123123);

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
	cgeteigen (dim, mat, eval, lev, rev, workdim, work, w2);

	/* ---------------------------------------------- */
	/* print the eigenvalues */
	printf("# Eigenvalues are:\n");
	for (i=0; i<dim; i++)
	{
		printf ("# eigen[%d]=%g +i %g\n", i, creal(eval[i]), cimag(eval[i]));
	}
	printf ("\n\n");
	
	int prtdim = 10;
	if (dim < prtdim) prtdim = dim;
	for (i=0; i<prtdim; i++)
	{
		for (j=0; j<prtdim; j++)
		{
			// printf ("# right %d'th eigenvector[%d]=%g (normalized=%g)\n", 
			//            i,j, rev[j+i*dim], rev[j+i*dim]/rev[i*dim]);
			// printf ("# right %d'th eigenvector[%d]=%g (term ratio=%g)\n", 
			//            i,j, rev[j+i*dim], rev[j+i*dim]/rev[j+1+i*dim]-2.0);
			// printf ("# right %d'th eigenvector[%d]=%g (vec ratio=%g)\n", 
			//            i,j, rev[j+i*dim],  tn*rev[j+i*dim] );
			//
			printf ("# right %d'th eigenvector[%d]=%g +i %g\n", 
			            i,j, creal(rev[j+i*dim]), cimag(rev[j+i*dim]));
		}
		printf ("#\n");
	}
	
	for (i=0; i<prtdim; i++)
	{
		for (j=0; j<prtdim; j++)
		{
			// printf ("# left %d'th eigenvector[%d]=%g (normalized=%g)\n", 
			//            i,j, lev[j+i*dim], lev[j+i*dim]/lev[i*dim]);
			printf ("# left %d'th eigenvector[%d]=%g +i %g)\n", 
			            i,j, creal(lev[j+i*dim]), cimag(lev[j+i*dim]));
		}
		printf ("#\n");
	}
	
#if VALIDATE
	/* ---------------------------------------------- */
	/* Verify i'th eigenvector -- multiply by the matrix, see that 
	 * we get the eigenvector back. */
	int validation_failed = 0;
	for (i=0; i<prtdim; i++)
	{
		/* The j'th element of the i'th eigenvector */
		for (j=0; j<prtdim; j++)
		{
			double sum = 0.0;
			for (k=0; k<dim; k++)
			{
				// sum += ache_mp(j,k) * rev[k+i*dim];
				// sum += mtm_svd(j,k) * rev[k+i*dim];
				sum += mmt_svd(j,k) * rev[k+i*dim];
			}
			sum /= rev[j+i*dim];
			sum /= ere[i];
			sum -= 1.0;
			if (1.0e-12 < fabs(sum))
			{
				validation_failed = 1;
				printf ("# Error: %d'th eigenvec validation failed (should be zero)  [%d]=%g\n", i, j, sum);
			}
		}
		printf ("#\n");
	}
	if (!validation_failed)
		printf("# eigenvec valdidation success\n");
#endif

exit(0);
	/* ---------------------------------------------- */
	/* Print graphable data */
	double y;
	// for (y=1.0; y>=0.0; y-=0.005)
	for (y=0.0; y<=1.0001; y+=0.005)
	{
		// double x = 1.0-y;
		double x = y;
		printf ("%g", x);
		// validate that the zeroth eignevec is 1/(1+x)
		// printf ("\t%g", 1.73205/(1.0+x));
		for (i=0; i<8; i++)
		{
			double yn = 1.0;
			double sum = 0.0;
			double fact = 1.0;
			double norm = 0.0;
			/* The j'th element of the i'th eigenvector */
			for (j=0; j<prtdim; j++)
			{
				norm += rev[j+i*dim];
			}
			for (j=0; j<prtdim; j++)
			{
				// Factorial not needed, ache already has the factorial folded in.
				// See, for example, the derives of seroth eigenvec. 
				// sum += yn * rev[j+i*dim] / fact;
				// sum += yn * rev[j+i*dim];
				sum += yn * rev[j+i*dim] / norm;
				// printf ("duuude j=%d fact=%g yn=%g\n", j, fact, yn);
				yn *= y;
				fact *= j+1;
			}
			printf("\t%g", sum);
		}
		printf ("\n");
	}

	/* ---------------------------------------------- */
	/* Print the eigenvectors again, normalized correctly, so we can graph. */
	/* This time, they are printed so that the i'th eigenvector is in the
	   i'th column, and each row is a component of the eigenvector */

exit(0);
	printf ("# ------------------------------------------- \n");
	double norm[50];
	for (i=0; i<prtdim; i++)
	{
		norm[i] = 0.0;
		/* The j'th element of the i'th eigenvector */
		for (j=0; j<prtdim; j++)
		{
			norm[i] += rev[j+i*dim];
		}
		// norm[i] = 31.0 * lev[30+i*dim];
		// norm[i] = (1<<25) *rev[25+i*dim];
	}

	for (j=0; j<prtdim; j++)
	{
		printf ("%d\t", j);
		for (i=0; i<13; i++)
		{
			printf ("%g\t", rev[j+i*dim]/norm[i]);
			// printf ("%g\t", lev[j+i*dim]/norm[i]);
		}
		printf ("\n");
	}
}
