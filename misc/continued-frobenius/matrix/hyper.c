/*
 * Eigenvalue finding
 * Find eigenvalues of the simple hyperbolic
 * composition operator
 *
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
 * Linas Vepstas December 2005
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "lapack.h"

int 
getworkdim (int dim, 
        double *matrix,  /* double precision complex */
        double *eigenvalues,
		  double *left_vec, double *right_vec,
		  double *workspace_a, double *workspace_b)
{
	char jobvl = 'V';
	char jobvr = 'V';
	int workdim;
	int info;

	workdim = -1;

	zgeev_ (&jobvl, &jobvr, &dim, matrix, &dim, 
	       eigenvalues, 
			 left_vec, &dim, right_vec, &dim, 
			 workspace_a, &workdim, workspace_b, &info);

	return (int) (workspace_a[0]+0.5);
}

int 
get_complex_eigen (int dim, 
        double *matrix,  /* double precision complex */
        double *eigenvalues,
		  double *left_vec, double *right_vec,
		  int workdim, double *workspace_a, double *workspace_b)
{

	char jobvl = 'V';
	char jobvr = 'V';
	int info;

	zgeev_ (&jobvl, &jobvr, &dim, matrix, &dim, 
	       eigenvalues, 
			 left_vec, &dim, right_vec, &dim, 
			 workspace_a, &workdim, workspace_b, &info);

	printf ("# duuude info = %d\n", info);
}

/* --------------------------------------------------- */

complex double in (int n)
{
	complex double ii = I;

	if (0 > n)
	{
		n = -n;
		ii = -I;
	}
	if (n%4 == 0) return 1.0;
	if (n%4 == 1) return 0.0 + ii;
	if (n%4 == 2) return -1.0;
	if (n%4 == 3) return 0.0 -ii;
	return 0.0;
}

complex double nin (int n)
{
	complex double ii = -I;

	if (0 > n)
	{
		n = -n;
		ii = I;
	}
	if (n%4 == 0) return 1.0;
	if (n%4 == 1) return 0.0 + ii;
	if (n%4 == 2) return -1.0;
	if (n%4 == 3) return 0.0 -ii;
	return 0.0;
}

complex double nn (int n)
{
	if (0 > n) n = -n;

	if (n%2 == 0) return 1.0;
	return -1.0;
}

complex double compose (int m, int n)
{
	complex double ret;

	if (m == n)
	{
		if (0 == m) return (complex double) 1.0;
		ret = 0.25 + (in(n) - 1.0)/ (0.0 + 3.0*M_PI*I*n);
		return ret;
	}
	if (m == 2*n)
	{
		ret = 0.25 + (nin(n) - 1.0)/ (0.0 + 6.0*M_PI*I*n);
		return ret;
	}
	if (2*m == n)
	{
		ret = 0.5 - (nin(n) - 1.0)/ (0.0 + 3.0*M_PI*I*n);
		return ret;
	}

	ret = 3.0 * nn(n+m) * (n-m) - in(m)*(n-2.0*m) - in(n) * (2.0*n-m);
	ret *= nn(m+n) * n;
	ret /= 2 * M_PI * I * (2*n-m)*(n-m)*(n-2.0*m);

	return ret;
}
  
/* --------------------------------------------------- */

int main (int argc, char * argv[]) 
{
	double *mat;
	double *val;
	double *lev;
	double *rev;
	double *worka, *workb;
	int dim;
	int workdim;
	int i,j, k;
	
	dim = 9;
	int nc = - (dim-1)/2;

	printf ("#\n#\n");
	printf ("# Eigenvectors of hyperbolic composition \n");
	printf ("# Numerically solved to rank=%d\n", dim);
	printf ("#\n#\n");

	mat = (double *) malloc (2*dim*dim*sizeof (double));
	val = (double *) malloc (2*dim*sizeof (double));
	lev = (double *) malloc (2*dim*dim*sizeof (double));
	rev = (double *) malloc (2*dim*dim*sizeof (double));
	workdim = 8*dim*dim;
	worka = (double *) malloc (workdim*sizeof (double));
	workb = (double *) malloc (workdim*sizeof (double));

	/* Insert values for the GKW operator at x=1 (y=1-x) */
	for (i=0; i<dim; i++)
	{
		for (j=0; j<dim; j++)
		{
			/* Note transposed matrix'ing for FORTRAN */
			complex double v = compose (i+nc,j+nc);
			mat[2*(i+j*dim)] = creal(v);
			mat[2*(i+j*dim)+1] = cimag(v);
			printf ("# m[%d, %d] = %g + i %g\n", i+nc, j+nc, creal (v), cimag (v));
		}
	}


	int wd = getworkdim (dim, mat, val, lev, rev, worka, workb);
	printf ("# recommended dim=%d actual dim=%d\n#\n", wd, workdim);
	// workdim = wd;
	
	// work = (double *) realloc (work, workdim*sizeof (double));
	get_complex_eigen (dim, mat, val, lev, rev, workdim, worka, workb);

	/* ---------------------------------------------- */
	/* print the eigenvalues */
	for (i=0; i<dim; i++)
	{
		double re = val[2*i];
		double im = val[2*i+1];
		double md = sqrt (re*re+im*im);
		printf ("# eigen[%d]=%20.15g +i %g = |%g|\n", i, re, im, md);
	}
	printf ("\n\n");
	
	int prtdim = 9;
	if (prtdim > dim) prtdim = dim;
	for (i=0; i<prtdim; i++)
	{
		for (j=0; j<prtdim; j++)
		{
			printf ("# right %d'th eigenvector[%d]=%g +i %g\n", 
			            i,j, rev[2*(j+i*dim)], rev[2*(j+i*dim)+1]);
		}
		printf ("#\n");
	}
	
#if 0
	for (i=0; i<prtdim; i++)
	{
		for (j=0; j<prtdim; j++)
		{
			printf ("# left %d'th eigenvector[%d]=%g +i %g\n", 
			            i,j, lev[2*(j+i*dim)], lev[2*(j+i*dim)+1]);
		}
		printf ("#\n");
	}
	
	/* ---------------------------------------------- */
	/* Verify i'th eigenvector */
	for (i=0; i<prtdim; i++)
	{
		/* The j'th element of the i'th eigenvector */
		for (j=0; j<prtdim; j++)
		{
			complex double sum = 0.0;
			for (k=0; k<dim; k++)
			{
				complex double v = rev[2*(k+i*dim)] + I * rev[2*(k+i*dim)+1];
				sum += compose(j+nc,k+nc) * v;
			}
			double re = val[2*i];
			double im = val[2*i+1];
			complex double eiv = re + I* im;
			sum /= eiv;
			printf ("# %d'th eigenvec validation [%d]=%g+ i%g\n", 
				i, j, creal(sum), cimag(sum));
		}
		printf ("#\n");
	}
#endif

	/* ---------------------------------------------- */
	/* Print graphable data */
	i = 1;  // which egenvector
	double x;
	for (x=0.0; x<1.0; x+=0.02)
	{
		printf ("%g", x);

		double complex sum = 0.0;
		/* The j'th element of the i'th eigenvector */
		for (j=0; j<dim; j++)
		{
			int k = j+nc;
			complex double ex = cexp (2.0*M_PI*I*k*x);
			complex double v = rev[2*(j+i*dim)] + I * rev[2*(j+i*dim)+1];
			sum += v*ex;
		}
		printf("\t%g\t%g\n", sum);
	}
}
