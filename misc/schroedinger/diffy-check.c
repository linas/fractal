
/*
 * Eigenvalue and eigenvector finding
 * Call lapack from C code.
 *
 * Use LAPACK DGEEV routine to find eigenvalues of tri-diagonal matrix.
 *
 * Apply to solve 1D Schroedinger equation
 * Used to sanity-echeck the trigiagonal solver
 *
 * Linas Vepstas September 2004
 * Linas Vepstas January 2006
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "lapack.h"


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

double 
kinetic (int i, int j, double step)
{
	if (i<j-1) return 0.0;
	if (i>j+1) return 0.0;
	if (i==j) return 1.0/(step*step);
	return -0.5/(step*step);
}

double 
potential (int n, double step)
{
	double y = n*step;
	double t = y*y+0.75;
	double s = y*y-0.25;
	return 0.5*t*t*t/(s*s);
	// return 0.5*y*y;  // plain old harmonic osc for comparison
}

main (int argc, char * argv[]) 
{
	double *mat;
	double *ere;
	double *eim;
	double *lev;
	double *rev;
	double *work;
	int workdim;
	int i,j, k;
	
	int kstep = 10;

	if (argc < 2)
	{
		fprintf (stderr, "Usage: %s <kstep>\n", argv[0]);
		exit (-1);
	}
	kstep = atoi (argv[1]);

	int Mprec = 6;

	double Npts = (2*kstep+1)*sqrt (2*Mprec*log(10.0));
	double delta = 1.0/(2*kstep+1);

	int dim = 2*Npts+1;

	printf ("#\n#\n");
	printf ("# Eigenvectors of the Galois Schroedinger\n");
	printf ("# Numerically solved on %d lattice points, step=%g\n", dim, delta);
	printf ("#\n#\n");

	mat = (double *) malloc (dim*dim*sizeof (double));
	ere = (double *) malloc (dim*sizeof (double));
	eim = (double *) malloc (dim*sizeof (double));
	lev = (double *) malloc (dim*dim*sizeof (double));
	rev = (double *) malloc (dim*dim*sizeof (double));
	workdim = 4*dim*dim;
	work = (double *) malloc (workdim*sizeof (double));

	/* Insert values for the GKW operator at x=1 (y=1-x) */
	for (i=0; i<dim; i++)
	{
		for (j=0; j<dim; j++)
		{
			/* Note transposed matrix'ing for FORTRAN */
			mat[i+j*dim] = kinetic(i,j, delta);
		}
		mat [i+i*dim] += potential (i-Npts, delta);
	}

	int wd = getworkdim (dim, mat, ere, eim, lev, rev, work);
	printf ("# recommended workdim=%d actual dim=%d\n#\n", wd, workdim);
	workdim = wd+10;
	
	work = (double *) realloc (work, workdim*sizeof (double));
	geteigen (dim, mat, ere, eim, lev, rev, workdim, work);

	/* ---------------------------------------------- */
	/* buble sort the eignes */
	for (i=0; i<dim; i++) {
		for (j=0; j<dim; j++) {
			if (ere[i]< ere[j]) {
				double tmp = ere[i];
				ere [i] = ere[j];
				ere[j] = tmp;
			}
		}
	}

	/* print the eigenvalues */
	for (i=0; i<dim; i++)
	{
		printf ("# eigen[%d]=%20.15g  diff=%20.15g\n", i, ere[i], ere[i+1]-ere[i]);
	}
	printf ("\n\n");
	
#if LATER
	int prtdim = 36;
	if (dim < prtdim) prtdim = dim;
	for (i=0; i<prtdim; i++)
	{
		double tn = 1.0;
		double thrn = 1.0;
		for (j=0; j<prtdim; j++)
		{
			// printf ("# right %d'th eigenvector[%d]=%g (normalized=%g)\n", 
			//            i,j, rev[j+i*dim], rev[j+i*dim]/rev[i*dim]);
			// printf ("# right %d'th eigenvector[%d]=%g (term ratio=%g)\n", 
			//            i,j, rev[j+i*dim], rev[j+i*dim]/rev[j+1+i*dim]-2.0);
			// printf ("# right %d'th eigenvector[%d]=%g (vec ratio=%g)\n", 
			//            i,j, rev[j+i*dim],  tn*rev[j+i*dim] );
			//
			// double r1 = 2.0 * rev[j+1+i*dim]/rev[j+i*dim] - 1.0;
			// double r2 = 2.0 * rev[j+2+i*dim]/rev[j+1+i*dim] - 1.0;
			// double r = r1/r2;
			// r *= j*j*j*j*log(log (log (log (j+1))));
			// r *= j;
		}
		printf ("#\n");
	}
	
	for (i=0; i<prtdim; i++)
	{
		for (j=0; j<prtdim; j++)
		{
			// printf ("# left %d'th eigenvector[%d]=%g (normalized=%g)\n", 
			//            i,j, lev[j+i*dim], lev[j+i*dim]/lev[i*dim]);
			printf ("# left %d'th eigenvector[%d]=%g (ratio=%g)\n", 
			            i,j, lev[j+i*dim], ((j+1)*lev[j+i*dim])/((j+2)*lev[j+1+i*dim]));
		}
		printf ("#\n");
	}
	
	/* ---------------------------------------------- */
	/* Verify i'th eigenvector -- multiply by the matrix, see that 
	 * we get the eigenvector back. */
	for (i=0; i<prtdim; i++)
	{
		/* The j'th element of the i'th eigenvector */
		for (j=0; j<prtdim; j++)
		{
			double sum = 0.0;
			for (k=0; k<dim; k++)
			{
				sum += ache_mp(j,k) * rev[k+i*dim];
			}
			sum /= rev[j+i*dim];
			sum /= ere[i];
			sum -= 1.0;
			printf ("# %d'th eigenvec validation [%d]=%g\n", i, j, sum);
		}
		printf ("#\n");
	}

	/* ---------------------------------------------- */
	/* Print graphable data */
	double y;
	for (y=1.0; y>=0.0; y-=0.005)
	{
		double x = 1.0-y;
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
		norm[i] = (1<<25) *rev[25+i*dim];
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
#endif
}
