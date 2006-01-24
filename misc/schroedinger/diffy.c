
/*
 * Eigenvalue and eigenvector finding
 * Call lapack from C code.
 *
 * Use LAPACK xSTEGR routine to find eigenvalues of tri-diagonal matrix.
 * specifically, use DSTEGR -- real double precision:
 *
 * Apply to solve 1D Schroedinger equation
 *
 * Linas Vepstas September 2004
 * Linas Vepstas January 2006
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "lapack.h"


#if 0
/* dry run -- get the working dimension */
int 
getworkdim (int dim, double *diags,
        double *subdiags, 
        double abstol,
		  double *workspace)
{
	char jobz = 'V';
	char range = 'A';
	double vu=0.0, vl=1.0E30;
	int il=0, iu=12345670;
	int m;
	int workdim;
	int info;

	workdim = -1;

	dstegr_ (&jobz, &range, &dim, diags, subdiags,
	         &vl, &vu, &il, &iu, &abstol,
            &m, eigenvals, eigenvecs, &dim, 
	         suppport, rwork, &rworkdim, iwork, &iworkdim,
			   &info);
}
#endif

int 
geteigen (int dim, double *diags,
        double *subdiags, 
        double abstol,
        double *eigenvals,
        double *eigenvecs,
		  int *support)
{
	char jobz = 'V';
	char range = 'A';
	double vu=0.0, vl=1.0E30;
	int il=0, iu=12345670;
	int m;
	int rworkdim, iworkdim;
	double *rwork;
	int * iwork;
	int info;

	rworkdim = 20*dim;
	iworkdim = 10*dim;

	rwork = (double *) malloc (rworkdim*sizeof (double));
	iwork = (int *) malloc (iworkdim*sizeof (int));

	dstegr_ (&jobz, &range, &dim, diags, subdiags,
	         &vl, &vu, &il, &iu, &abstol,
            &m, eigenvals, eigenvecs, &dim, 
	         support, rwork, &rworkdim, iwork, &iworkdim,
			   &info);

	free (iwork);
	free (rwork);

	return info;
}


double 
kinetic_diag (double step)
{
	return 1.0/(step*step);
}

double 
kinetic_subdiag (double step)
{
	return -0.5/(step*step);
}

/* --------------------------------------------------- */

int centered_quadratic_dim (int kstep, int decimal_prec)
{
	double Npts = (2*kstep+1)*sqrt (2*decimal_prec*log(10.0));
	int dim = 2*Npts+1;
	return dim;
}

/* --------------------------------------------------- */
double 
galois_potential (int n, double step)
{
	double y = n*step;
	double t = y*y+0.75;
	double s = y*y-0.25;
	// return t*t*t/(s*s);
	return 0.5*t*t*t/(s*s);
	// return 0.5*y*y;  // plain old harmonic osc for comparison
}

double 
double_potential (int n, double step)
{
	double y = (n+0.5)*step;
	return 0.5*(y*y + 1.0/(y*y));
}

void galois_data (int kstep, int dim, 
           double *diags, double *subdiags,
           int *ncenter, double *pdelta)
{
	int i;
	int Npts = (dim-1)/2;
	double delta = 1.0/(2*kstep+1);
	*ncenter = Npts;
	*pdelta = delta;

	printf ("#\n#\n");
	printf ("# Eigenvectors of the Galois Schroedinger\n");
	printf ("# Numerically solved on %d lattice points, step=%g\n", dim, delta);
	printf ("#\n#\n");

	/* Insert values  */
	for (i=0; i<dim; i++)
	{
		diags[i] = kinetic_diag(delta);
		diags[i] += galois_potential (i-Npts, delta);
		// diags[i] += galois_potential (i+kstep+1, delta);
		// diags[i] += potential (i, delta);

		subdiags[i] = kinetic_subdiag(delta);
	}
}

/* --------------------------------------------------- */
double 
elliptic_potential (int n, double step)
{
	double y = fabs(n*step);
	double g2 =3.0;
	double g3=1.0;
	double x = sqrt (y*y*y +g2*y+g3);
	return x;
}

void elliptic_data (int kstep, int dim, 
           double *diags, double *subdiags,
           int *ncenter, double *pdelta)
{
	int i;
	int Npts = (dim-1)/2;
	double delta = 1.0/(2*kstep+1);
	*ncenter = Npts;
	*pdelta = delta;

	printf ("#\n#\n");
	printf ("# Eigenvectors of Elliptic\n");
	printf ("# Numerically solved on %d lattice points, step=%g\n", dim, delta);
	printf ("#\n#\n");

	/* Insert values  */
	for (i=0; i<dim; i++)
	{
		diags[i] = kinetic_diag(delta);
		diags[i] += elliptic_potential (i-Npts, delta);
		subdiags[i] = kinetic_subdiag(delta);
	}
}

/* --------------------------------------------------- */
main (int argc, char * argv[]) 
{
	double *diags;
	double *subdiags;
	double *eigenvals;
	double *eigenvecs;
	int * support;
	int i,j, k;
	
	int kstep = 10;

	if (argc < 3)
	{
		fprintf (stderr, "Usage: %s <kstep> <eigen_prt>\n", argv[0]);
		exit (-1);
	}
	kstep = atoi (argv[1]);
	int eigen_prt = atoi(argv[2]);

	int Mprec = 9;
	double abstol=1.0e-6;
	int dim = centered_quadratic_dim (kstep, Mprec);

	diags = (double *) malloc (dim*sizeof (double));
	subdiags = (double *) malloc (dim*sizeof (double));
	eigenvals = (double *) malloc (dim*sizeof (double));
	eigenvecs = (double *) malloc (dim*dim*sizeof (double));
	support = (int *) malloc (2*dim*sizeof (int));

	double delta;
	int ncenter;
	// galois_data (kstep, dim, diags, subdiags, &ncenter, &delta);
	elliptic_data (kstep, dim, diags, subdiags, &ncenter, &delta);

	geteigen (dim, diags, subdiags, abstol,
        eigenvals, eigenvecs, support);

	/* ---------------------------------------------- */

	/* print the eigenvalues */
	int prtdim = 96;
	if (dim < prtdim) prtdim = dim;
	for (i=0; i<prtdim; i++)
	{
		// double eddie = eigenvals[i+1]/eigenvals[i] - 1.0;
		double eddie = eigenvals[i+1] -eigenvals[i];
		printf ("# eigen[%d]=%20.15g  diff=%8.5f\n", 
		       i, eigenvals[i], eddie);
	}
	printf ("\n\n");
	
	/* Print the eigenvectors */
	/* Note transposed matrix'ing for FORTRAN */
	for (j=0; j<dim; j++)
	{
		printf ("%d	%g", j, (j-ncenter)*delta);
		for (i=eigen_prt; i<9+eigen_prt; i++)
		{
			printf ("\t%g", eigenvecs[j+i*dim]);
		}
		printf ("\n");
	}
	
}
