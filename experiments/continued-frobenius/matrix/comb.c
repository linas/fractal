
/* 
 * build and solve comb-like eignevalue equation
 *
 * Linas Vepstas September 2004
 */

#include <math.h>
#include <stdio.h>
#include "lapack.h"

// LAPACK API wrappers -----------------------------------------
int 
getworkdim (int dim, double *matrix,
        double *eigenvalues_re, double *eigenvalues_im, 
		  double *left_eigen, double *right_eigen,
		  double *workspace)
{

	char jobvl = 'N';
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

	char jobvl = 'N';
	char jobvr = 'V';
	int info;

	dgeev_ (&jobvl, &jobvr, &dim, matrix, &dim, 
	       eigenvalues_re, eigenvalues_im, 
			 left_eigen, &dim, right_eigen, &dim, 
			 workspace, &workdim, &info);

	if (info) printf ("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxx info=%d\n", info);
	return info;
}

// ---------------------------------------------------------

/*
 * Build a matrix corresponding to the differential equation
 * y^4r'' + y^3r' (2-4ay(y-b)) + y^3r[4a^2y(y-b)^2 -2ay -4a(y-b)] + rsin(2piy)
 * where
 * the matrix expresses this in terms of the taylor expansion, 
 * with r(y) = sum a_n y^n
 */

double
get_matrix_elt (int p, int n, double a, double b)
{
	double elt = 0.0;
	int k;

	// multiply by sin (2piy)
	k = p-n;
	// whoops, this if statment is for cos (2pi y)
	// if ((k>=0) && (0 == k%2))
	if ((k>0) && (1 == k%2))
	{
		int i;
		double term = 1.0;
		if (1 == k%2) term = -1.0;
		for (i=0; i<k; i++)
		{
			term *= 2.0*M_PI / ((double) (i+1));
		}
		elt += term;
	}
	
	if (p == n+2)
	{
		elt += n*(n-1);   // y^4 r''
		elt += 2*n;       // y^3 r' 2		
	}
	else if (p == n+3)
	{
		elt += 4.0*a*b* ((double) n);  // y^3r' 4ayb
		elt += 4.0*a*b;                // y^3r (-4ab)
	}
	else if (p == n+4)
	{
		elt -= 4.0*a* ((double) n);  // y^3r' (-4ay^2)
		elt += 4.0*a*a*b*b;          // y^3 r  4a^2y b^2
		elt -= 2.0*a;                // y^3 r  (-2ay)
		elt -= 4.0*a;                // y^3 r  (-4ay)
	}
	else if (p == n+5)
	{
		elt -= 8.0*a*a*b;       // y^3 r 4a^2y (-2by)
	}
	else if (p == n+6)
	{
		elt += 4.0*a*a;       // y^3 r 4a^2y y^2
	}

	return elt;
}

void
fill_matrix (double a, double b, int dim, double *mat)
{
	int i,j;

	double regulate;
	regulate = (double) dim;
	regulate = log(1.0e-20) / (regulate*regulate);

	// build the matrix
	for (i=0; i<dim; i++)
	{
		double cutoff = exp (regulate *((double) i)*((double) i));
		for (j=0; j<dim; j++)
		{
			double v = get_matrix_elt (i,j, a, b);
			v *= cutoff;
			// use the fortran-style matrix conventions
			mat [i+dim*j] = v;
			printf ("%g\t", v);
		}
		printf ("\n");
	}
}

main ()
{
	int i, j;
	int dim;

	double *taylor;   // the matrix
	double *revec;    // right eigenvector
	double *levec;    // left eigenvector
	double *reev, *imev;  // real, imaginary eigenvalues

	int workdim;
	double *work;

	dim = 13;

	taylor = (double *) malloc (dim*dim*sizeof(double));
	revec = (double *) malloc (dim*dim*sizeof(double));
	levec = (double *) malloc (dim*dim*sizeof(double));
	reev = (double *) malloc (dim*sizeof(double));
	imev = (double *) malloc (dim*sizeof(double));
	
	workdim = 4*dim*dim;
	work = (double *) malloc (workdim*sizeof(double));
	
	// wild values are a=0.16 b=6.0
	double a = 0.16;
	double b = 6.0;

// for (b=6.0; b<16.0; b+=1.0) {
for (a=0.1; a<6.0; a *= 1.1) {
	// build the matrix
	fill_matrix (a,b,dim,taylor);

	// workdim = getworkdim (dim, taylor, reev, imev, levec, revec, work);
	geteigen (dim, taylor, reev, imev, levec, revec, workdim, work);

#if 1
	// print the eigenvalues and the eigenvectors
	int prtdim=dim;
	for (i=0; i<prtdim; i++)
	{
		printf ("ev %g  %g  \n", reev[i], imev[i]);
#if 0
		for (j=0; j<prtdim; j++)
		{
			printf ("%g\n", revec[j+dim*i]);
		}
		printf ("\n");
#endif
	}

	// validate the eigenvalues and vectors
	fill_matrix (a,b,dim,taylor);

	// the i'th eigenvector
	for (i=dim-1; i> dim-prtdim; i--)
	{
		printf ("eigenvalue = %g %g\n", reev[i], imev[i]);
		// check the j'th element of the i'th eigenvector
		for (j=0; j<prtdim; j++)
		{
			double sum=0.0;
			int k;
			for (k=0; k<dim; k++)
			{
				sum += taylor[j+dim*k] * revec[k+dim*i];
			}
			sum /= reev[i] * revec[j+dim*i];
			printf ("vec %d comp%d evec=%g  \t m*evec=%g\n", i, j, revec[j+dim*i], sum);
		}
		printf ("\n");
	}
#endif

	// find the largest eigenvalue
	int ibig = 0;
	double ebig = -1.0e100;
	for (i=0; i<dim; i++) 
	{
		if (reev[i] > ebig) { ebig = reev[i]; ibig = i; }
	}
	printf ("biggest eigenvalue is at %d val=%g\n", ibig, ebig);

	// now evaluate the first eigenvector for 'large'  values of x:
	double y = 50.0;
	double yn = 1.0;
	double sum = 0.0;
	for (j=0; j<dim; j++)
	{
		sum += revec[j+ibig*dim] * yn;
		yn *= y;
	}
	sum *= exp (-a*(y-b)*(y-b));

	printf ("duude a=%g b=%g sum=%g\n", a,b,sum);
}
	
}
