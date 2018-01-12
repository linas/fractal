/*
 * psiblas.c
 *
 * Eigenvalues of downshift from Hessenberg matrix functions.
 * Construct wave functions which put the downshift operator
 * into Hessenberg form.
 *
 * January 2018
 */
#define NOMAIN 1
#include "psi.c"

/* 
 * Eigenvalue and eigenvector problem.
 * Lapack, Eigenvalues only, use: HSEQR
 * Eigenvalues an eigenvectors: HSEQR TREVC
 * DHSEQR for double precision, ZHSEQR for complex double precision.
 *
 * We expect eigenvalues to be real. So don't screw with the complex.
 */
void eigen(double K, int dim)
{
	int mxi = MAXN-1;
mxi = dim;
	double* matrix = malloc(mxi*mxi*sizeof(double));

	// Do we need this, or the transpose ??? Does it matter?
	int idx = 0;
	for (int i=0; i< mxi; i++)
	{
		for (int j=0; j< mxi; j++)
		{
			double h = hess(K, i, j);
			matrix[idx] = h;
			idx++;
		}
	}
}


int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s K dim\n", argv[0]);
		exit(1);
	}
	double K = atof(argv[1]);
	int dim = atoi(argv[2]);
	printf("#\n# K=%g\n#\n", K);

#if 0
#define NPTS 201
	double s = 0.0;
	for (int i=0; i< NPTS; i++)
	{
		double x = ((double) i + 0.5) / ((double) NPTS);
		double (*f)(double, double) = psi_1;
		double ef = f(x, K);
		double p0 = part(x, K, f, 0);
		double p1 = part(x, K, f, 1);
		double y = xfer(x, K, f);
		s += y / ((double) NPTS);
		printf("%d	%g	%g %g	%g	%g	%g\n", i, x, ef, p0, p1, y, s);
	}
#endif

	find_midpoints(K);
	verify_ortho();
	eigen(K, dim);
}
