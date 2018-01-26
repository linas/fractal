/*
 * psiblas.c
 *
 * Eigenvalues of downshift from Hessenberg matrix functions.
 * Uses LAPACK to find these numerically.
 *
 * Construct wave functions which put the downshift operator
 * into Hessenberg form.
 *
 * January 2018
 */
#define NOMAIN 1
#include "psi.c"
#include "psibig.c"

#include <lapacke.h>
#include <cblas.h>

/*
 * Eigenvalue and eigenvector problem. -- LAPACK version.
 * This computes the eigenvalues, only, for the dim x dim submatrix
 * of the downshift operator.
 *
 * One might think that the eigenvalues are going to be real. But they
 * are not. They are complex, and seem to lie on a unit circle of
 * radius 1/2K.  Wow!
 *
 * This version (LAPACK) gives the same results as the GSL version.
 * (to at least 6 decimal places). (For the ROW_MAJOR version.)
 * So this implies that GSL is not insane (like I thought it was).
 *
 * Lapack, Eigenvalues only, use: HSEQR
 * Eigenvalues an eigenvectors: HSEQR TREVC
 * DHSEQR for double precision, ZHSEQR for complex double precision.
 *
 * This does eigenvalues only, right now.
 */
void eigen(double K, int dim)
{
	int mxi = MAXN-1;
mxi = dim;
	double* matrix = malloc(mxi*mxi*sizeof(double));

	// This fills matrix in row-major order.
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

	/*
	 * lapack_int LAPACKE_dhseqr(int matrix_layout,
	 *                           char job, char  compz,
	 *                           lapack_int n,
	 *                           lapack_int ilo, lapack_int ihi, double* h,
	 *                           lapack_int ldh, double* wr, double* wi,
	 *                           double* z, lapack_int ldz );
	 */

	double* wr = malloc(mxi* sizeof(double));
	double* wi = malloc(mxi* sizeof(double));
	double* z = malloc(mxi*mxi* sizeof(double));
	lapack_int info;
	info = LAPACKE_dhseqr(LAPACK_ROW_MAJOR, 'E', 'N', mxi, 1, mxi, matrix,
	                      mxi, wr, wi, z, mxi);
	printf("# info=%d  K=%g\n", info, K);
	double avg = 0.0;
	double cnt = 0.0;
	for (int i=0; i< mxi; i++)
	{
		double mag = sqrt(wr[i]*wr[i] + wi[i]*wi[i]);
		printf("%d	%g	%g	%g\n", i, wr[i], wi[i], mag);
		double a = avg/cnt;
		if (i<20 && 0.98*a < mag && mag < 1.02*a) {avg += mag; cnt += 1; }
		else if (0.5 < mag && mag < 0.999) {avg += mag; cnt += 1; }
	}
	printf("#\n# average radius= %g cnt= %g\n#\n", avg/cnt, cnt);
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

	// find_midpoints(K);
	big_midpoints(K, 9000, midpoints, MAXN);
	sequence_midpoints(K);
	verify_ortho();
	eigen(K, dim);
}
