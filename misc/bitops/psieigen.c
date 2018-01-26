/*
 * psieigen.c
 *
 * Eigenvalues of downshift from Hessenberg matrix functions.
 * Construct wave functions which put the downshift operator
 * into Hessenberg form.
 *
 * January 2018
 */
#define NOMAIN 1
#include "psi.c"
#include "psibig.c"

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

/*
 * Eigenvalue and eigenvector problem for the downshift operator.
 * GSL version.
 *
 * GSL has a simpler API than Lapack, so try that first. Both seem
 * to agree to six decimal places or better, so I guess its both
 * numerically stable and GSL is returning valid results.
 *
 * Surprisingly, the eigenvalue spectrum lies on a circle in the
 * complex plane.
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

	// Magic incantation to diagonalize the matrix.
	gsl_matrix_view m = gsl_matrix_view_array (matrix, mxi, mxi);
	gsl_vector_complex *eval = gsl_vector_complex_alloc (mxi);
	gsl_matrix_complex *evec = gsl_matrix_complex_alloc (mxi, mxi);
	gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (mxi);
	gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);
	gsl_eigen_nonsymmv_free (w);


#define FP_EIGENFN
#ifdef FP_EIGENFN
	/* Get the eigenvector for the largest eigennvalue, which
	 * should be 1.0.  Graph it to make sure it matches what
	 * we already know it should be.
	 */
	gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);

	gsl_complex eval_0 = gsl_vector_complex_get (eval, 0);

	double x = GSL_REAL(eval_0);
	double y = GSL_IMAG(eval_0);
	double mag = sqrt(x*x+y*y);
	printf ("# eigenvalue = %20.18g + i %20.18g Magnitude=%g\n", x, y, mag);

	/* Get the eigenvector */
	gsl_vector_complex_view evec_0 = gsl_matrix_complex_column (evec, 0);
	double *recof = malloc(mxi*sizeof(double));
	double *imcof = malloc(mxi*sizeof(double));

	gsl_complex z = gsl_vector_complex_get(&evec_0.vector, 0);
	double sgn = GSL_REAL(z);
	if (0 < sgn) sgn = 1.0; else sgn = -1.0;

	mag = 0.0;
	for (int j=0; j< mxi; j++)
	{
		gsl_complex z = gsl_vector_complex_get(&evec_0.vector, j);
		printf("# %d	%g	%g\n", j, GSL_REAL(z), GSL_IMAG(z));
		recof[j] = sgn * GSL_REAL(z);
		imcof[j] = sgn * GSL_IMAG(z);

		mag += recof[j] * recof[j] + imcof[j] * imcof[j];
	}
	mag = sqrt(mag);
	printf("# Vector length = %g\n", mag);

	/* Now graph the eigenvector */
#define NPTS 1201
	for (int n=0; n< NPTS; n++)
	{
		double ex = (((double) n) + 0.5) / ((double) NPTS);
		double y = 0;
		for (int j=0; j< mxi; j++)
		{
			double f = psi_n(ex, K, j);
			y += f * recof[j];
		}
		printf("%d	%g	%g\n", n, ex, y);
	}
#endif

	
#ifdef DECAYING_EXPERIMENT
	for (int i=0; i< mxi; i++)
	{
		gsl_complex eval_i = gsl_vector_complex_get (eval, i);

		double x = GSL_REAL(eval_i);
		double y = GSL_IMAG(eval_i);
		double mag = sqrt(x*x+y*y);
		// printf ("%d	%g	%g	%g\n", i, x, y, mag);
#if 1
		// if (0 == i)
		if (0 !=i && fabs(0.5/K-mag)< 0.01 && 
			// x < 0 && 0 < y && fabs(y) < 0.012) // OK for K=0.6 and 861
			// x < 0 && 0 < y && fabs(y) < 0.0056) // OK for K=0.6 and 461
			// x< 0 && 0 < y && fabs(y) < 0.0063) // OK for K=0.7 and 461
			// x < 0 && 0 < y && fabs(y) < 0.0068) // OK for K=0.8 and 461
			x < 0 && 0 < y && fabs(y) < 0.0047) // OK for K=0.8 and 861
			// x < 0 && 0 < y && fabs(y) < 0.015) // OK for K=0.9 and 461
			// x < 0 && 0 < y && fabs(y) < 0.0016) // OK for K=0.9 and 861
			// Look for period-three eigenvalue
			// 0<y && fabs(y-0.5*sqrt(3.0)/(2.0*K)) < 0.0015)
		{
			gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, i);
			printf ("#\n# eigenvalue = %g + i %g mag=%g\n#\n", x, y, mag);
			double *recof = malloc(mxi*sizeof(double));
			double *imcof = malloc(mxi*sizeof(double));
			for (int j=0; j< mxi; j++)
			{
				gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
				printf("# %d	%g	%g\n", j, GSL_REAL(z), GSL_IMAG(z));
				recof[j] = GSL_REAL(z);
				imcof[j] = GSL_IMAG(z);
			}

#define NPTS 1201
			double sre = 0;
			double sim = 0;
			double sabs = 0;
			double sdot = 0;
			for (int n=0; n< NPTS; n++)
			{
				double ex = (((double) n) + 0.5) / ((double) NPTS);
				double rey = 0;
				double imy = 0;
				double trey = 0;
				double timy = 0;
				double ttrey = 0;
				double ttimy = 0;
				for (int j=0; j< mxi; j++)
				{
					double f = psi_n(ex, K, j);
					rey += f * recof[j];
					imy += f * imcof[j];

					double tex = ex / (2.0*K);
					double tf = psi_n(tex, K, j) + psi_n(0.5+tex, K, j);
					if (K < ex) tf = 0;
					trey += tf * recof[j];
					timy += tf * imcof[j];

					// Supposed to be doubled but its fu'ed up.
					double ttex = tex / (2.0*K);
					double ttf = psi_n(ttex, K, j);
					ttf += psi_n(0.5+ttex, K, j);
					ttf += psi_n(ttex+0.5/(2.0*K), K, j);
					ttf += psi_n(ttex+(0.5/(2.0*K))+0.5, K, j);
					if (K < ex) ttf = 0;
					ttrey += ttf * recof[j];
					ttimy += ttf * imcof[j];

				}
				printf("%d	%g	%g	%g	%g	%g	%g	%g\n",
					n, ex, rey, imy, trey, timy, ttrey, ttimy);

				sre += rey;
				sim += imy;
				sabs += sqrt(rey*rey+imy*imy);
				sdot += rey*imy;
			}
			sre /= NPTS;
			sim /= NPTS;
			sabs /= NPTS;
			sdot /= NPTS;
			printf("#\n# gral-re=%g gral-im=%g gral-abs=%g gral-dot=%g\n#\n",
			     sre, sim, sabs, sdot);
		}
#endif
	}
#endif // DECAYING_EXPERIMENT

	gsl_vector_complex_free(eval);
	gsl_matrix_complex_free(evec);
	free(matrix);
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
	printf("#\n# GSL variant K=%g\n#\n", K);

	// find_midpoints(K);
	big_midpoints(K, 4000, midpoints, MAXN);
	sequence_midpoints(K);
	verify_ortho();
	eigen(K, dim);
}
