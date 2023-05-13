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
 * specifically, real double precision:
 * DGEEV
 *
 * Apply to find eigenvalues of the Gauss-Kuz'min-Wirsing (GKW) operator 
 *
 * Linas Vepstas September 2004
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ache.h"
#include "zetafn.h"
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

/* kinetic part only */
double 
kino (int m, int n)
{
	if (m != n-2) return 0.0;
	return (double) (n*(n+1));
}

double sst (int i, int j)
{
	double lam;
	if (i>j) return 0.0;
	lam = pow (2.0, (double) (j+1)) -1.0;
	double bin = binomial (j,i);
	printf ("mat(%d, %d) = %g/%g\n", i,j,bin,lam);
	return bin /lam;
}

/* Singular Value Decomposition M^transpose M
 * where M is the binomial transform
 */
double mtm_svd (int m, int k)
{
	if (k == m)
	{
		double mtm = binomial (2*m,m);
		mtm += -1.0 + 1.0/pow(4.0, m);
		return mtm;
	}
	if (k>m) {int j = m; m=k; k=j; }
	double mtm = binomial (m+k, k);
	mtm -= binomial (m,k) * (1.0 - 1.0/ pow(2.0, m+k+1));
	return mtm;
}

/** Bernoulli process, doubling */
double bern (int m, int k)
{
	if (k<m) return 0.0;
	if (k == m) return pow(2.0, -k);
	return binomial(k,m) * pow(2.0, -k-1);
}

double mmt_svd (int m, int k)
{
	int p=0;
	double acc=0.0;
	double term;
	do
	{
		term = bern(m,p) * bern(k,p);
		acc += term;
		// printf ("yo m=%d k=%d p=%d term=%g acc=%g\n", m,k,p,term,acc);
		p++;
	} while ((1.0e-14 < term) || (acc == 0.0));
	return acc;
}

void prtbin (double x)
{
	if (0.0 > x) x = -x;
	int i;
	for (i=0; i<45; i++)
	{
		if (0.5 < x) {
			printf ("1");
			x -= 0.5;
		} else {
			printf ("0");
		}
		x *= 2.0;
	}
}

main (int argc, char * argv[]) 
{
	double *mat;
	double *ere;
	double *eim;
	double *lev;
	double *rev;
	double *work;
	int dim;
	int workdim;
	int i,j, k;
	
	dim = 28;

	if (argc < 2)
	{
		fprintf (stderr, "Usage: %s <dim>\n", argv[0]);
		exit (-1);
	}
	dim = atoi (argv[1]);

	printf ("#\n#\n");
	printf ("# Eigenvectors of the GKW (Gauss Kuz'min Wirsing) Operator\n");
	printf ("# or maybe something else, depending on the version of this code\n");
	printf ("# Numerically solved to rank=%d\n", dim);
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
			// mat[i+j*dim] = ache_mp(i,j);
			// mat[i+j*dim] = sst(i,j);
			// mat[i+j*dim] = binomial(i,j) * exp (-(i+j)*0.2/dim);
			// mat[i+j*dim] = mtm_svd(i,j);
			mat[i+j*dim] = mmt_svd(i,j);
			printf ("mat(%d, %d) = %g\n", i,j,mat[i+j*dim]);
		}
		printf("\n");
	}

#if POORLY_CONDITIONED_BOMB_OUT_KINETIC
	/* Add in the kinetic term ! */
	for (i=0; i<dim; i++)
	{
		for (j=0; j<dim; j++)
		{
			/* Note transposed matrix'ing for FORTRAN */
			mat[i+j*dim] -= kino(i,j);
		}
	}
#endif
	
	int wd = getworkdim (dim, mat, ere, eim, lev, rev, work);
	printf ("# recommended dim=%d actual dim=%d\n#\n", wd, workdim);
	// workdim = wd;
	
	work = (double *) realloc (work, workdim*sizeof (double));
	geteigen (dim, mat, ere, eim, lev, rev, workdim, work);

	/* ---------------------------------------------- */
	/* print the eigenvalues */
	printf("# Eigenvalues are:\n");
	for (i=0; i<dim; i++)
	{
		printf ("# eigen[%d]=%20.15g +i %g  ratio=%20.15g\n", i, ere[i], eim[i], ere[i]/ere[i+1]);
	}
	printf ("\n\n");
	
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
#if RECURRSION_RELATION
			double r = rev[j+1+i*dim]- 0.5*rev[j+i*dim];
			r /= rev[j-1+i*dim];
			// r *= j*j*j* log(j+1) * log(j+1) * log(j+1);
			r = log (r);
			r /= j;
			printf ("# right %d'th eigenvector[%d]=%g (term log ratio=%g)\n", 
			           i,j, rev[j+i*dim], r);
#endif
#if WHO_KNOWS_WHAT
			double r = rev[j+i*dim] / (rev[30+i*dim] * (1<<30));
			r *= tn;
			r -= 1.0;
			r *= thrn;
			printf ("# right %d'th eigenvector[%d]=%g (guh %g\n", 
			            i,j, rev[j+i*dim], r);
#endif
#if BINARY
			printf ("# right %d'th eigenvector[%d]=%g (binary ", 
			            i,j, rev[j+i*dim]);
			r = rev[j+i*dim] / (rev[25+i*dim] * (1<<25));
			r *= tn;
			prtbin (r);
			printf (" )\n");
#endif 
			printf ("# right %d'th eigenvector[%d]=%g \n", 
			            i,j, rev[j+i*dim]);
			            
			tn *= 2.0;
			thrn *= sqrt (3.0);
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
