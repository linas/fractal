
/*
 * diagon.c
 *
 * Try to work out the eigenvalues and eigenvectors for 
 * the (regulated) frobenius-perron operator of the 
 * continued fraction map.
 *
 * Do this by concatenating diagonalization of 2x2 sub-pieces.
 *
 * Linas Dec 2003
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "zetafn.h"

// =================================================================

#define MS 100 // matrix dimension
#include "matrix.h"

// =================================================================

#define MT 25  // number of tee values
matrix  *gfp; // frobenius perron for continued fraction

matrix *rl; // left-rotation
matrix *rr; // right-rotation aka rl-inverse

int nrow[MT];
#define TSTEP 0.8

void initialize (void)
{
	int i,j,k;

	printf ("sizeof matrix = %d\n", sizeof (matrix));
	
	gfp = (matrix *) malloc (MT * sizeof (matrix));
	rl = (matrix *) malloc (MT * sizeof (matrix));
	rr = (matrix *) malloc (MT * sizeof (matrix));
	printf ("g=%p rl=%p rr=%p\n", gfp, rl, rr);
	
	// initialize the diagonalization matrices
	for (k=0; k<MT; k++)
	{
		matrix *mrl = &rl[k];
		matrix *mrr = &rr[k];
		
		identmatrix (mrl);
		identmatrix (mrr);
	}
	
	// initialize the frobenius-perron operator
	matrix *mg = &gfp[0];
	long double fi = 1.0;    // factorial i;
	long double sign = 1.0;  // sign
	for (i=0; i<MS; i++)
	{
		if (0 != i) fi *= (double) i;
		long double fij = fi;
		long double fjp1 = 1.0;  // factorial j+1
		for (j=0; j<MS; j++)
		{
			fjp1 *= (long double) (j+1);
			fij *= (long double) (i+j+1);
			double z = zetam1(i+j+2);
			(*mg)[i][j] = sign * fij * (1.0+z) / (fjp1 * fi);
		}
		sign = -sign;
	}

	// now regulate if so that the matrix elements are reasonable
	long double tee = 1.0;
	for (k=1; k<MT; k++)
	{
		tee *= TSTEP;
		matrix *mgt = &gfp[k];
		for (i=0; i<MS; i++)
		{
			for (j=0; j<MS; j++)
			{
				long double reg = i+j;
				// double reg = j;
				reg = expl (-reg*reg*tee*tee);
				(*mgt)[i][j] = reg * (*mg)[i][j];
			}
		}

		// max for this regulator
		nrow[k] = sqrtl (-logl (1.0e-9L) / (tee*tee));
		if (MS < nrow[k]) nrow[k] = MS;
		printf ("k=%d tee=%Lf nrow=%d \n", k, tee, nrow[k]);
	}
	nrow[0] = MS;
	
	printf ("\ng = \n");
	for (i=0; i<6; i++)
	{
		printf ("%Lf %Lf %Lf %Lf %Lf %Lf %Lf\n", gfp[0][i][0],
				gfp[0][i][1], gfp[0][i][2], gfp[0][i][3], 
				gfp[0][i][4], gfp[0][i][5], gfp[0][i][6]);
	}
#if 0
	printf ("\ng for tee=xxx \n");
	for (i=0; i<6; i++)
	{
		printf ("%f %f %f %f %f %f %f\n", gfp[1][i][0],
				gfp[1][i][1], gfp[1][i][2], gfp[1][i][3], 
				gfp[1][i][4], gfp[1][i][5], gfp[1][i][6]);
	}
#endif
	printf ("\n\n");
}

// =================================================================
// diagonalize a 2x2 subelement

int 
diag (matrix *mg, matrix *mrl, matrix *mrr, int p, int n)
{
	double g00 = (*mg)[p][p];
	double g01 = (*mg)[p][n];
	double g10 = (*mg)[n][p];
	double g11 = (*mg)[n][n];

	double disc = (g00-g11)*(g00-g11) + 4.0 * g01 * g10;
	if (0.0 >= disc)
	{
		printf ("XXX bad desc! disc=%g p=%d n=%d\n", disc, p, n);
		printf ("g00=%g g01=%g\n", g00, g01);
		printf ("g10=%g g11=%g\n", g10, g11);
		return 1;
	}
	disc = sqrt (disc);
	double lplus = 0.5 *(g00+g11 + disc);
	double lminus = 0.5 *(g00+g11 - disc);
	
	// if (fabs(1.0-lminus) < fabs(1.0-lplus))
	if (fabs(lminus) < fabs(lplus))
	{
		double tmp = lplus;
		lplus = lminus;
		lminus = tmp;
	}

	double xplus = lplus - g00;
	double xminus = lminus - g00;

	double aplus = 1.0 / sqrt (g01*g01 + xplus*xplus);
	double aminus = 1.0 / sqrt (g01*g01 + xminus*xminus);

	double vv = aplus*aminus *(g01*g01 +xplus*xminus);
	double c = 1.0 / (1.0-vv*vv);
	double s = vv *c;

	// accumulate the right rotation 
	double rr00 = aplus * g01;
	double rr01 = aminus * g01;
	double rr10 = aplus * xplus;
	double rr11 = aminus * xminus;

	matrix r, tmp;
	identmatrix (&r);
	r[p][p] = rr00;
	r[p][n] = rr01;
	r[n][p] = rr10;
	r[n][n] = rr11;

	multmatrix (&tmp, mrr, &r);
	copymatrix (mrr, &tmp);

	// accumulate the left rotation
	double rl00 = (c*aplus - s*aminus) *g01;
	double rl01 = c*aplus*xplus - s*aminus*xminus;
	double rl10 = (-s*aplus + c*aminus) * g01;
	double rl11 = -s *aplus*xplus + c*aminus*xminus;

	identmatrix (&r);
	r[p][p] = rl00;
	r[p][n] = rl01;
	r[n][p] = rl10;
	r[n][n] = rl11;

	multmatrix (&tmp, &r, mrl);
	copymatrix (mrl, &tmp);

	// diagonalize
	multmatrix (&tmp, mg, mrr);
	copymatrix (mg, &tmp);
	multmatrix (&tmp, mrl, mg);
	copymatrix (mg, &tmp);
	
#if DEBUG
	printf ("finished diag row=%d lam= %f %f eigen=%f\n", 
		 n, lplus, lminus, (*mg)[0][0]);
		// printf ("\tg00=%g \tg01=%g\n", g00, g01);
		// printf ("\tg10=%g \tg11=%g\n", g10, g11);
#endif

	return 0;
}

// =================================================================

void solve (void)
{
	int i, j, k, n, p;

#ifdef DO_NON_REGULATED
	/* do the non-regulated case */
	matrix *mg = &gfp[0];
	matrix *mrl = &rl[0];
	matrix *mrr = &rr[0];
		
	for (n=1; n<MS; n++)
	{
		diag (mg, mrl, mrr, 0, n);
		
		for (i=0; i<MS; i++)
		{
			for (j=0; j<MS; j++)
			{
				double acc = 0.0;
				for (k=0; k<MS; k++)
				{
					acc += (*mrl)[i][k] * (*mrr)[k][j];
				}
				if ((i!=j) && (fabs(acc) > 1.0e-12))
				{
					printf (" bad off-diag %d %d %g\n", i,j,acc);
				}
				else if ((i==j) &&(fabs (acc-1.0) > 1.0e-12))
				{
					printf (" bad diag %d %g\n", i,acc);
				}
										
			}
		}

		printf ("vector=%f %f %f %f %f\n", (*mrr)[0][0], (*mrr)[1][0],
				(*mrr)[2][0], (*mrr)[3][0], (*mrr)[4][0]);
		printf ("  %f %f %f %f %f\n\n", (*mrr)[5][0], (*mrr)[6][0],
				(*mrr)[7][0], (*mrr)[8][0], (*mrr)[9][0]);
	}
#endif
	
	long double tee = 1.0;
	for (k=1; k<MT; k++)
	{
		tee *= TSTEP;

		matrix *mg = &gfp[k];
		matrix *mrl = &rl[k];
		matrix *mrr = &rr[k];
		
#if DEBUG
		printf ("-----\n");
		printf ("g for tee=%f \n", tee);
		for (i=0; i<6; i++)
		{
			printf ("%f %f %f %f %f %f %f\n", (*mg)[i][0],
					(*mg)[i][1], (*mg)[i][2], (*mg)[i][3], 
					(*mg)[i][4], (*mg)[i][5], (*mg)[i][6]);
		}
#endif

		printf ("======== start k=%d tee=%g\n", k, (double) tee);
		for (p=0; p<nrow[k]; p++)
		{
			for (n=p+1; n<nrow[k]; n++)
			{
				int rc;
				// printf ("-----\n");
				printf ("start p=%d n=%d\n", p, n);
				rc = diag (mg, mrl, mrr, p, n);
				if (rc) break;
#if 0 
				printf ("v=%f %f %f %f %f\n", (*mrr)[0][0], (*mrr)[1][0],
						(*mrr)[2][0], (*mrr)[3][0], (*mrr)[4][0]);
				printf ("  %f %f %f %f %f\n\n", (*mrr)[5][0], (*mrr)[6][0],
						(*mrr)[7][0], (*mrr)[8][0], (*mrr)[9][0]);
#endif
				if (fabs(gfp[k][n][p] > 1.0e-9))
				{
					printf (" fresh bad column (%d %d) = %g\n", n,p,(double) gfp[k][n][p]);
				}
				if (fabs(gfp[k][p][n] > 1.0e-9))
				{
					printf (" fresh bad row (%d %d) = %g\n", p,n, (double) gfp[k][p][n]);
				}
			}
			
			for (n=p+1; n<nrow[k]; n++)
			{
				if (fabs(gfp[k][n][p] > 1.0e-9))
				{
					printf (" bad column (%d %d) = %g\n", n,p, (double) gfp[k][n][p]);
				}
				if (fabs(gfp[k][p][n] > 1.0e-9))
				{
					printf (" bad row (%d %d) = %g\n", p,n, (double) gfp[k][p][n]);
				}
			}
		}
		printf ("k=%d tee=%f nrow=%d first eigen=%f\n", k, tee, nrow[k], (double) gfp[k][0][0]);
		for (p=0; p<nrow[k]; p++)
		{
			printf ("\t%d	%g\n", p, (double) gfp[k][p][p]);
		}
		printf ("\n");

		// check for off-diagonals (should be none)
		for (i=0; i<MS; i++)
		{
			for (j=0; j<MS; j++)
			{
				if ((i!=j) && (fabs(gfp[k][i][j]) > 1.0e-9))
				{
					printf (" bad off-diag %d %d %g\n", i,j, (double) gfp[k][i][j]);
				}
			}
		}

	}
}

// =================================================================

void trial (void)
{
	vector vv;
	
	long double acc = 1.75L;
	acc = 1.0;
	long double pw = 2.5L;
	long double sign = 1.0L;
	int i;
	for (i=0; i<MS; i++)
	{
		vv[i]= sign * acc;
		acc *= pw;
		sign = -sign;
		pw += 1.0L;
	}
	// vv[0] = 1.0L;
	
	for (i=0; i<5; i++)
	{
		long double tee = 1.0;
		int k;
		for (k=1; k<9; k++)
		{
			tee *= TSTEP;

			matrix *mg = &gfp[k];
			int nr = nrow[k] +15;
			if (MS < nr) nr = MS;

			int j;
			long double acc = 0.0;
			for (j=0; j<MS; j++)
			{
				acc += (*mg)[i][j] * vv[j];
				printf ("tee=%Lf row=%d col=%d acc=%Lf rat=%Lf \n", 
									 tee, i,j,acc, acc/vv[i]);
			}
			printf ("================\n");
		}
	}
}

// =================================================================

main ()
{
	initialize();
	// solve();
	trial ();
}
