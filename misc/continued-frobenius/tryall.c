
/* tryall.c
 *
 * Try as many guesses at the eigenvectors as possible.
 * try the combinatorial explosion of basic combinations
 * of vectors.
 *
 * failed tries:
 * range = 0 to 10, NLIN = 3, no prods
 * range = -9 to 9, NLIN = 2, 2x prods
 * 
 *
 * Linas Dec 2003
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ache.h"
#include "zetafn.h"


#define MS 30 // vector length

typedef double dubya;
// typedef long double dubya;

typedef dubya vector[MS];
typedef dubya matrix[MS][MS];

// ============================================================
// initialize a static copy of the H matrix
matrix *
init_ache (void)
{
	matrix * mat = (matrix *) malloc (sizeof (matrix));

	int m,p;
	for (m=0; m<MS; m++)
	{
		for (p=0; p<MS; p++)
		{
			(*mat)[m][p] = ache_mp(m,p);
		}
	}
	return mat;
}

vector zeroth_eigenvec;

void
init_zeroth (void)
{
	int k;
	dubya acc = sqrtl(3.0L) / 2.0L;
	for (k=0; k<MS; k++)
	{
		zeroth_eigenvec[k] = acc;
		acc *= 0.5;
	}
}

// ============================================================
// misc utility funcs

long double sum_over_k (int n)
{
	int k;
	long double acc = 0.0L;
	for (k=1; k<n; k++)
	{
		acc += 1.0L / ((long double) k);
	}
	return acc;
}

long double sum_over_k_over_two (int n)
{
	int k;
	long double acc = 0.0L;
	long double tp = 1.0L;
	for (k=1; k<n; k++)
	{
		acc += 1.0L / ((long double) k);
		tp *= 0.5L;
	}
	acc *= tp;
	return acc;
}

long double zeta_binomial_div_k (int n)
{
	int k;
	long double sign = -1.0L;
	long double acc = 0.0L;
	for (k=1; k<n; k++)
	{
		long double term = zetam1 (k+1);
		term /= ((long double) k+1);
		term *= binomial (n,k);
		term *= sign;
		acc += term;
		sign = -sign;
	}
	return acc;
}

long double zeta_binomial_alt (int n)
{
	int k;
	long double sign = -1.0L;
	long double acc = 0.0L;
	for (k=1; k<n; k++)
	{
		long double term = zetam1 (k+1);
		term *= binomial (n,k);
		term *= sign;
		acc += term;
		sign = -sign;
	}
	return acc;
}

long double zeta_binomial (int n)
{
	int k;
	long double acc = 0.0L;
	for (k=1; k<n; k++)
	{
		long double term = zetam1 (k+1);
		term *= binomial (n,k);
		acc += term;
	}
	return acc;
}

long double zeta_times_k_alt (int n)
{
	int k;
	long double sign = -1.0L;
	long double acc = 0.0L;
	for (k=n; k<n+60; k++)
	{
		long double term = zetam1 (k+2);
		term *= (long double) (k+1);
		term *= sign;
		acc += term;
		sign = -sign;
	}
	return acc;
}

long double zeta_times_k_alt_p (int n)
{
	int k;
	long double acc = zeta_times_k_alt (n);
	acc += (2.0L/3.0L) * ((long double) (n+1)) * zetam1(n+2);
	return acc;
}

// ============================================================
// set up a set of basis vectors from which guesses will be 
// combined.
//

#define NBASIS 80
int do_combo[NBASIS];

#define VINIT(expr) {  \
	printf ("\nInitial basis vector %d: %s:\n", ib, #expr);  \
	int k; \
	dubya sign = 1.0L;  \
	for (k=0; k<MS; k++)  \
	{  \
		basis[ib][k] = expr; \
		printf ("%d %d %g\n", ib,k, basis[ib][k]);   \
		sign = - sign;  \
	}  \
	ib ++; \
	if (NBASIS <= ib) { printf ("NOT ENOUGH!\n"); exit (1); } \
}

vector *
init_basis_vectors (int *len)
{
	vector * basis = (vector *) malloc (NBASIS * sizeof (vector));

	int ib = 0;
	
#ifdef FULL_SET
	VINIT (sqrtl((long double)(k+4)));
	VINIT (zetam1 (k+4));
	VINIT (sqrtl((long double)(k+3)));
	VINIT (zetam1 (k+3));
	VINIT (sqrtl((long double)(k+2)));
	VINIT (zetam1 (k+2));
	VINIT (sqrtl((long double)(k+1)));
	VINIT (1.0L/ (long double)(k+1));
	VINIT (1.0L/ (long double) (k+2));
	VINIT (1.0L/ powl (2.0L, k));
	VINIT (sign);
	VINIT (1.0L);
	VINIT (M_GAMMA);
	VINIT (1.0L/sqrtl(2.0L));
	VINIT (1.0L/sqrtl(3.0L));
	VINIT (1.0L/sqrtl(5.0L));
	VINIT (sum_over_k(k));
	VINIT (sum_over_k(k+1));
	VINIT (sum_over_k(k+2));
	VINIT (zeta_binomial_div_k(k));
	VINIT (zeta_binomial_div_k(k+1));
	VINIT (zeta_binomial_div_k(k+2));
	VINIT (zeta_binomial(k));
	VINIT (zeta_binomial(k+1));
	VINIT (zeta_binomial(k+2));
	VINIT (zeta_binomial_alt(k));
	VINIT (zeta_binomial_alt(k+1));
	VINIT (zeta_binomial_alt(k+2));
	VINIT (a_sub_n(k));
	VINIT (a_sub_n(k+1));
	VINIT (a_sub_n(k+2));
	VINIT (k);
	VINIT (k+1);
	VINIT (k+2);
	VINIT (1.0L/ (long double) (k+3));
	VINIT (1.0L/ powl (3.0, k));
	VINIT (factorial (k+1));
	VINIT (1.0L/factorial (k+1));
	VINIT (1.0L/sqrtl((long double)(k+1)));
	VINIT (k+3);
	VINIT (1.0L/(long double) (k+4));
	VINIT (binomial (k+2,k));
	VINIT (binomial (k+3,k));
	VINIT (binomial (k+4,k));
	VINIT (1.0L / binomial (k+2,k));
	VINIT (1.0L / binomial (k+3,k));
	VINIT (1.0L / binomial (k+4,k));
	VINIT ((k>0) ? 0.0L : 1.0L);
	VINIT ((k>0) ? 1.0L : 0.0L);
	VINIT ((k>1) ? 1.0L : 0.0L);
	VINIT ((k>2) ? 1.0L : 0.0L);
	VINIT ((k>3) ? 1.0L : 0.0L);
	VINIT (zetam1 (k+5));
	VINIT (zetam1 (k+1));
	VINIT (zetam1 (k));
	VINIT (zetam1 (k-1));
	VINIT (powl (2.0L, k));
	VINIT (powl (3.0L, k));
	VINIT (1.0L/ (((long double)k)-2.5L));
	VINIT (1.0L/ (((long double)k)-1.5L));
	VINIT (1.0L/ (((long double)k)-0.5L));
	VINIT (1.0L/ (((long double)k)+0.5L));
	VINIT (1.0L/ (((long double)k)+1.5L));
	VINIT (1.0L/ (((long double)k)+2.5L));
	VINIT (fbinomial(((long double) k)-1.5L,k));
	VINIT (fbinomial(((long double) k)-0.5L,k));
	VINIT (fbinomial(((long double) k)+0.5L,k));
	VINIT (fbinomial(((long double) k)+1.5L,k));
	VINIT (fbinomial(((long double) k)+2.5L,k));
	VINIT (fbinomial(((long double) k)+3.5L,k));
#endif

#define TRIM_SET
#ifdef TRIM_SET
	do_combo[ib] = 1; VINIT (1.0L);
	do_combo[ib] = 0; VINIT (zeta_times_k_alt_p (k));
	do_combo[ib] = 0; VINIT (zeta_times_k_alt_p (k+1));
	do_combo[ib] = 0; VINIT (zeta_times_k_alt_p (k+2));
	do_combo[ib] = 0; VINIT (zeta_times_k_alt (k));
	do_combo[ib] = 0; VINIT (zeta_times_k_alt (k+1));
	do_combo[ib] = 0; VINIT (zeta_times_k_alt (k+2));
	do_combo[ib] = 1; VINIT (k+1);
	do_combo[ib] = 1; VINIT (k+2);
	do_combo[ib] = 1; VINIT (k+3);
	do_combo[ib] = 1; VINIT (k+4);
	do_combo[ib] = 1; VINIT (zetam1 (k+2));
	do_combo[ib] = 0; VINIT (zetam1 (k+3));
	do_combo[ib] = 0; VINIT (zetam1 (k+4));
	do_combo[ib] = 0; VINIT (zetam1 (k+5));
	do_combo[ib] = 0; VINIT (sum_over_k_over_two(k+2));
	do_combo[ib] = 1; VINIT (M_GAMMA);
	do_combo[ib] = 1; VINIT (1.0L/sqrtl(2.0L));
	do_combo[ib] = 1; VINIT (sqrtl((long double)(k+2)));
	do_combo[ib] = 1; VINIT (sqrtl((long double)(k+1)));
	do_combo[ib] = 1; VINIT (1.0L/ (long double)(k+1));
	do_combo[ib] = 1; VINIT (1.0L/ (long double) (k+2));
	do_combo[ib] = 1; VINIT (1.0L/ powl (2.0L, k));
	do_combo[ib] = 0; VINIT (sum_over_k(k+1));
	do_combo[ib] = 0; VINIT (sum_over_k(k+2));
	do_combo[ib] = 0; VINIT (sum_over_k(k+3));
	do_combo[ib] = 0; VINIT (zeta_binomial_div_k(k));
	do_combo[ib] = 0; VINIT (zeta_binomial_div_k(k+1));
	do_combo[ib] = 0; VINIT (zeta_binomial_div_k(k+2));
	do_combo[ib] = 0; VINIT (zeta_binomial(k));
	do_combo[ib] = 0; VINIT (zeta_binomial(k+1));
	do_combo[ib] = 0; VINIT (zeta_binomial(k+2));
	do_combo[ib] = 0; VINIT (zeta_binomial_alt(k));
	do_combo[ib] = 0; VINIT (zeta_binomial_alt(k+1));
	do_combo[ib] = 0; VINIT (zeta_binomial_alt(k+2));
	do_combo[ib] = 0; VINIT (a_sub_n(k));
	do_combo[ib] = 0; VINIT (a_sub_n(k+1));
	do_combo[ib] = 0; VINIT (a_sub_n(k+2));
	do_combo[ib] = 1; VINIT (factorial (k+1));
	do_combo[ib] = 1; VINIT (1.0L/factorial (k+1));
	do_combo[ib] = 1; VINIT (1.0L/factorial (k+2));
	do_combo[ib] = 1; VINIT (1.0L/sqrtl((long double)(k+1)));
	do_combo[ib] = 0; VINIT (binomial (k+2,k));
	do_combo[ib] = 0; VINIT (binomial (k+3,k));
	do_combo[ib] = 0; VINIT ((k>0) ? 0.0L : 1.0L);
	do_combo[ib] = 0; VINIT ((k>0) ? 1.0L : 0.0L);
	do_combo[ib] = 0; VINIT ((k>1) ? 1.0L : 0.0L);
	do_combo[ib] = 0; VINIT ((k>2) ? 1.0L : 0.0L);
	do_combo[ib] = 1; VINIT (1.0L/ (((long double)k)+0.5L));
	do_combo[ib] = 1; VINIT (1.0L/ (((long double)k)+1.5L));
	do_combo[ib] = 0; VINIT (fbinomial(((long double) k)-0.5L,k));
	do_combo[ib] = 0; VINIT (fbinomial(((long double) k)+0.5L,k));
	do_combo[ib] = 0; VINIT (fbinomial(((long double) k)+1.5L,k));
	do_combo[ib] = 0; VINIT (fbinomial(((long double) k)+2.5L,k));
#endif
	fflush (stdout);
	
	*len = ib;
	
	return basis;
}

// ============================================================
// create multiplicative combinations of basis vectors

vector *
init_double_prod (vector *basis, int baselen, int *retlen)
{
	*retlen = baselen * (baselen+1);
	*retlen /= 2;
	vector * prod = (vector *) malloc ((*retlen) * sizeof (vector));

	printf ("double combo \n");
	int jdx =0;
	int i,j;
	for (i=0; i<baselen; i++)
	{
		vector *ba = &basis[i];
		for (j=i; j<baselen; j++)
		{
			// we don't want to pair up certain pairs
			if (0 == do_combo[i]+do_combo[j]) continue;

			vector *bb = &basis[j];
			vector *v = &prod[jdx];
			int k;
			for (k=0; k<MS; k++)
			{
				(*v)[k] = ((*ba)[k]) * ((*bb)[k]);
			}
			printf ("product [%d] = [%d] * [%d]\n", jdx, i, j);
			jdx ++;
		}
	}
	*retlen = jdx;
	
	return prod;
}

vector *
init_triple_prod (vector *basis, int baselen, int *retlen)
{
	*retlen = baselen * (baselen+1)*(baselen+2);
	*retlen /= 6;
	vector * prod = (vector *) malloc ((*retlen) * sizeof (vector));

	printf ("triple combo \n");
	int jdx =0;
	int i,j,m;
	for (i=0; i<baselen; i++)
	{
		vector *ba = &basis[i];
		for (j=i; j<baselen; j++)
		{
			vector *bb = &basis[j];
			for (m=j; m<baselen; m++)
			{
				vector *bc = &basis[m];
				vector *v = &prod[jdx];
				int k;
				for (k=0; k<MS; k++)
				{
					(*v)[k] = ((*ba)[k]) * ((*bb)[k]) * ((*bc)[k]);
				}
				printf ("product [%d] = [%d] * [%d] * [%d]\n", jdx, i, j, m);
				jdx ++;
			}
		}
	}
	if (jdx != (*retlen))
	{
		printf ("Ooooops!! %d %d \n", jdx, *retlen);
		exit (1);
	}
	
	return prod;
}

// ============================================================
// make bais vectors orthogonal to zeroth eigenvector

void
orthogonalize (vector *basis, int baselen)
{
	printf ("orthogonalizing the basis vectors !\n");
	int ib;
	for (ib=0; ib<baselen; ib++)
	{
		vector *w = &basis[ib];
		// orthogonalize to eigenvec
		dubya dot = 0.0;
		int n;
		for (n=0; n<MS; n++)
		{
			dot += zeroth_eigenvec[n] * (*w)[n];
		}
		for (n=0; n<MS; n++)
		{
			(*w)[n] -= zeroth_eigenvec[n] * dot;
		}
	}
}

// ============================================================
// evaluate a trial vector for eignehood

static unsigned long long int nattempts = 0;
static unsigned long long int ntoobig = 0;
static unsigned long long int ntrivial = 0;
static unsigned long long int ndive = 0;
static unsigned long long int nlate3 = 0;
static unsigned long long int nlate4 = 0;
static unsigned long long int nlate5 = 0;
static unsigned long long int nlate6 = 0;
static unsigned long long int nlate7 = 0;
static unsigned long long int nlate8 = 0;
static unsigned long long int nlater = 0;
static unsigned long long int nqnan = 0;
static unsigned long long int nrediscover = 0;
static unsigned long long int nclosecall = 0;
static unsigned long long int ngood = 0;

int evaluate (matrix *H, vector *vec, dubya goodness, int nterms)
{
	nattempts ++;
	if (0 == nattempts%4000100)
	{
		double pi = 100.0 *(double) ndive / (double) nattempts;
		double pt = 100.0 *(double) ntoobig / (double) nattempts;
		double pr = 100.0 * (double) ntrivial / (double) nattempts;
		double pl = 100.0 * (double) (nlate3+nlate4+nlate5+nlate6+nlate7+nlate8+nlater) / (double) nattempts;
		double pd = 100.0 * (double) nrediscover / (double) nattempts;
		double pc = 100.0 * (double) nclosecall / (double) nattempts;
		double pq = 100.0 * (double) nqnan / (double) nattempts;
		printf ("n=%lld %5.2f %5.2f %5.2f l3=%lld %lld %lld l6=%lld %lld %lld %lld %5.2f re=%lld %5.2f cl=%lld %5.2f nan=%lld %5.2f good=%lld\n", 
				 nattempts, pi, pt, pr, nlate3, nlate4, nlate5, nlate6, nlate7, nlate8, nlater, pl, nrediscover, pd, nclosecall, pc, nqnan, pq, ngood);
		fflush (stdout);
	}

	int k = MS-1;
	if ((1.0L < (*vec)[k]) || (-1.0L > (*vec)[k]))
	{
		ndive ++;
		return 99;
	}
	
#define JORDANIZE
#ifdef JORDANIZE
	vector w;
	int n;
	for (n=0; n<MS; n++)
	{
		dubya acc = 0.0L;
		for (k=0; k<MS; k++)
		{
			acc += (*H)[n][k] * (*vec)[k]; 
		}
		w[n] = acc;
	}

	// now orthogonalize to eigenvec
	dubya dot = 0.0;
	for (n=0; n<MS; n++)
	{
		dot += zeroth_eigenvec[n] * w[n];
	}
	for (n=0; n<MS; n++)
	{
		w[n] -= zeroth_eigenvec[n] * dot;
	}
	
	dubya w0 = w[0] / (*vec)[0];

	// resume simple tests
	if (isnan (w0))
	{
		nqnan ++;
		return 1;
	}

	// eigenvalue must be less than 1
	// (eventually, let first guesses go high)
	if ((1.3L <= w0) || (-2.0L >= w0))
	{
		ntoobig ++;
		return 2;
	}

	// sample 'widely'
	double step = ((double) MS-2) / ((double) nterms-2);
	int nstep = (int) step;

	// eigenvalue must be eigen-like
	dubya w1 = w[1] /(*vec)[1];

	dubya wd = w0-w1;
#define SCALE_GOOD
#ifdef SCALE_GOOD
	goodness *= w0;
	if (0.0L > goodness) goodness = -goodness;
#endif
	if ((goodness <= wd) || (-goodness >= wd))
	{
		ntrivial ++;
		return 3;
	}

	// try higher up the eigenvector
	int ntest = 3;
	for (n=2; n<MS; n+= nstep) 
	{
		dubya wn = w[n]/ (*vec)[n];
		dubya wd = w0-wn;
		if ((goodness <= wd) || (-goodness >= wd))
		{
			if (3 == ntest) { nlate3 ++; return 4; }
			if (4 == ntest) { nlate4 ++; return 5; }
			if (5 == ntest) { nlate5 ++; return 6; }
			if (6 == ntest) { nlate6 ++; return 7; }
			if (7 == ntest) { nlate7 ++; return 8; }
			if (8 == ntest) { nlate8 ++; return 9; }
			nlater ++;
			return 10;
		}
		ntest ++;
	}
	
	// did we rediscover the zeroth eigenval ??
	if (fabs (w0-1.0) < 10.e-9) 
	{
		// printf ("Rediscover !\n");
		nrediscover ++;
		return 11;
	}
	
	if (fabs (w0-1.0) < 0.001) 
	{
		// printf ("Rediscover !\n");
		nclosecall ++;
		return 12;
	}
	
	// if (8 < nterms)
	{
		// if we got to here, we have what smells like
		// an eignevector!
		dubya avg = 0.0L;
		for (n=0; n<MS; n++)
		{
			dubya term = w[n] / (*vec)[n];
			avg += term;
		}
		avg /= ((dubya) MS);
		
		dubya dev = 0.0L;
		for (n=0; n<MS; n++)
		{
			dubya term = w[n] / (*vec)[n] - avg;
			dev += term*term;
		}
		dev /= (dubya) MS;
		dev = sqrt(dev);
		
		printf ("===================================\n");
		printf ("CANDidate eigen=%g    DEVIATION=%g\n", avg, dev);

		printf ("candidate eigenvalue w[0]=%g   w[1]=%g   ", w0, w1);
		for (n=2; n<MS; n += nstep) 
		{
			double wn = w[n] / (*vec)[n];
			printf ("w[%d]=%g   ",n, (double)wn);
		}
		printf ("\n");
		for (n=0; n<MS; n++)
		{
			dubya eval = w[n] / (*vec)[n];
			printf ("cand eigenvec v[%d]=%g Hv=%g eigval=%g diff=%g\n", n, (*vec)[n], w[n], eval, w0-eval); 
		}
		fflush (stdout);
		ngood ++;
	}
#endif 


#if PLAIN_EIGENVEC
	dubya w0 = 0.0L;
	for (k=0; k<MS; k++)
	{
		w0 += (*H)[0][k] * (*vec)[k]; 
	}
	w0 /= (*vec)[0];

	if (isnan (w0))
	{
		nqnan ++;
		return 1;
	}

	// eignvalue must be less than 1
	// (eventually, let first guesses go high)
	if ((0.7L <= w0) || (-2.0L >= w0))
	{
		ntoobig ++;
		return 2;
	}

	// sample 'widely'
	double step = ((double) MS-2) / ((double) nterms-2);
	int nstep = (int) step;

	dubya w1 = 0.0L;
	for (k=0; k<MS; k++)
	{
		w1 += (*H)[1][k] * (*vec)[k]; 
	}
	
	// eigenvalue must be eigen-like
	w1 /= (*vec)[1];
	dubya wd = w0-w1;
#if SCALE_GOOD
	goodness *= w0;
	if (0.0L > goodness) goodness = -goodness;
#endif
	if ((goodness <= wd) || (-goodness >= wd))
	{
		ntrivial ++;
		return 3;
	}

	// try higher up the eigenvector
	int n;
	int ntest = 3;
	for (n=2; n<MS; n+= nstep) 
	{
		dubya w = 0.0L;
		for (k=0; k<MS; k++)
		{
			w += (*H)[n][k] * (*vec)[k]; 
		}
		w /= (*vec)[n];
		dubya wd = w0-w;
		if ((goodness <= wd) || (-goodness >= wd))
		{
			if (3 == ntest) { nlate3 ++; return 4; }
			if (4 == ntest) { nlate4 ++; return 5; }
			if (5 == ntest) { nlate5 ++; return 6; }
			if (6 == ntest) { nlate6 ++; return 7; }
			if (7 == ntest) { nlate7 ++; return 8; }
			if (8 == ntest) { nlate8 ++; return 9; }
			nlater ++;
			return 10;
		}
		ntest ++;
	}
	
	// did we rediscover the zeroth eigenval ??
	if (fabs (w0-1.0) < 10.e-9) 
	{
		// printf ("Rediscover !\n");
		nrediscover ++;
		return 11;
	}
	
	if (fabs (w0-1.0) < 0.001) 
	{
		// printf ("Rediscover !\n");
		nclosecall ++;
		return 12;
	}
	
	// if (8 < nterms)
	{
		// if we got to here, we have what smells like
		// an eignevector!
		printf ("plain candidate eigenvalue w[0]=%g   w[1]=%g   ", w0, w1);
		for (n=2; n<MS; n += nstep) 
		{
			dubya wn = 0.0L;
			for (k=0; k<MS; k++)
			{
				wn += (*H)[n][k] * (*vec)[k]; 
			}
			wn /= (*vec)[n];
			printf ("w[%d]=%g   ",n, (double)wn);
		}
		printf ("\n");
		for (n=0; n<MS; n++)
		{
			dubya wn = 0.0L;
			for (k=0; k<MS; k++)
			{
				wn += (*H)[n][k] * (*vec)[k]; 
			}
			dubya eval = wn / (*vec)[n];
			printf ("plain cand eigenvec v[%d]=%g Hv=%g eigval=%g diff=%g\n", n, (*vec)[n], wn, eval, w0-eval); 
		}
		fflush (stdout);
		ngood ++;
	}
#endif

	return 0;
}
					 
// ============================================================
// try the various different combinatorial possibilities

static vector *combo_basis = NULL;
static int combo_basis_len = 0;

#define NLIN 20
int num[NLIN];
int deno[NLIN];
dubya alt[NLIN];

int idx[NLIN];

#define MAXCOFF 9

void 
init_combo_basis (vector *basis, int nvec)
{
	combo_basis = basis;
	combo_basis_len = nvec;
}

int idx_end[NLIN];

void 
init_combo (int nlow)
{
	int i;
	for (i=0; i<NLIN; i++)
	{
		num[i] = -MAXCOFF;
		deno[i] = 1;
		alt[i] = 1.0L;
		idx_end[i] = 99999999;
	}
	for (i=nlow; i<NLIN; i++)
	{
		idx[i] = i-nlow;
	}
	printf ("init combo ");
	for (i=0; i<nlow+1; i++)
	{
		printf ("%d ", idx[i]);
	}
	printf (" ========\n");
	fflush (stdout);
}


vector * 
next_combo (int nlin, int nlow)
{
	static vector combo_vec;
	dubya frac[NLIN];
	dubya sign[NLIN];

	int i;
	for (i=nlin-1; i>=0; i--)
	{
		alt[i] = -alt[i];
		if (0.0 > alt[i]) goto nextvec;
		
		// only explore relatively-prime combos
		// int reject=0;
another:
		// if (reject) {printf ("reject %d/%d   %d/%d   %d/%d\n", num[0], deno[0], num[1], deno[1], num[2], deno[2]); fflush (stdout); }
		// reject = 1;
		num[i] += 1;
		if (num[i] > MAXCOFF) goto nextdeno;
		if ((1 != deno[i]) && (0 == (num[i]%deno[i]))) goto another;
		if (((1<num[i]) || (-1>num[i])) && (0 == (deno[i]%num[i]))) goto another;
		if ((0 == (num[i]%2)) && (0 == deno[i]%2)) goto another;
		if ((0 == (num[i]%3)) && (0 == deno[i]%3)) goto another;
		if ((2==nlin) && (0==deno[0]%2) && (0==deno[1]%2)) goto another;
		if ((2==nlin) && (0==deno[0]%3) && (0==deno[1]%3)) goto another;
		if ((3==nlin) && (0==deno[0]%2) && (0==deno[1]%2) && (0==deno[2]%2)) goto another;
		if ((3==nlin) && (0==deno[0]%3) && (0==deno[1]%3) && (0==deno[2]%3)) goto another;
		if ((4==nlin) && (0==deno[0]%2) && (0==deno[1]%2) && (0==deno[2]%2) &&(0==deno[3]%2)) goto another;
		if ((4==nlin) && (0==deno[0]%3) && (0==deno[1]%3) && (0==deno[2]%3) &&(0==deno[3]%3)) goto another;
		
		if (num[i] <= MAXCOFF) goto nextvec;

nextdeno:
		num[i] = -MAXCOFF;
		if (0 == i) num[i] = 1;  // only need to search positive side of first one
		deno[i] += 1;
		if (deno[i] <= MAXCOFF) goto nextvec;

		deno[i] = 1;
	}

	for (i=nlin-1; i>=nlow; i--)
	{
		idx[i] ++;
		if (idx[i] < combo_basis_len-(nlin-i-1)) 
		{
			if (idx[i] >= idx_end[i]) 
			{
				printf ("ending, because next combo too far ");
				int j;
				for (j=0; j<nlin; j++) printf ("%d ", idx[j]);
				printf (" =============================\n");
				fflush (stdout);
				return NULL;
			}
			if (0 == (idx[nlin-1] %120))
			{
				printf ("try combo ");
				int j;
				for (j=0; j<nlin; j++) printf ("%d ", idx[j]);
				printf (" =============================\n");
				fflush (stdout);
			}
			goto nextvec;
		}
		if (0 == i) return NULL;
		idx[i] = idx[i-1] +2;
		int j;
		for (j=i+1; j<nlin; j++)
		{
			idx[j] = idx[j-1] +1;
		}
	}

	return NULL;
	
nextvec:

	// premultiply -- optimization
	for (i=0; i<nlin; i++)
	{
  		frac[i] = ((dubya) num[i]) / ((dubya) deno[i]);
		sign[i] = 1.0L;
	}
	
	// don't compute all the terms if the highest term is too big.
	// this should save some cpu time.
	int k = MS-1;
	combo_vec[k] = 0.0L;
	for (i=0; i<nlin; i++)
	{
		dubya s = 1.0L;
		if (k%2 == 1) s *= alt[i];
		dubya term = s * frac[i] * combo_basis[idx[i]][k];
		combo_vec[k] += term;
	}
	if ((1.0L > combo_vec[k]) && (-1.0L < combo_vec[k]))
	{ 
		for (k=0; k<MS-1; k++)
		{
			combo_vec[k] = 0.0L;
			for (i=0; i<nlin; i++)
			{
				dubya term = sign[i] * frac[i] * combo_basis[idx[i]][k];
				combo_vec[k] += term;
				sign[i] = alt[i] * sign[i];
			}
		}
	}

	return &combo_vec;
}

void
set_combo_state (int i0, int i1, int i2, int i3)
{
	idx[0] = i0;
	idx[1] = i1;
	idx[2] = i2;
	idx[3] = i3;
	
	int nlin = 4;
	printf ("setup combo ");
	int j;
	for (j=0; j<nlin; j++) printf ("%d ", idx[j]);
	printf (" =============================\n");
	fflush (stdout);
}

void
set_end_state (int i0, int i1, int i2, int i3)
{
	idx_end[0] = i0;
	idx_end[1] = i1;
	idx_end[2] = i2;
	idx_end[3] = i3;
	
	int nlin = 4;
	printf ("setup end state: ");
	int j;
	for (j=0; j<nlin; j++) printf ("%d ", idx[j]);
	printf (" =============================\n");
	fflush (stdout);
}

void 
dump_combo_state (int nlin)
{
	printf ("Current combinatorics state:\n");
	int i;
	
	for (i=0; i<nlin; i++)
	{
		printf ("%d num=%d deno=%d   alt=%d  idx=%d\n", i, num[i], deno[i], (int) alt[i], idx[i]);
	}
	fflush (stdout);
}

// ============================================================
// run the search algo
// 
void 
run (void)
{
	matrix *h = init_ache ();
	init_zeroth ();

	int nbasis = 0;
	vector *basis = init_basis_vectors (&nbasis);
	
#if 0
	// tried these, didn't get anything 
	int nc = nbasis;
	vector *cv = basis;
#else
	
	int nprod = 0;
	vector *prod = init_double_prod (basis, nbasis, &nprod);
	// vector *prod = init_triple_prod (basis, nbasis, &nprod);
	
	int nc = nprod;
	vector *cv = prod;
#endif
	init_combo_basis (cv, nc);
	orthogonalize (cv, nc);

	init_combo (0);

	
#if DO_DEPTH_FIRST_PROMOTION_SEARCH
	printf ("Progressive attack \n");
	int nlin = 2;
	int nterms = 5;
	dubya accuracy = 0.3L;
	
	vector *vec = next_combo(nlin, 0);
	while (1)
	{
		int rc = evaluate (h, vec, accuracy, nterms);
		if (0 == rc)
		{
			accuracy *= 0.25L;
			nterms += 3;
			nlin ++;
			if (8 < nterms) 
			{
				printf ("YESSSS !!!!!!   got one !!! !!!!!!!!!!!!!!!!!!!!!!!!!\n");
				dump_combo_state (nlin);
			}
			printf ("escalate %d linear combo accurate to %d terms at %g accuracy\n", nlin, nterms, (double) accuracy);
			init_combo (nlin-1);
		}
		vec = next_combo(nlin, nlin-1);
		while (!vec)
		{
			nlin --;
			nterms -= 3;
			accuracy *= 4.0L;
			if (1 > nlin) goto done;
			printf ("fallback %d linear combo accurate to %d terms at %g accuracy\n", nlin, nterms, (double) accuracy);
			vec = next_combo(nlin, nlin-1);
		}
	}
done:
	return;
#endif

	// try a breadth-first search instead.
	int nlin = 2;
	int nterms = 9;
	dubya accuracy = 0.002L;

	printf ("linear combo of %d demand %d terms match to %f accuracy\n", 
						 nlin, nterms, accuracy);
	
	// set_combo_state (4,48,67,86);
	vector *vec = next_combo(nlin, 0);
	while (vec)
	{
		int rc = evaluate (h, vec, accuracy, nterms);
		if (0 == rc)
		{
			printf ("YESSSS !!!!!!   got one !!! !!!!!!!!!!!!!!!!!!!!!!!!!\n");
			dump_combo_state (nlin);
		}
		vec = next_combo(nlin, 0);
	}
}

// ============================================================

int
main (int argc, char *argv[])
{
	double p;
	run();

	printf ("Number of Attempts = %lld\n", nattempts);
	
	p = 100.0 *(double) ndive / (double) nattempts;
	printf ("Number that were divergent = %lld (%g %%)\n", ndive, p);

	p = 100.0 *(double) ntoobig / (double) nattempts;
	printf ("Number that were too big = %lld (%g %%)\n", ntoobig, p);

	p = 100.0 * (double) ntrivial / (double) nattempts;
	printf ("Number that fail at first attempt = %lld (%g %%)\n", ntrivial, p);
	
	p = 100.0 * (double) nlater / (double) nattempts;
	printf ("Number that fail at later stages = %lld (%g %%)\n", nlater, p);
	
	p = 100.0 * (double) nrediscover / (double) nattempts;
	printf ("Number that were rediscovered = %lld (%g %%)\n", nrediscover, p);

	p = 100.0 * (double) nclosecall / (double) nattempts;
	printf ("Number that were near the zeroth eigenvec = %lld (%g %%)\n", nclosecall, p);

	p = 100.0 * (double) nqnan / (double) nattempts;
	printf ("Number that are NaN = %lld (%g %%)\n", nqnan, p);

	long long int nlost = nattempts;
	nlost -= ntoobig+ndive+ntrivial+nlate3+nlate4+nlate5+nlate6+nlate7+nlate8+nlater+nrediscover+nclosecall+nqnan+ngood;
	if (0 != nlost)
	{
		printf ("bad total!!!! =%lld\n", nlost);
	}
	return 0;
}
