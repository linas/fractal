/*
 * ache.c
 *
 * routines that compute various commonly required 
 * things, such as matrix elements for the continued fraction
 *
 * Linas  Vepstas
 * December 2003
 */

#include <math.h>
#include <stdio.h>

#include "ache.h"
#include "zetafn.h"

// return the matrix element for H_mp
//
long double
ache_mp(int m, int p)
{
	int k;

	long double acc = 0.0L;
	long double sign = 1.0L;
	for (k=0; k<=p; k++)
	{
		long double term = zetam1 (k+m+2);
		term *= binomial (m+k+1,m);
		term *= binomial (p,k);
		term *= sign;
		acc += term;
		sign = -sign;
	}
	return acc;
}

// ==========================================================
// actually, this one fails to improve by much (or by any?)
// over the so-called "a_sub_n_poor_convergence" below.
// oh well
long double 
a_sub_n (int n)
{
	int k;
	long double val = 0.0L;
	
	// minimize roundoff errors by doing this sum first
	for (k=1; k<=n; k++)
	{
		val -= 1.0L/((long double) (k+1));
	}
	val += 1.0L - M_GAMMA;
	val -= 0.5L/((long double) (n+1));

	// the following sum is still badly behaved
	long double acc = 0.0L;
	long double sign = -1.0L;
	for (k=1; k<=n; k++)
	{
		long double term = zetam1 (k+1)/ (long double) (k+1);
		term *= binomial (n,k);
		term *= sign;
		acc += term;
		// printf ("duuude a_sub_n k=%d term=%Lg, acc=%Lg\n", k, term, acc);
		sign = -sign;
	}
	// printf ("finally asub_n=%Lg+%Lg\n",val, -acc);
	return val-acc;
}


long double 
a_sub_n_poor_convergence (int n)
{
	long double val = 0.0L;
	int k;
	long double sk = -1.0L;
	long double acc = 1.0L - M_GAMMA;
	for (k=1; k<=n; k++)
	{
		long double term = 1.0L/ (long double) k;
		term -= zetam1(k+1);
		term /= (long double) (k+1);
		// term -= zetam1(k+1)/((long double) (k+1));
		term *= binomial (n,k);
		term *= sk;
		acc += term;
		// printf ("duuude a_n k=%d term=%Lg, acc=%Lg\n", k, term, acc);
		sk = -sk;
	}
	val = acc - 0.5/((long double) (n+1));
	return val;
}

// ==========================================================
// return the p'th element of the zeroth eigenvector.
// this is the eigenvector with eigenvalue 1
// Note that zer_m = sum_p h_mp zer_p and that
// zer is normalized to unit length.  
// zer is just the taylor coeffs of 1/(2-y)
long double zer_p (int p)
{
	int k;
	long double term = 0.5 * sqrt (3.0L);
	for (k=0; k<p; k++) term *= 0.5L;
	return term;
}

// ==========================================================
// t_sub_ne for general expanstion point "a"
// The a_sub_n is this, for a=1, and leading fator subtracted.
//
long double 
t_sub_n (int n, long double a)
{
	int k;
	long double val = 0.0L;
	
	// minimize roundoff errors by doing this sum first
	val += 1.0L - M_GAMMA;

	// the following sum is still badly behaved
	long double acc = 0.0L;
	long double an = -a;
	for (k=1; k<=n; k++)
	{
		long double term = 1.0L/ (long double) (k);
		term -= zetam1 (k+1)/ (long double) (k+1);
		term *= binomial (n,k);
		term *= an;
		acc += term;
		// printf ("duuude a_sub_n k=%d term=%Lg, acc=%Lg\n", k, term, acc);
		an *= -a;
	}
	// printf ("finally asub_n=%Lg+%Lg\n",val, -acc);
	return val+acc;
}

