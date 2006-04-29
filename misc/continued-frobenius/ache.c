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
#include <stdlib.h>

#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_psi.h>

#include "ache.h"
#include "zetafn.h"

// Return the matrix element for H_mp aka the matrix element of GKW.
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
		term *= binomial (n,k);
		term *= sk;
		acc += term;
		// printf ("duuude a_n k=%d term=%Lg, acc=%Lg\n", k, term, acc);
		sk = -sk;
	}
	val = acc - 0.5/((long double) (n+1));
	return val;
}

long double 
b_sub_n (int n)
{
	if (n==0) return 0.5L;
	return ((long double) n) * a_sub_n (n-1);
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
// The a_sub_n is this, for a=1, and leading factor subtracted.
//
// That is, a_sub_n (n) == t_sub_n (n, 1.0) - 1/2(n+1)
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
		term -= zetam1 (k+1);
		term /= (long double) (k+1);
		term *= binomial (n,k);
		term *= an;
		acc += term;
		// printf ("duuude a_sub_n k=%d term=%Lg, acc=%Lg\n", k, term, acc);
		an *= -a;
	}
	// printf ("finally asub_n=%Lg+%Lg\n",val, -acc);
	return val+acc;
}

// ==========================================================
// return the harmonic numbers

long double harmonic_n (int n)
{
	int k;
	long double sum = 0.0L;
	for (k=1; k<=n; k++)
	{
		sum += 1.0L/((long double) k);
	}
	return sum;
}

long double harmonic_n2p1 (int n)
{
	return n/(n+1.0L);
}

// Return a_sub_n but for Dirichlet L-function
// 
long double 
dirichlet_L_function_a_sub_n (int n)
{
	int k;
	long double val = 0.0L;
	
	val = (2.0L*M_GAMMA) - 1.0L;
	val = M_GAMMA;
	val = 1.0 - 0.25L*M_PI;
	val *= 0.666666666;
	// printf ("duude val=%Lg\n", val);

	//
	val = 0.0L;
	// val -= M_GAMMA;
	// val -= 0.25L*M_PI;
	// val -= 0.135181;
	val -= 0.5L/((long double) (n+1));
	// val += 1.0L /((long double) 6*n*(n+1));
	// val -= 0.5*harmonic_n (n+1);

#ifdef CHARACTER_3_1
	val = 1.0L;
	val -= 0.5L/((long double) (n+1));
	val -= (2.0/3.0)*harmonic_n (n+1);
	val -= (1.0/6.0)*harmonic_n2p1 (n);
#endif

#ifdef CHARACTER_3_2
	val = 0.0L;
	val -= 0.5L/((long double) (n+1));
	val += (1.0/6.0)*harmonic_n2p1 (n);
#endif

#ifdef CHARACTER_5_1
	val = 1.0L;
	val -= 0.5L/((long double) (n+1));
	val -= (4.0/5.0)*harmonic_n (n+1);
	val -= 0.3*harmonic_n2p1 (n);
#endif

#ifdef CHARACTER_5_2_OOPS
	val = 1.0L;
	val -= 0.5L/((long double) (n+1));
	val -= (2.0/5.0)*harmonic_n (n+1);
	val += 0.1*harmonic_n2p1 (n);
#endif

#ifdef CHARACTER_5_X_1
	val = 1.0L;
	val -= 0.5L/((long double) (n+1));
	val -= (1.0/5.0)*harmonic_n (n+1);
#endif

#ifdef CHARACTER_7_X_1
	val = 1.0L;
	val -= 0.5L/((long double) (n+1));
	val -= (1.0/7.0)*harmonic_n (n+1);
#endif

#define CHARACTER_9_X_1
#ifdef CHARACTER_9_X_1
	val = 1.0L;
	val -= 0.5L/((long double) (n+1));
	val -= (1.0/9.0)*harmonic_n (n+1);
#endif

	// the following sum is patterned on a sub n
	long double acc = 0.0L;
	long double sign = -1.0L;
	for (k=1; k<=n; k++)
	{
		// long double term = zetam1 (k+1)/ (long double) (k+1);
#ifdef CHARACTER_FOUR
		// long double four = gsl_sf_hzeta (k+1, 0.25) + gsl_sf_hzeta (k+1, 0.75);
		long double four = gsl_sf_hzeta (k+1, 0.25) - gsl_sf_hzeta (k+1, 0.75);
		four *= powl (4.0, -(k+1));
		long double term = (four -1.0L)/ ((long double) (k+1));
#endif 
#ifdef CHARACTER_3_1
		long double three = gsl_sf_hzeta (k+1, 1.0/3.0) + gsl_sf_hzeta (k+1, 2.0/3.0);
		three *= powl (3.0, -(k+1));
		long double term = (three -1.0L)/ ((long double) (k+1));
#endif
#ifdef CHARACTER_3_2
		long double three = gsl_sf_hzeta (k+1, 1.0/3.0) - gsl_sf_hzeta (k+1, 2.0/3.0);
		three *= powl (3.0, -(k+1));
		long double term = (three -1.0L)/ ((long double) (k+1));
#endif
#ifdef CHARACTER_5_1
		long double five = gsl_sf_hzeta (k+1, 0.2);
		five += gsl_sf_hzeta (k+1, 0.4);
		five += gsl_sf_hzeta (k+1, 0.6);
		five += gsl_sf_hzeta (k+1, 0.8);
		five *= powl (5.0, -(k+1));
		long double term = (five -1.0L)/ ((long double) (k+1));
#endif
#ifdef CHARACTER_5_2_OOPS
		long double five = gsl_sf_hzeta (k+1, 0.2);
		five += gsl_sf_hzeta (k+1, 0.8);
		five *= powl (5.0, -(k+1));
		long double term = (five -1.0L)/ ((long double) (k+1));
#endif
#ifdef CHARACTER_5_X_1
		long double five = gsl_sf_hzeta (k+1, 0.2);
		five *= powl (5.0, -(k+1));
		long double term = (five -1.0L)/ ((long double) (k+1));
#endif
#ifdef CHARACTER_7_X_1
		long double seven = gsl_sf_hzeta (k+1, 1.0/7.0);
		seven *= powl (7.0, -(k+1));
		long double term = (seven -1.0L)/ ((long double) (k+1));
#endif
#ifdef CHARACTER_9_X_1
		long double nine = gsl_sf_hzeta (k+1, 1.0/9.0);
		nine *= powl (9.0, -(k+1));
		long double term = (nine -1.0L)/ ((long double) (k+1));
#endif
		// long double term = four/ ((long double) (k+1));
		term *= binomial (n,k);
		term *= sign;
		acc += term;
		// printf ("duuude a_sub_n k=%d term=%Lg, acc=%Lg\n", k, term, acc);
		sign = -sign;
	}
	// printf ("finally asub_n=%Lg+%Lg\n",val, -acc);
	return val-acc;
}

// ==========================================================
// finite terms, per the Norlund-Rice integral
static long double norlund_b_sub_n (int n, int m_idx, int k_order)
{
	long double val = 0.0;
	
	long double x = ((long double) m_idx)/((long double) k_order);
	val = gsl_sf_psi(x) + logl(k_order) + 1.0L;
	val -= harmonic_n (n-1);
	val *= ((long double) n)/((long double) k_order);
	val = -val;

	// val = - hurwitz_zero (m_idx, k_order) / ((long double) n+1);
	val += (x-0.5);

	return val;
}

// Return b_sub_n but for Hurwitz zeta
// 
long double 
hurl_b_sub_n (int n, int m, int k)
{
	int p;
	long double val = 0.0L;
	val = norlund_b_sub_n (n,m,k);
	
	// the following sum is patterned on b sub n
	long double acc = 0.0L;
	long double sign = 1.0L;
	for (p=2; p<=n; p++)
	{
		long double term = gsl_sf_hzeta (p, ((double) m)/ ((double) k));
		term *= powl (k, -p);

		term *= binomial (n,p);
		term *= sign;
		acc += term;
		// printf ("duuude b_sub_n k=%d term=%Lg, acc=%Lg\n", p, term, acc);
		sign = -sign;
	}
	// printf ("finally asub_n=%Lg+%Lg\n",val, -acc);
	return acc-val;
}

// ==========================================================
// Hurwitz zeta at s=0, should equal 0.5-m/k
long double hurwitz_zero (int m, int k)
{
	long double val = 0.0;
	int n;
	long double ok = 1.0L/ ((long double) k);
	long double theta = 2.0L*M_PI*m * ok;
	for (n=1; n<= k; n++)
	{
		val += sinl (theta*n) * gsl_sf_psi (n*ok);
	}
	val = -val / (M_PI*k);
	return val;
}

// finite terms, per the Norlund-Rice integral
static long double norlund_a_sub_n (int n, int m_idx, int k_order)
{
	long double val = 0.0;
	
	long double x = ((long double) m_idx)/((long double) k_order);
	val = gsl_sf_psi(x) + logl(k_order) + 1.0L;
	val -= harmonic_n (n);
	val /= (long double) k_order;
	val = -val;

	// val = - hurwitz_zero (m_idx, k_order) / ((long double) n+1);
	val += (x-0.5) / ((long double) n+1);

	return val;
}

// Return a_sub_n but for Hurwitz zeta
// 
long double 
hurwitz_a_sub_n (int n, int m_idx, int k_order)
{
	int k;

	long double kay_order = k_order;	
	long double em_idx = m_idx;
	long double val;

	val = -norlund_a_sub_n (n, m_idx, k_order);

	// the following sum is patterned on a sub n
	long double acc = 0.0L;
	long double sign = -1.0L;
	for (k=1; k<=n; k++)
	{
		long double nine = gsl_sf_hzeta (k+1, em_idx/kay_order);
		nine *= powl (kay_order, -(k+1));
		long double term = nine/ ((long double) (k+1));

		term *= binomial (n,k);
		term *= sign;
		acc += term;
		// printf ("duuude a_sub_n k=%d term=%Lg, acc=%Lg\n", k, term, acc);
		sign = -sign;
	}
	// printf ("finally asub_n=%Lg+%Lg\n",val, -acc);

	return val-acc;
}

// Return a_sub_n but for Dirichlet eta
// 
long double 
eta_a_sub_n (int n)
{
	int k;

	long double val = 1.0L;
	val -= 0.5L/((long double) (n+1));
	val -= logl (2.0L);

	// the following sum is patterned on a sub n
	long double acc = 0.0L;
	long double sign = -1.0L;
	for (k=1; k<=n; k++)
	{
		long double eta = gsl_sf_zeta (k+1);
		eta *= 1.0 - pow (0.5, k);
		long double term = (eta -1.0L)/ ((long double) (k+1));

		term *= binomial (n,k);
		term *= sign;
		acc += term;
		// printf ("duuude a_sub_n k=%d term=%Lg, acc=%Lg\n", k, term, acc);
		sign = -sign;
	}
	// printf ("finally asub_n=%Lg+%Lg\n",val, -acc);
	return val-acc;
}

// Return a_sub_n but for Dirichlet eta
// 
long double 
x_sub_n (int n)
{
	int k;

	// the following sum is patterned on a sub n
	long double acc = 0.0L;
	long double sign = -1.0L;
	for (k=1; k<=n; k++)
	{
		long double term = 1.0/ ((long double) (k+1));

		term *= binomial (n,k);
		term *= sign;
		acc += term;
		sign = -sign;
	}
	return acc;
}

