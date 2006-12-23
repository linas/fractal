
/* polylog.c
 *
 * Implement Borwein-style polylogarithm 
 * on unit circle.
 *
 * Incomplte, partly functional.
 *
 * Linas November 2006
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sf_zeta.h>

#include "binomial.h"
#include "bernoulli.h"

typedef struct {
	double re;
	double im;
} cplex;

static inline cplex cplex_zero (void)
{
	cplex z; z.re=0.0; z.im=0.0; return z;
}
static inline cplex cplex_one (void)
{
	cplex z; z.re=1.0; z.im=0.0; return z;
}

static inline cplex cplex_neg (cplex z)
{
	cplex rv;
	rv.re = -z.re;
	rv.im = -z.im;
	return rv;
}

static inline cplex cplex_add (cplex a, cplex b)
{
	cplex rv;
	rv.re = a.re + b.re;
	rv.im = a.im + b.im;
	return rv;
}

static inline cplex cplex_sub (cplex a, cplex b)
{
	cplex rv;
	rv.re = a.re - b.re;
	rv.im = a.im - b.im;
	return rv;
}

static inline cplex cplex_scale (double x, cplex z)
{
	cplex rv;
	rv.re = x * z.re;
	rv.im = x * z.im;
	return rv;
}

static inline cplex cplex_mult (cplex a, cplex b)
{
	cplex rv;
	rv.re = a.re * b.re - a.im * b.im;
	rv.im = a.re * b.im + b.re * a.im;
	return rv;
}

static inline cplex cplex_recip (cplex z)
{
	cplex rv;
	double mag = 1.0/ (z.re*z.re + z.im*z.im);
	rv.re = z.re * mag;
	rv.im = -z.im * mag;
	return rv;
}

static inline double cplex_modulus (cplex z)
{
	return sqrt (z.re*z.re + z.im*z.im);
}

static inline double cplex_phase (cplex z)
{
	return atan2 (z.im, z.re);
}

static inline cplex cplex_pow (cplex z, int n)
{
	cplex rv;
	if (350>n)
	{
		int k;
		rv = cplex_one ();
		for (k=0; k<n; k++)
		{
			rv = cplex_mult (rv, z);
		}
	}
	else
	{
		fprintf (stderr, "pow=%d unimplemented\n", n);
	}
	return rv;
}

static inline cplex cplex_exp (cplex z)
{
	cplex rv;
	rv.re = rv.im = exp (z.re);
	rv.re *= cos (z.im);
	rv.im *= sin (z.im);
	return rv;
}

/* 
 * bee_k() 
 * Return value of sum_{j=0}^k (n j) oz^j
 *
 * where (n j) is binomial coefficient 
 */
static cplex bee_k (int n, int k, cplex oz)
{
	int j;
	cplex pz = cplex_one();
	cplex acc = cplex_zero();

	for (j=0; j<=k; j++)
	{
		cplex term = pz;
		double bin = binomial (n,j);

		term = cplex_scale (bin, term);
		acc = cplex_add (acc, term);
		pz = cplex_mult(pz, oz);
	}

	return acc;
}

/* polylog_est() 
 * Return estimate of polylog, Li_s(z) for estimator n.
 *
 * Appears to work well. Suggest n=31 for most cases,
 * should return answers accurate to 1e-16
 */
static cplex polylog_est (cplex s, cplex z, int n)
{
	int k;

	/* oz = 1/z   whereas moz = -1/z */
	cplex oz = cplex_recip(z);
	cplex moz = cplex_neg(oz);

	/* ska = [z/(z-1)]^n */
	cplex ska = z;
	ska.re -= 1.0;
	ska = cplex_mult (z, cplex_recip(ska));
	ska = cplex_pow (ska,n);
	
	cplex pz = z;
	cplex acc = cplex_zero();
	cplex cacc = cplex_zero();

	for (k=1; k<=n; k++)
	{
		/* The inverse integer power */
		cplex dir = cplex_scale (log(k), s);
		dir = cplex_exp (cplex_neg(dir));

		/* Put it together */
		cplex term = cplex_mult (pz, dir);
		acc = cplex_add (acc, term);

		pz = cplex_mult(pz, z);
	}

	for (k=n+1; k<=2*n; k++)
	{
		/* The inverse integer power */
		cplex dir = cplex_scale (log(k), s);
		dir = cplex_exp (cplex_neg(dir));

		/* Put it together */
		cplex term = cplex_mult (pz, dir);
		acc = cplex_add (acc, term);

		/* The coefficient c_k of the polynomial */
		cplex ck = bee_k (n,k-n-1,moz);
		term = cplex_mult (ck, term);
		cacc = cplex_add (cacc, term);

		pz = cplex_mult(pz, z);
	}

	cacc = cplex_mult (cacc, ska);
	acc = cplex_sub (acc, cacc);

	return acc;
}

int nt = 31;

/**
 * periodic_zeta -- Implement periodic zeta function
 *
 * Periodic zeta function is defined as F(s,q) by Tom Apostol, chapter 12
 *
 * Incomplete implementation, n3eeds to 
 */
static cplex periodic_zeta (cplex s, double q)
{
	int nterms =nt;
	cplex z;

	if ((0.0>q) || (1.0<q))
	{
		q -= floor (q);
	}
	
	if ((1.0e-10 > q) || (1.0e-10 > 1.0-q))
	{
		// XXX this is wrong, it should return riemann zeta 
		return cplex_zero();
	}
	else if (0.25 > q) 
	{
		/* use the duplication formula to get into convergent region */
		cplex sm = cplex_neg(s);
		sm.re += 1.0;
		cplex ts = cplex_exp(cplex_scale (log(2), sm));
		cplex bt = periodic_zeta (s, 2.0*q);
		bt = cplex_mult (bt, ts);

		cplex z = periodic_zeta (s, q+0.5);
		z = cplex_sub(bt, z);
		return z;

	}
	else if (0.75 < q) 
	{
		/* use the duplication formula to get into convergent region */
		cplex sm = cplex_neg(s);
		sm.re += 1.0;
		cplex ts = cplex_exp(cplex_scale (log(2), sm));
		cplex bt = periodic_zeta (s, 2.0*q-1.0);
		bt = cplex_mult (bt, ts);

		cplex z = periodic_zeta (s, q-0.5);
		z = cplex_sub(bt, z);
		return z;

	}
	else
	{
		/* Normal case; within the convergence region */
		double ph = 2.0*M_PI*q;
		z.re = cos (ph);
		z.im = sin (ph);
		return polylog_est (s, z, nterms);
	}
}

/* ============================================================= */

/** 
 * test_zeta_values() -- compare periodic zeta to reiman zeta
 * 
 * As of 22 December 2006, this test is passing, with flying colors
 * Explores value of hurwitz zeta on s=real line, for 
 * q=1/2, where it can be compared to the Riemann zeta.
 * Passes, very nicely and cleanly, (i.e. error of order 1e-16)
 * although starts to get rough for the large negative s.
 */
void test_zeta_values (double max)
{
	cplex zl, zh;

	cplex s;
	s.im = 0.0;
	for (s.re = -2; s.re < max; s.re += 0.1)
	{
		if (1 == s.re) continue;

		zl = periodic_zeta (s, 0.5);
		
		/* sm = 1-s */
		cplex sm = cplex_neg(s);
		sm.re += 1.0;

		/* ts = 2^(1-s) */
		cplex ts = cplex_exp(cplex_scale (log(2), sm));
		
		/* ots = -1/(1-2^(1-s)) */
		cplex ots = cplex_neg(ts);
		ots.re += 1.0;
		ots = cplex_recip(ots);
		ots = cplex_neg (ots);
		
		zl = cplex_mult (zl, ots);
		
		double zeta = gsl_sf_zeta (s.re);
		
		printf ("s=%5.3g	algo=%12.10g	exact=%12.10g	diff=%6.3g\n", s.re, zl.re, zeta, zl.re-zeta);
	}
}

/* ============================================================= */

/** 
 * test_bernoulli_poly - compare periodic zeta to the Bernoulli poly's
 *
 * The Bernoulli polynomials have a fourier transform that is the 
 * periodic zeta function. 
 *
 * Test is mostly right except for n%4 == zero, where a sign goes 
 * crazy. WTF !!??
 */
int test_bernoulli_poly (int n)
{
	cplex zl, zh;

	cplex s, z;
	s.im = 0.0;
	s.re = n;
	double q;
	for (q=-0.2; q<=1.2; q+=0.02)
	{
		zl = periodic_zeta (s, q);
		zh = periodic_zeta (s, 1.0-q);
		if (n%2) {
			z = cplex_sub (zl,zh);
		} else {
			z = cplex_add (zl,zh);
		}
		
		double bs = z.re;
		if (n%2)
		{
			bs = z.im;
			if (n%4 == 3) bs = -bs;
		}
		if (n%2 == 0) bs = -bs;

		bs *= factorial (n) * pow (2.0*M_PI, -n);
		
		// double b = q*q-q+1.0/6.0;
		double b = bernoulli_poly (n,q);
		
		printf ("q=%5.3g	bs=%13.10g	bernoulli=%13.10g	reldiff=%6.3g\n", q, bs, b, (bs-b)/b);
	}
}

/* ============================================================= */

main (int argc, char * argv[])
{
	int n;

	if (argc != 2)
	{
		fprintf (stderr, "Usage: %s <n>\n", argv[0]);
		exit (1);
	}
	n = atoi (argv[1]);

	// test_zeta_values (n);
	test_bernoulli_poly (n);

// #define HURWITZ_ZETA
#ifdef HURWITZ_ZETA
	cplex s;
	s.im = 0.0;
	double q=0.4;
	for (s.re = 1.1; s.re < 8; s.re += 0.1)
	{
		zl = periodic_zeta (s, q);
		zh = periodic_zeta (s, 1.0-q);
		
		/* wrong ... */
		double zeta = gsl_sf_hzeta (s.re, q);
		
		printf ("%7.3g	%15.10g	%15.10g	%6.3g\n", s.re, zl.re, zeta, zl.re-zeta);
	}
#endif

#if 0
	cplex s;
	s.re = 6.0;
	s.im = 0.0;

	double q;
	for (q=0.0; q<1.0; q+=0.02)
	{
		nt=41;
		zl = periodic_zeta (s, q);
		double b = q*q-q+1.0/6.0;
		printf ("%7.3g	%15.10g	%15.10g\n", q, zl.re, zl.im);
	}
#endif
}
