
/** 
 * mp-polylog.c
 *
 * Implement Borwein-style polylogarithm.
 * Also implement the "periodic zeta" and 
 * the Hurwitz zeta function.
 *
 * As of 22 December 2006, seems to be fully functional
 * and correct, and passes tests. The range of convergence
 * is rather limited because of precision/rounding errors.
 *
 * Linas November 2006
 */

#include <math.h>
#include <stdio.h>

#include "mp-binomial.h"
#include "mp-complex.h"
#include "mp-consts.h"
#include "mp-gamma.h"
#include "mp-polylog.h"
#include "mp-trig.h"

/* 
 * bee_k() 
 * Return value of sum_{j=0}^k (n j) oz^j
 *
 * where (n j) is binomial coefficient 
 */
static void bee_k (cpx_t bee, int n, int k, const cpx_t oz)
{
	int j;
	cpx_t pz, z, term;
	mpz_t bin;
	mpf_t binom;

	mpz_init (bin);
	mpf_init (binom);
	cpx_init (z);
	cpx_init (pz);
	cpx_init (term);

	/* Make a copy if input arg now! */
	cpx_set (z, oz);
	cpx_set_ui (pz, 1,0);
	cpx_set_ui (bee, 0,0);

	for (j=0; j<=k; j++)
	{
		cpx_set (term, pz);

		i_binomial (bin, n,j);
		mpf_set_z (binom, bin);
		cpx_mul_mpf (term, term, binom);
		cpx_add (bee, bee, term);
		cpx_mul(pz, pz, z);
	}

	mpz_clear (bin);
	mpf_clear (binom);
	cpx_clear (z);
	cpx_clear (pz);
	cpx_clear (term);
}

/* polylog_est() -- Return polylog, Li_s(z) for estimator n.
 *
 * Appears to work well. Suggest n=31 for most cases,
 * should return answers accurate to 1e-16
 */
static void polylog_est (cpx_t plog, const cpx_t ess, const cpx_t zee, int norder, int prec)
{
	cpx_t s, z, oz, moz, ska, pz, acc, term, ck;
	int k;

	cpx_init (s);
	cpx_init (z);
	cpx_init (oz);
	cpx_init (moz);
	cpx_init (ska);
	cpx_init (pz);
	cpx_init (acc);
	cpx_init (term);
	cpx_init (ck);

	/* s = -ess */
	cpx_neg (s, ess);
	cpx_set (z, zee);
	
	/* oz = 1/z   whereas moz = -1/z */
	cpx_recip (oz, z);
	cpx_neg (moz, oz);

	/* ska = [z/(z-1)]^n */
	cpx_set (ska, z);
	cpx_sub_ui (ska, ska, 1, 0);
	cpx_recip (ska, ska);
	cpx_mul (ska, ska, z);
	cpx_pow_ui (ska, ska, norder);
	
	cpx_set (pz, z);
	cpx_set_ui (acc, 0, 0);
	cpx_set_ui (plog, 0, 0);

	for (k=1; k<=norder; k++)
	{
		/* The inverse integer power */
		mpf_set_ui (term[0].re, k);
		cpx_pow_mpf (term, term[0].re, s, prec);

		/* Put it together */
		cpx_mul (term, term, pz);
		cpx_add (acc, acc, term);

		cpx_mul(pz, pz, z);
	}

	for (k=norder+1; k<=2*norder; k++)
	{
		/* The inverse integer power */
		mpf_set_ui (term[0].re, k);
		cpx_pow_mpf (term, term[0].re, s, prec);

		/* Put it together */
		cpx_mul (term, term, pz);
		cpx_add (acc, acc, term);

		/* The coefficient c_k of the polynomial */
		bee_k (ck, norder, k-norder-1, moz);
		cpx_mul (term, ck, term);
		cpx_add (plog, plog, term);

		cpx_mul(pz, pz, z);
	}

	cpx_mul (plog, plog, ska);
	cpx_sub (plog, acc, plog);
	
	cpx_clear (s);
	cpx_clear (z);
	cpx_clear (oz);
	cpx_clear (moz);
	cpx_clear (ska);
	cpx_clear (pz);
	cpx_clear (acc);
	cpx_clear (term);
	cpx_clear (ck);
}

void cpx_polylog (cpx_t plog, const cpx_t ess, const cpx_t zee, int prec)
{
	// XXX this is a really crappy estimate
	int nterms = 20 + 0.7*prec;
	polylog_est (plog, ess, zee, nterms, prec);
}

/**
 * cpx_periodic_zeta -- Periodic zeta function 
 *
 * F(s,q) = sum_{n=1}^infty exp(2pi iqn)/ n^s
 *        = Li_s (exp(2pi iq))
 * where 
 * Li_s(z) is the polylogarithm
 *
 * Periodic zeta function is defined as F(s,q) by Tom Apostol, chapter 12
 */
void cpx_periodic_zeta (cpx_t z, const cpx_t ess, const mpf_t que, int prec)
{
	mpf_t q, qf;
	mpf_init (q);
	mpf_init (qf);
	
	cpx_t s, sm;
	cpx_init (s);
	cpx_init (sm);

	mpf_set (q, que);
	mpf_floor (qf, q);
	mpf_sub (q, q, qf);

	cpx_set (s, ess);
	
	// if ((1.0e-10 > q) || (1.0e-10 > 1.0-q))
	if (0)
	{
		// XXX should be more precise with the next order
		// q correction ... 
		// riemann_zeta (s.re, s.im, &z.re, &z.im);
	}
	else if (mpf_cmp_d (q, 0.25) < 0) 
	{
		/* Use the duplication formula to get into convergent region */
		cpx_t ts, bt;
		cpx_init (ts);
		cpx_init (bt);

		/* sm = 1-s */
		cpx_neg (sm, s);
		cpx_add_ui (sm, sm, 1, 0);

		/* ts = 2^{1-s} */
		fp_log2 (qf, prec);
		cpx_mul_mpf (sm, sm, qf);
		cpx_exp (ts, sm, prec);
		
		/* bt = pzeta (2q) * 2^{1-s} */
		mpf_mul_ui (qf, q, 2);
		cpx_periodic_zeta (bt, s, qf, prec);
		cpx_mul (bt, bt, ts);

		/* pzeta (q+0.5) */
		mpf_set_ui (qf, 1);
		mpf_div_ui (qf, qf, 2);
		mpf_add (qf, q, qf);
		cpx_periodic_zeta (z, s, qf, prec);
		cpx_sub (z, bt, z);
		
		cpx_clear (ts);
		cpx_clear (bt);
	}
	else if (mpf_cmp_d (q, 0.75) > 0) 
	{
		/* Use the duplication formula to get into convergent region */
		cpx_t ts, bt;
		cpx_init (ts);
		cpx_init (bt);

		/* sm = 1-s */
		cpx_neg (sm, s);
		cpx_add_ui (sm, sm, 1, 0);

		/* ts = 2^{1-s} */
		fp_log2 (qf, prec);
		cpx_mul_mpf (sm, sm, qf);
		cpx_exp (ts, sm, prec);
		
		/* bt = pzeta (2q-1) * 2^{1-s} */
		mpf_mul_ui (qf, q, 2);
		mpf_sub_ui (qf, qf, 1);
		cpx_periodic_zeta (bt, s, qf, prec);
		cpx_mul (bt, bt, ts);

		/* pzeta (q-0.5) */
		mpf_set_ui (qf, 1);
		mpf_div_ui (qf, qf, 2);
		mpf_sub (qf, q, qf);
		cpx_periodic_zeta (z, s, qf, prec);
		cpx_sub (z, bt, z);

		cpx_clear (ts);
		cpx_clear (bt);
	}
	else
	{
		/* Normal case; within the convergence region */
		fp_two_pi (qf, prec);
		mpf_mul (qf, qf, q);

		fp_cosine (z[0].re, qf, prec);
		fp_sine (z[0].im, qf, prec);
		
		cpx_polylog (z, s, z, prec);
	}
	
	mpf_clear (q);
	mpf_clear (qf);

	cpx_clear (s);
	cpx_clear (sm);
}

/* ============================================================= */
/**
 * cpx_periodic_beta -- Periodic beta function 
 *
 * Similar to periodic zeta, but with different normalization
 *
 * beta = 2 Gamma(s+1) (2\pi)^{-s} F(s,q)
 *
 * As of 22 December, seems to be passing the tests -- 
 * that is, it gives the Bernoulli polynomials for integer s,
 * with all the right scale factors and signs, etc. Yay!
 */
void cpx_periodic_beta (cpx_t zee, const cpx_t ess, const mpf_t que, int prec)
{
	mpf_t twopi;
	mpf_init (twopi);
	fp_two_pi (twopi, prec);
	
	cpx_t s, tps;
	cpx_init (s);
	cpx_init (tps);

	cpx_set (s, ess);
	cpx_periodic_zeta (zee, s, que, prec);

	/* times (2pi)^{-s} */
	cpx_neg (s, s);
	cpx_pow_mpf (tps, twopi, s, prec); 
	cpx_mul (zee, zee, tps);
	cpx_neg (s, s);

	/* times gamma(s+1) */
	cpx_add_ui (s, s, 1,0);
	cpx_gamma (tps, s, prec);
	cpx_mul (zee, zee, tps);
	cpx_mul_ui (zee, zee, 2);
	
	cpx_clear (s);
	cpx_clear (tps);
	mpf_clear (twopi);
}

/* ============================================================= */
/**
 * cpx_hurwitz_zeta -- Hurwitz zeta function
 *
 * Built up from the periodic beta
 */
void cpx_hurwitz_zeta (cpx_t zee, const cpx_t ess, const mpf_t que, int prec)
{
	mpf_t t;
	mpf_init (t);
	
	cpx_t s, tps, zm, piss;
	cpx_init (s);
	cpx_init (tps);
	cpx_init (zm);
	cpx_init (piss);

	/* s = 1-ess */
	cpx_neg (s, ess);
	cpx_add_ui (s, s, 1, 0);

	/* F(s,q) and F(s, 1-q) */
	cpx_periodic_zeta (zee, s, que, prec);
	mpf_ui_sub (t, 1, que);
	cpx_periodic_zeta (zm, s, t, prec);

	/* exp (i pi s/2) */
	fp_pi_half (t, prec);
	cpx_mul_mpf (piss, s, t);
	cpx_times_i (piss, piss);
	cpx_exp (piss, piss, prec);

	/* assemble the thing */
	cpx_mul (zm, zm, piss);
	cpx_div (zee, zee, piss);
	cpx_add (zee, zee, zm);
	
	/* times gamma(s) */
	cpx_gamma (tps, s, prec);
	cpx_mul (zee, zee, tps);
	cpx_mul_ui (zee, zee, 2);
	
	/* times (2pi)^{-s} */
	fp_two_pi (t, prec);
	cpx_neg (s, s);
	cpx_pow_mpf (tps, t, s, prec); 
	cpx_mul (zee, zee, tps);

	cpx_clear (s);
	cpx_clear (tps);
	cpx_clear (zm);
	cpx_clear (piss);
	mpf_clear (t);
}

/* ============================================================= */

/** 
 * test_bernoulli_poly - compare periodic zeta to the Bernoulli poly's
 *
 * The Bernoulli polynomials have a fourier transform that is the 
 * periodic zeta function. 
 *
 * Test is now passing with flying colors
 */
#if 0
int test_bernoulli_poly (int n)
{
	cplex zl, zh;

	cplex s, z;
	s.im = 0.0;
	s.re = n;
	double q;
	for (q=-0.2; q<=1.2; q+=0.02)
	{
		// zl = periodic_zeta (s, q);
		// zh = periodic_zeta (s, 1.0-q);
		zl = periodic_beta (s, q);
		zh = periodic_beta (s, 1.0-q);
		if (n%2) {
			z = cplex_sub (zl,zh);
		} else {
			z = cplex_add (zl,zh);
		}
		
		double bs;
		if (0 == n%2)
		{
	 		bs = z.re;
			if (n%4 == 0) bs = -bs;
		}
		if (n%2)
		{
			bs = -z.im;
			if (n%4 == 3) bs = -bs;
		}

		// bs *= factorial (n) * pow (2.0*M_PI, -n);
		bs *= 0.5;

		/* short-circuit the above, test directly */
		cplex ess;
		ess.re = 1-n;
		ess.im = 0;
		cplex hz = hurwitz_zeta(ess,q);
		bs = -n * hz.re;
		
		// double b = q*q-q+1.0/6.0;
		double b = bernoulli_poly (n,q);
		
		printf ("q=%5.3g	bs=%13.10g	bernoulli=%13.10g	reldiff=%6.3g\n", q, bs, b, (bs-b)/b);
	}
}

/* ============================================================= */
#include <stdlib.h>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_zeta.h>

main (int argc, char * argv[])
{
	int n;
	double en;

	if (argc != 2)
	{
		fprintf (stderr, "Usage: %s <n>\n", argv[0]);
		exit (1);
	}
	n = atoi (argv[1]);
	en = atof (argv[1]);

	// test_zeta_values (n);
	// test_bernoulli_poly (n);

// #define HURWITZ_ZETA
#ifdef HURWITZ_ZETA
	/* As of 22 December 2006, this test seems to be passing, 
	 * with decent accuracy, for anything with real part less than about 8
	 */
	cplex s;
	s.im = 0.0;
	double q=0.5;
	for (s.re = 1.05; s.re < n; s.re += 0.1)
	{
		cplex hz= hurwitz_zeta (s, q);
		
		double zeta = gsl_sf_hzeta (s.re, q);
		
		printf ("s=%5.3g	algo=%12.10g	exact=%12.10g	reldiff=%6.3g\n", s.re, hz.re, zeta, (hz.re-zeta)/zeta);
	}
#endif

}
#endif
