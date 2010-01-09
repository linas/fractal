
/** 
 * polylog-test.c
 *
 * Use/test polylogarithm, &the Hurwitz zeta function.
 *
 * Linas November 2006
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_zeta.h>

#include "binomial.h"
#include "cplex.h"
#include "bernoulli.h"
#include "harmonic.h"
#include "polylog.h"

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
		
		zl = cplex_mul (zl, ots);
		
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
 * Test is now passijng with flying colors
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

#if 1
	cplex s;
	s.re = 1.5;
	s.im = en;

	printf ("#\n# periodic zeta for s=%g +i %g\n#\n", s.re, s.im);
	fflush (stdout);

	double q;
	for (q=0.002; q<1.0; q+=0.004)
	{
		// cplex hz = hurwitz_zeta (s, q);
		cplex hz = periodic_zeta_sum (s, q);
		printf ("%7.3g	%15.10g	%15.10g\n", q, hz.re, hz.im);
		fflush (stdout);
	}
#endif
}
