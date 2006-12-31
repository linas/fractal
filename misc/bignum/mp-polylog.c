
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

/* ============================================================= */
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

/**
 * polylog_borwein() -- Return polylog, Li_s(z) for estimator n.
 *
 * Uses the Borwein algorithm to estimate the polylog value.
 * A polynomial of order 2n is used to perform the estimation.
 * The value of n must be suitably large.
 * Muse have |z^2/(z-1)| < 3.7 in order for the thing to converge
 * to an accurate answer. 
 *
 * Appears to work well. Suggest n=31 for most cases,
 * should return answers accurate to 1e-16
 */
static void polylog_borwein (cpx_t plog, const cpx_t ess, const cpx_t zee, int norder, int prec)
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
		cpx_ui_pow_cache (term, k, s, prec);

		/* Put it together */
		cpx_mul (term, term, pz);
		cpx_add (acc, acc, term);

		cpx_mul(pz, pz, z);
	}

	for (k=norder+1; k<=2*norder; k++)
	{
		/* The inverse integer power */
		cpx_ui_pow_cache (term, k, s, prec);

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

/* ============================================================= */

/* polylog_get_zone -- return | z^2 / (z-1) |^2
 *
 * The value of | z^2 / (z-1) | is used to determine
 * the convergence zone.
 */
inline static double polylog_get_zone (double zre, double zim)
{
	double den = 1.0 / ((zre-1.0)*(zre-1.0) + zim*zim);
	double sre = zre*zre - zim*zim;
	double sim = 2.0*zre*zim;
	double fre = sre * (zre-1.0) + zim*sim;
	double fim = sim * (zre-1.0) - zim*sre;
	den = (fre*fre + fim*fim)*den*den;

	return den;
}

inline static double polylog_modsq (const cpx_t zee)
{
	double zre = mpf_get_d (zee[0].re);	
	double zim = mpf_get_d (zee[0].im);	

	double den = zre*zre + zim*zim;

	return den;
}

/*
 * polylog_terms_est() -- estimate number of terms needed 
 * in the polylog summation in order to keep the error
 * to be less than 10^-prec.
 *
 * The estimation is based on the  development in the 
 * paper, but is broken for values of ess which are 
 * non-positive integers, since then, the gamma explodes
 */
static int polylog_terms_est (const cpx_t ess, const cpx_t zee, int prec)
{
#if 1
	double fterms = 2.302585 * prec;  /* log(10) */

	/* Estimate for the gamma. A slightly better estimate
	 * can be obtains for sre negative but still small. 
	 */
	double gamterms;
	double sre = mpf_get_d (ess[0].re);	
	double sim = mpf_get_d (ess[0].im);	
	if (0.0 > sim) sim = -sim;
	if (0.0 < sre) {
		gamterms = 0.5*M_PI*sim;
	} else {
		gamterms = M_PI*sim;
	}
	gamterms -= lgamma(sre);
	fterms += gamterms;

	/* This part of the estimate not needed, since
	 * we assume that zre =< 0.0 always; if not then
	 * there are other, more severe problems.
	 *
	 * double zre = mpf_get_d (zee[0].re);	
	 * double zim = mpf_get_d (zee[0].im);	
	 * if ((0.0 < zre) && (zre < 1.0)) fterms += log (zim);
	 */

	/* den = | z^2/(z-1)|^2 */
	double zre = mpf_get_d (zee[0].re);	
	double zim = mpf_get_d (zee[0].im);	
	double den = polylog_get_zone (zre, zim);

	fterms /= -0.5*log(den) + 1.345746719;  /* log (0.260345491) */

	int nterms = (int) (fterms+1.0);

	if (nterms > prec*10+150)
	{
		fprintf (stderr, "Truncate nterms=%d at z=%g + i %g\n", nterms, zre, zim);
		nterms = 10*prec;
	}
// gamterms /=  -0.5*log(den) +1.345746719;
// printf ("# duude z= %g +i %g den=%g  nterms = %d gam=%g\n", zre, zim, sqrt(den), nterms, gamterms);
#endif

	// XXX this is a really crappy estimate
// int nterms = 30+0.8*prec;
	return nterms;
}

static void recurse_polylog (cpx_t plog, const cpx_t ess, const cpx_t zee, int prec, int depth);
		  
static inline void polylog_recurse_sqrt (cpx_t plog, const cpx_t ess, const cpx_t zee, int prec, int depth)
{
	cpx_t zroot, s, pp, pn;
	cpx_init (zroot);
	cpx_init (s);
	cpx_init (pp);
	cpx_init (pn);

	cpx_sqrt (zroot, zee, prec);
	cpx_set (s, ess);

#if 0
cpx_prt ("zee= ", zee);
printf ("\n");
cpx_prt ("zroot= ", zroot);
printf ("\n");
#endif
	recurse_polylog (pp, s, zroot, prec, depth);

	cpx_neg (zroot, zroot);
	recurse_polylog (pn, s, zroot, prec, depth);

	cpx_add (plog, pp, pn);

	/* now, compute 2^{s-1} in place */
	cpx_sub_ui (s, s, 1, 0);
	cpx_ui_pow (s, 2, s, prec);
	cpx_mul (plog, plog, s);

	cpx_clear (s);
	cpx_clear (pp);
	cpx_clear (pn);
	cpx_clear (zroot);
}

static inline void polylog_recurse_duple (cpx_t plog, const cpx_t ess, const cpx_t zee, int prec, int depth)
{
	cpx_t zsq, s, pp, pn;
	cpx_init (zsq);
	cpx_init (s);
	cpx_init (pp);
	cpx_init (pn);

	cpx_mul (zsq, zee, zee);
	cpx_set (s, ess);

#if 0
cpx_prt ("dupl zee= ", zee);
printf ("\n");
cpx_prt ("zsq= ", zsq);
printf ("\n");
#endif
	recurse_polylog (pp, s, zsq, prec, depth);

	cpx_neg (zsq, zee);
	recurse_polylog (pn, s, zsq, prec, depth);

	/* now, compute 2^{1-s} in place */
	cpx_sub_ui (s, s, 1, 0);
	cpx_neg (s, s);
	cpx_ui_pow (s, 2, s, prec);
	cpx_mul (plog, pp, s);

	cpx_sub (plog, plog, pn);
	
	cpx_clear (s);
	cpx_clear (pp);
	cpx_clear (pn);
	cpx_clear (zsq);
}

static inline void polylog_recurse_triple (cpx_t plog, const cpx_t ess, const cpx_t zee, int prec, int depth)
{
	cpx_t zcu, s, tr, pp, pu, pd;
	cpx_init (zcu);
	cpx_init (s);
	cpx_init (tr);
	cpx_init (pp);
	cpx_init (pu);
	cpx_init (pd);

	cpx_set (s, ess);

	cpx_mul (zcu, zee, zee);
	cpx_mul (zcu, zcu, zee);
	recurse_polylog (pp, s, zcu, prec, depth);

	/* tr = exp (i 2pi/3) = -1/2  + i sqrt(3)/2 */
	mpf_set_ui (tr[0].re, 1);
	mpf_div_ui (tr[0].re, tr[0].re, 2);
	mpf_neg (tr[0].re, tr[0].re);
	fp_half_sqrt_three (tr[0].im);
	
	cpx_mul (zcu, tr, zee);
	recurse_polylog (pu, s, zcu, prec, depth);

	cpx_mul (zcu, tr, zcu);
	recurse_polylog (pd, s, zcu, prec, depth);

	/* now, compute 3^{1-s} in place */
	cpx_sub_ui (s, s, 1, 0);
	cpx_neg (s, s);
	cpx_ui_pow (s, 3, s, prec);
	cpx_mul (plog, pp, s);

	cpx_sub (plog, plog, pu);
	cpx_sub (plog, plog, pd);
	
	cpx_clear (s);
	cpx_clear (tr);
	cpx_clear (pp);
	cpx_clear (pu);
	cpx_clear (pd);
	cpx_clear (zcu);
}

inline static int accept (double zre, double zim)
{
	double mod = zre*zre + zim*zim;
	if (0.3>mod) return 1;
	if ((0.85>mod) && (0.0 <zre)) return 1;
	return 0;
}

static void recurse_polylog (cpx_t plog, const cpx_t ess, const cpx_t zee, int prec, int depth)
{
	double zre = mpf_get_d (zee[0].re);	
	double zim = mpf_get_d (zee[0].im);	

	if (5 < depth)
	{
		fprintf (stderr, "ecessive recusrsion at z=%g+ i%g\n", zre, zim);
		cpx_set_ui (plog, 0,0);
		return;
	}
	depth ++;

	/* The zone of convergence for the Borwein algorithm is 
	 * |z^2/(z-1)| < 3.  If z is within this zone, then all is 
	 * well, otherwise, use the duplication formula to make 
	 * recusrsive calls, until the leaves of the recursion 
	 * in in this zone.
	 *
	 * Two types of recursion to be applied: 
	 * If |z| > 1, use square-root to move to Borwein.
	 * If |z| < 0.9, use the simple series summation
	 * If 1.01 > |z| > 0.9, square away from z=1
	 */
	if (accept (zre,zim))
	{
		cpx_polylog_sum (plog, ess, zee, prec);
		return;
	}

	/* den = | z^2/(z-1)|^2 */
	double den = polylog_get_zone (zre, zim);

#if 0
	if (den > 10.0)
	{
// printf ("splitsville, den=%g\n", den);
double dre = zre*zre - zim*zim;
double dim = 2.0 *zre*zim;
double d1 = polylog_get_zone (dre, dim);
double d2 = polylog_get_zone (-dre, -dim);
		if (((d1 < den)||accept(dre, dim)) && ((d2 < den)||accept(-dre,-dim)))
		// if ((d1 < den) && (d2 < den))
		{
			polylog_recurse_dupl (plog, ess, zee, prec, depth);
			return;
		}
double mod = 1.0 / sqrt(zre*zre+zim*zim);
zre *= mod;
zim *= mod;
dre = sqrt (0.5*(zre + 1.0));
dim = 0.5*zim / dre;
d1 = polylog_get_zone (dre, dim);
d2 = polylog_get_zone (-dre, -dim);
		if (((d1 < den)||accept(dre, dim)) && ((d2 < den)||accept(-dre,-dim)))
		// if ((d1 < den) && (d2 < den))
		{
			polylog_recurse_sqrt (plog, ess, zee, prec, depth);
			return;
		}
		fprintf (stderr, "no convergence at z=%g+ i%g\n", zre, zim);
		cpx_set_ui (plog, 0,0);
		return;
	}
#endif

	if (den > 10.0)
	{
// printf ("splitsville, den=%g\n", den);
double mod = 1.0 / sqrt(zre*zre+zim*zim);
zre *= mod;
zim *= mod;
double dre = sqrt (0.5*(zre + 1.0));
double dim = 0.5*zim / dre;
double d1 = polylog_get_zone (dre, dim);
double d2 = polylog_get_zone (-dre, -dim);
		if (((d1 < den)||accept(dre, dim)) && ((d2 < den)||accept(-dre,-dim)))
		// if ((d1 < den) && (d2 < den))
		{
			polylog_recurse_sqrt (plog, ess, zee, prec, depth);
			return;
		}
		polylog_recurse_triple (plog, ess, zee, prec, depth);
		return;
	}
	
	int nterms = polylog_terms_est (ess, zee, prec);
	polylog_borwein (plog, ess, zee, nterms, prec);
}

void cpx_polylog (cpx_t plog, const cpx_t ess, const cpx_t zee, int prec)
{
	recurse_polylog (plog, ess, zee, prec, 0);
}

/* ============================================================= */
/**
 * cpx_periodic_zeta -- Periodic zeta function 
 *
 * F(s,q) = sum_{n=1}^infty exp(2pi iqn)/ n^s
 *        = Li_s (exp(2pi iq))
 * where 
 * Li_s(z) is the polylogarithm
 *
 * Periodic zeta function is defined as F(s,q) by Tom Apostol, chapter 12
 *
 * For 0<q<1/4 or for 3/4<q<1, this algo applies the duplication
 * formula, and calls itself recusrively, until 1/4<q<3/4, which
 * can then be evaluated at a single shot using teh Borwein algorithm.
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
		
		// cpx_polylog (z, s, z, prec);
		int nterms = polylog_terms_est (s, z, prec);
		polylog_borwein (z, s, z, nterms, prec);
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
	static cpx_t cache_s, scale;
	static int init = 0;

	if (!init)
	{
		init = 1;
		cpx_init (cache_s);
		cpx_init (scale);
	}

	if (!cpx_eq (ess, cache_s, prec*3.322))
	{
		cpx_set (cache_s, ess);

		mpf_t two_pi;
		mpf_init (two_pi);

		cpx_t s, tps;
		cpx_init (s);
		cpx_init (tps);

		/* 2 gamma(s+1)/ (2pi)^s */
		cpx_add_ui (s, ess, 1,0);
		cpx_gamma (scale, s, prec);

		/* times (2pi)^{-s} */
		fp_two_pi (two_pi, prec);
		cpx_neg (s, ess);
		cpx_mpf_pow (tps, two_pi, s, prec); 
		cpx_mul (scale, scale, tps);

		/* times two */
		cpx_mul_ui (scale, scale, 2);
		cpx_clear (tps);
		cpx_clear (s);
		mpf_clear (two_pi);
	}

	cpx_periodic_zeta (zee, ess, que, prec);
	cpx_mul (zee, zee, scale);
}

/* ============================================================= */
/**
 * cpx_hurwitz_zeta -- Hurwitz zeta function
 *
 * Built up from the periodic zeta. Caches intermediate terms, and so
 * performance is much better if s is held const, while q is varied.
 */

void cpx_hurwitz_zeta (cpx_t zee, const cpx_t ess, const mpf_t que, int prec)
{
	static cpx_t cache_s, piss, niss, scale;
	static int init = 0;

	if (!init)
	{
		init = 1;
		cpx_init (cache_s);
		cpx_init (piss);
		cpx_init (niss);
		cpx_init (scale);
	}

	mpf_t t;
	mpf_init (t);

	cpx_t s, zm;
	cpx_init (s);
	cpx_init (zm);

	/* s = 1-ess */
	cpx_neg (s, ess);
	cpx_add_ui (s, s, 1, 0);

	if (!cpx_eq (s, cache_s, prec*3.322))
	{
		cpx_set (cache_s, s);

		cpx_t tps;
		cpx_init (tps);

		/* exp (i pi s/2) */
		fp_pi_half (t, prec);
		cpx_mul_mpf (piss, s, t);
		cpx_times_i (piss, piss);
		cpx_exp (piss, piss, prec);
		cpx_recip (niss, piss);

		/* 2 gamma(s)/ (2pi)^s */
		cpx_gamma (scale, s, prec);

		/* times (2pi)^{-s} */
		fp_two_pi (t, prec);
		cpx_neg (s, s);
		cpx_mpf_pow (tps, t, s, prec); 
		cpx_neg (s, s);
		cpx_mul (scale, scale, tps);

		/* times two */
		cpx_mul_ui (scale, scale, 2);
		cpx_clear (tps);
	}
	
	/* F(s,q) and F(s, 1-q) */
	cpx_periodic_zeta (zee, s, que, prec);
	mpf_ui_sub (t, 1, que);
	cpx_periodic_zeta (zm, s, t, prec);

	/* assemble the thing */
	cpx_mul (zm, zm, piss);
	cpx_mul (zee, zee, niss);
	cpx_add (zee, zee, zm);
	cpx_mul (zee, zee, scale);

	cpx_clear (s);
	cpx_clear (zm);
	mpf_clear (t);
}

/* ============================================================= */
/**
 * cpx_polylog_sum -- compute the polylogarithm by direct summation
 *
 * Caches intermediate results, so that overall performance is
 * considerably better if z is varied while s is held fixed.
 *
 * The magnitude of z must be less than one.
 */

void cpx_polylog_sum (cpx_t plog, const cpx_t ess, const cpx_t zee, int prec)
{
	int n;

	cpx_t s, z, zp, term;
	cpx_init (s);
	cpx_init (z);
	cpx_init (zp);
	cpx_init (term);

	cpx_set (s, ess);
	cpx_set (z, zee);
	cpx_set (zp, zee);
	cpx_set_ui (plog, 0, 0);

	/* Estimate the number of terms needed to sum over */
	double mag = polylog_modsq (zee);
	
	/* Domain error, should be less than one */
	if (1.0 <= mag)
	{
		fprintf (stderr, "cpx_polylog_sum(): Domain error, |z|=%g\n", sqrt(mag)); 
		return;
	}
	
	int nterms = -2.0 * prec *2.302585093 / log(mag);
	for (n=1; n<nterms; n++)
	{
		cpx_ui_pow_cache (term, n, s, prec);
		cpx_div (term, zp, term);

		cpx_add (plog, plog, term);
		cpx_mul (zp, zp, z);
	}

	cpx_clear (s);
	cpx_clear (z);
	cpx_clear (zp);
	cpx_clear (term);
}

/* ============================================================= */
/**
 * cpx_polylog_nint -- compute the polylogarithm at negetive integers
 *
 * At the negative integers, the polylog is a rational function,
 * meromorphic everywhere except for multiple poles at z=1.
 */

void cpx_polylog_nint (cpx_t plog, unsigned int negn, const cpx_t zee)
{
	int k;

	mpz_t stir, fac;
	mpz_init (stir);
	mpz_init (fac);

	cpx_t z, zp, term;
	cpx_init (z);
	cpx_init (zp);
	cpx_init (term);

	cpx_set (z, zee);
	cpx_sub_ui (zp, zee, 1, 0);
	cpx_recip (zp, zp);

	if (0 == negn)
	{
		cpx_mul (plog, z, zp);
		cpx_neg (plog, plog);
	}
	else
	{
		cpx_set_ui (plog, 0, 0);
		mpz_set_ui (fac, 1);
		cpx_set (z, zp);
		for (k=1; k<= negn+1; k++)
		{
			i_stirling_second (stir, negn+1, k);
			mpz_mul (stir, stir, fac);
			mpf_set_z (term[0].re, stir);
			mpf_set_ui (term[0].im, 0);

			cpx_mul (term, term, zp);
	
			cpx_add (plog, plog, term);
			cpx_mul (zp, zp, z);
			mpz_mul_ui (fac, fac, k);
		}

		if (0==negn%2)
		{
			cpx_neg (plog, plog);
		}
	}

	cpx_clear (z);
	cpx_clear (zp);
	cpx_clear (term);
	mpz_clear (stir);
	mpz_clear (fac);
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
