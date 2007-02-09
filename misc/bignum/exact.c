/**
 * exact.c
 *
 * Look at hurtwitz zeta at the location of the Riemann zeta zeros
 *
 * December 2006
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mp-complex.h"
#include "mp-consts.h"
#include "mp-gamma.h"
#include "mp-misc.h"
#include "mp-polylog.h"
#include "mp-trig.h"

void polylog_sheet(cpx_t delta, const cpx_t ess, const cpx_t zee, int sheet, int prec);

#ifdef BORKOKEN_AND_DONT_WORK

#include "mp-zeta.h"

/* ============================================================================ */
/* The cpx_zeta_series() is supposed to sum up to 
 * the periodic zeta function. However, its does not,
 * not sure why -- maybe the formula is wrong, maybe 
 * the evaluation of the reimann zeta at large values in the 
 * second quadrant is wrong. At any rate, the sum
 * doesn't even converge.
 */
void cpx_zeta_series (cpx_t result, const cpx_t ess, const mpf_t que, int prec)
{
	int k;
	cpx_t mu, pmu, term, s, sk;

	mpf_t twopi, fact, arg;
	mpf_init (twopi);
	mpf_init (fact);
	mpf_init (arg);
	fp_two_pi (twopi, prec);
	mpf_set_ui (fact, 1);

	cpx_init (mu);
	cpx_init (pmu);
	cpx_init (term);
	cpx_init (s);
	cpx_init (sk);

	cpx_set (s, ess);

	mpf_mul (arg, que, twopi);
	cpx_set_ui (mu, 0, 1);
	cpx_mul_mpf (mu, mu, arg);
	cpx_set_ui (pmu, 1, 0);

	cpx_set_ui (result, 0, 0);

	for (k=0; k<70; k++)
	{
printf ("duude k=%d\t", k);
		/* zeta (s-k) * mu^k / k! */
		cpx_sub_ui (sk, s, k, 0);
		cpx_borwein_zeta (term, sk, prec);
cpx_prt ("zeta=", term);
printf ("\n");
		cpx_mul (term, term, pmu);
		cpx_mul_mpf (term, term, fact);
		
cpx_prt ("term=", term);
printf ("\n");
		cpx_add (result, result, term);
cpx_prt ("sum=", result);
printf ("\n");
printf ("\n");

		/* pmu = mu^k and fact = 1/k! */
		cpx_mul (pmu, pmu, mu);
		mpf_div_ui (fact, fact, k+1);
	}
cpx_prt ("sum=", result);
printf ("\n");

	/* (-mu)^{s-1} = (2pi q)^{s-1} exp(-i pi (s-1)/2) */
	cpx_sub_ui (sk, s, 1, 0);
	cpx_mpf_pow (term, arg, sk, prec);

	cpx_mul_mpf (sk, sk, twopi);
	cpx_div_ui (sk, sk, 4);
	cpx_neg (sk, sk);
	cpx_exp (sk,sk,prec);
	cpx_mul (term, term, sk);

	/* Gamma (1-s) */
	cpx_sub_ui (sk, s, 1, 0);
	cpx_neg (sk, sk);
	cpx_gamma_cache (pmu, sk, prec);
	cpx_mul (term, term, pmu);

	cpx_add (result, result, term);

	cpx_clear (mu);
	cpx_clear (pmu);
	cpx_clear (term);
	cpx_clear (s);
	cpx_clear (sk);
	mpf_clear (twopi);
	mpf_clear (fact);
	mpf_clear (arg);
}
#endif

/* ============================================================================ */

void sumpos (cpx_t coho, const cpx_t lambda, const cpx_t zee)
{
	int n;

	cpx_t sum, term, zp;
	cpx_init (sum);
	cpx_init (term);
	cpx_init (zp);
	
	cpx_set_ui (sum, 0, 0);
	cpx_set (zp, zee);

	/* sum_{n=1} nz^n / (n-lambda) */
	for (n=1; n<50; n++)
	{
		cpx_sub_ui (term, lambda, n, 0);
		cpx_div (term, zp, term);
		cpx_mul_ui (term, term, n);
		cpx_add (sum, sum, term);
		cpx_mul (zp, zp, zee); 
	}

	cpx_neg (coho, sum);
	
	cpx_clear (sum);
	cpx_clear (term);
	cpx_clear (zp);
}

void cohere (cpx_t coho, const cpx_t lambda, const cpx_t zee, int prec)
{
	int n;

	cpx_t sum, term, lamp, om;
	cpx_init (om);
	cpx_init (sum);
	cpx_init (term);
	cpx_init (lamp);
	
	cpx_set_ui (sum, 0, 0);
	sumpos (sum, lambda, zee);
			  
	cpx_recip (om, lambda);
	cpx_set_ui (lamp, 1, 0);
	for (n=1; n<80; n++)
	{
		cpx_polylog_nint (term, n, zee);
		cpx_mul (term, term, lamp);
printf("n=%d ", -n);
cpx_prt ("neg term= ", term);
printf("\n");
		cpx_add (sum, sum, term);
		cpx_mul (lamp, lamp, om); 
	}

	cpx_set (coho, sum);
	
	cpx_clear (om);
	cpx_clear (sum);
	cpx_clear (term);
	cpx_clear (lamp);
}

/* ============================================================================ */
/* branch discontinuity */

void disco (cpx_t disco, const cpx_t ess, int dir, int prec)
{
	mpf_t p;
	mpf_init (p);

	cpx_t s, term;
	cpx_init (s);
	cpx_init (term);
	
	/* s = (1-s) */
	cpx_set_ui (s, 1, 0);
	cpx_sub (s, s, ess);

	cpx_gamma (disco, s, prec);

	/* times exp(i pi (1-s)/2) */
	fp_pi_half (p, prec);
	cpx_mul_mpf (term, s, p);
	cpx_times_i (term, term);
	if (dir) cpx_neg (term, term);
	cpx_exp (term, term, prec);
	cpx_mul (disco, disco,term);

	/* times (2pi)^{s-1} */
	fp_two_pi (p, prec);
	cpx_neg (s, s);
	cpx_mpf_pow (term, p, s, prec);
	cpx_mul (disco, disco,term);
	
	cpx_clear (s);
	cpx_clear (term);
	mpf_clear (p);
}


/* ============================================================================ */

int
main (int argc, char * argv[])
{
	int prec = 40;
	prec = 20;
	double q;

	/* Set the precision (number of binary bits) */
	int nbits = 3.3*prec+100;
	mpf_set_default_prec (nbits);

	if (argc != 4)
	{
		fprintf (stderr, "Usage: %s <sre> <sim> <zre>\n", argv[0]);
		exit (1);
	}
	double sre = atof (argv[1]);
	double sim = atof (argv[2]);
	double zre = atof (argv[3]);
	double zim = zre;
	
	cpx_t ess, zeta, z2, zee, plog, ph;
	cpx_init (ess);
	cpx_init (zeta);
	cpx_init (z2);
	cpx_init (zee);
	cpx_init (plog);
	cpx_init (ph);

	mpf_t que, tp;
	mpf_init (que);
	mpf_init (tp);
			  
	// cpx_set_d (ess, 1.5, 14.134725);
	cpx_set_d (ess, -1.563331235, sim);
	cpx_set_d (ess, 0.5, sim);
	cpx_set_d (ess, sre, sim);

	mpf_t twopi;
	mpf_init (twopi);
	fp_two_pi (twopi, prec);

	cpx_mul_mpf (ph, ess, twopi);
	cpx_times_i (ph, ph);
	cpx_neg (ph,ph);
	cpx_exp (ph, ph, prec);
	cpx_neg (ph,ph);

#if 0
	char * zero;
	zero = "14.134725141734693790457251983562470270784257115699243175685567460149 \
	        9634298092567649490103931715610127792029715487974367661426914698822545 \
	        8250536323944713778041338123720597054962195586586020055556672583601077";

	zero = "21.022039638771554992628479593896902777334340524902781754629520403587 \
	        5985860688907997136585141801514195337254736424758913838650686037313212 \
	        6211882162437574166925654471184407119403130672564622779261488733743555";
					
	zero = "25.010857580145688763213790992562821818659549672557996672496542006745 \
	        0920984416442778402382245580624407504710461490557783782998515227308011 \
	        8813393358267168958722516981043873551292849372719199462297591267547869";

	mpf_set_str (ess[0].im, zero, 10);
#endif

#if 1

	cpx_t facup, facdown, ms, z3, z4, cq;
	cpx_init (facup);
	cpx_init (facdown);
	cpx_init (ms);
	cpx_init (z3);
	cpx_init (z4);
	cpx_init (cq);
	disco (facup, ess, 0, prec);
	disco (facdown, ess, 1, prec);
	
	// printf ("#\n# graph of Hurwitz zeta as function of q, at \n#\n");
	printf ("#\n# graph of periodic zeta as function of q, at \n#\n");
	fp_prt ("# at s= ", ess[0].re);
	fp_prt (" +i ", ess[0].im);
	printf ("\n#\n# prec=%d nbits=%d\n#\n", prec, nbits);
	fflush (stdout);
	// for (q=0.02; q<0.991; q+=0.008)
	// for (zre=-0.53; zre<1.00; zre+=0.05)
	for (zim=-0.00999999; zim<0.00999999; zim+=0.0000811)
	{
		// mpf_set_d (que, q);
		// cpx_set_d (cq, q, 0.0);
		// cpx_hurwitz_zeta (zeta, ess, que, prec);
		// cpx_hurwitz_zeta (z3, ess, que, prec);
		// cpx_hurwitz_taylor (zeta, ess, cq, prec);
		// cpx_pade_hurwitz_zeta (z2, ess, que, prec);
		// cpx_periodic_beta (zeta, ess, que, prec);
		// cpx_periodic_zeta (zeta, ess, que, prec);
		// cpx_polylog_sum (plog, ess, zee, prec);
		// cpx_set_d (zee, q, 0.002);
		// cpx_set_ui (zeta, 0, 0);
		// cpx_polylog (zeta, ess, zee, prec);
		//
#if 0
		cpx_set (ms,ess);
		cpx_sub_ui (ms,ms,1, 0);
		cpx_mpf_pow (z2, que, ms, prec);
		cpx_mul (z2, z2, facup);

		mpf_ui_sub (que, 1, que);
		cpx_mpf_pow (z3, que, ms, prec);
		cpx_mul (z3, z3, facdown);
		
		// cpx_sub (z2, zeta, z2);
		// cpx_add (z3, z2, z3);
#endif

#if 0
		fp_two_pi (que, prec);
		cpx_mul_mpf (z3, ms, que);
		cpx_times_i (z3, z3);
		cpx_neg (z3, z3);
		cpx_exp (z3, z3, prec);
		cpx_sub_ui (z3,z3, 1, 0);
		cpx_mul (z2, z2, z3);
		cpx_add (z2, zeta, z2);
#endif


#if 0
		cpx_conj (ess, ess);
		mpf_ui_sub (que, 1, que);
		cpx_periodic_zeta (z2, ess, que, prec);
		cpx_conj (ess, ess);
#endif

// #define PARALLEL_TO_CUT 1
#ifdef PARALLEL_TO_CUT
		printf ("%g\t",zre);
		cpx_set_d (zee, zre, zim);
		cpx_polylog (zeta, ess, zee, prec);
		double zetare = mpf_get_d(zeta[0].re);
		double zetaim = mpf_get_d(zeta[0].im);
		printf ("%g\t%g\t", zetare, zetaim);

		polylog_sheet(z2, ess, zee, 1, prec);
		zetare = mpf_get_d(z2[0].re);
		zetaim = mpf_get_d(z2[0].im);
		printf ("%g\t%g\t", zetare, zetaim);

		polylog_sheet(z3, ess, zee, 2, prec);
		zetare = mpf_get_d(z3[0].re);
		zetaim = mpf_get_d(z3[0].im);
		printf ("%g\t%g\t", zetare, zetaim);
#endif

#define CROSS_CUT 1
#ifdef CROSS_CUT
		printf ("%g\t",zim);
		cpx_set_d (ess, sre, sim);
		cpx_set_d (zee, zre, zim);
#if 0
		// cpx_polylog (zeta, ess, zee, prec);
		cpx_set_ui (zeta, 0, 0);
		int m = 1;
		if (zim > 0)
		{
			cpx_polylog_sheet(z2, ess, zee, m, 1, prec);
			cpx_sub(zeta, zeta, z2);
		} else {
			/* N is less here */
			cpx_polylog_sheet(z2, ess, zee, m-1, -1, prec);
			cpx_sub(zeta, zeta, z2);
		}
#endif
#if 0
		// cpx_polylog (zeta, ess, zee, prec);
		cpx_set_ui (zeta, 0, 0);
		if (zim > 0)
		{
			cpx_polylog_sheet_g1_action(z2, ess, zee, 0, 1, prec);
			cpx_mul (z2,z2,ph);
			cpx_sub(zeta, zeta, z2);
		} else {
			/* N is less here */
			cpx_polylog_sheet_g1_action(z2, ess, zee, 0, -1, prec);
			cpx_sub(zeta, zeta, z2);
		}
#endif
#if 0
		// cpx_polylog (zeta, ess, zee, prec);
		cpx_set_ui (zeta, 0, 0);
		if (zim > 0)
		{
			cpx_polylog_sheet_g1_action(z2, ess, zee, 0, 1, prec);
			cpx_mul (z2,z2,ph);
			cpx_sub(zeta, zeta, z2);
			cpx_polylog_sheet_g1_action(z2, ess, zee, 0, -1, prec);
			cpx_sub(zeta, zeta, z2);
		} else {
			/* below the cut -- this is the two g_1's down */
			cpx_polylog_sheet_g1_action(z2, ess, zee, 0, -1, prec);
			cpx_sub(zeta, zeta, z2);
			cpx_polylog_sheet_g1_action(z2, ess, zee, -1, -1, prec);
			cpx_sub(zeta, zeta, z2);
		}
#endif
#if 1
		cpx_polylog (zeta, ess, zee, prec);
		cpx_set_ui (zeta, 0, 0);
		if (zim > 0)
		{
			cpx_polylog_sheet_g1_action(z2, ess, zee, 0, 1, prec);
			cpx_mul (z2,z2,ph);
			cpx_sub(zeta, zeta, z2);
			cpx_polylog_sheet_g1_action(z2, ess, zee, 0, -1, prec);
			cpx_sub(zeta, zeta, z2);
			cpx_polylog_sheet_g1_action(z2, ess, zee, -1, -1, prec);
			cpx_sub(zeta, zeta, z2);
		} else {
			/* below the cut -- this is the two g_1's down */
			cpx_polylog_sheet_g1_action(z2, ess, zee, 0, -1, prec);
			cpx_sub(zeta, zeta, z2);
			cpx_polylog_sheet_g1_action(z2, ess, zee, -1, -1, prec);
			cpx_sub(zeta, zeta, z2);
			cpx_polylog_sheet_g1_action(z2, ess, zee, -2, -1, prec);
			cpx_sub(zeta, zeta, z2);
		}
#endif
#if 0
		cpx_polylog (zeta, ess, zee, prec);
		if (q < 0)
		{
			cpx_polylog_sheet(z2, ess, zee, -2, prec);
			cpx_sub(zeta, zeta, z2);
		}
		else
		{
			cpx_polylog_sheet(z2, ess, zee, -1, prec);
			cpx_sub(zeta, zeta, z2);
		}
#endif
		double zetare = mpf_get_d(zeta[0].re);
		double zetaim = mpf_get_d(zeta[0].im);
		printf ("%g\t%g\t", zetare, zetaim);
#endif

		// fp_prt ("\t", zeta[0].re);
		// fp_prt ("\t", z2[0].re);
		// fp_prt ("\t", z3[0].re);
		// fp_prt ("\t", zeta[0].im);
		printf ("\n");
		fflush (stdout);
	}
#endif 

#ifdef COHERENT_BORKEN
	/* remnants of failed experiment */
	cpx_t lam;
	cpx_init (lam);
	cpx_set_d (lam, sim, 0.0);
	cpx_set_d (lam, 0.5, 14.134725);
	cpx_set_d (lam, 0.5, sim);
	
	double z;
	for (z=-0.6; z< 0.6; z+= 0.03)
	{
		cpx_set_d (zee, z, 0);
		cohere (zeta, lam, zee, prec);

		printf ("%g",z);
		fp_prt ("\t", zeta[0].re);
		// fp_prt ("\t", zeta[0].im);
		printf ("\n");
		fflush (stdout);
	}
#endif 
	return 0;
}
