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
#include "mp-misc.h"
#include "mp-polylog.h"


#ifdef BORKOKEN_AND_DONT_WORK

#include "mp-consts.h"
#include "mp-gamma.h"
#include "mp-trig.h"
#include "mp-zeta.h"

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

int
main (int argc, char * argv[])
{
	int prec = 20;
	double q;

	/* Set the precision (number of binary bits) */
	int nbits = 3.3*prec+100;
	mpf_set_default_prec (nbits);

	if (argc != 2)
	{
		fprintf (stderr, "Usage: %s <sim>\n", argv[0]);
		exit (1);
	}
	double sim = atof (argv[1]);
	
	cpx_t ess, zeta, zee, plog;
	cpx_init (ess);
	cpx_init (zeta);
	cpx_init (zee);
	cpx_init (plog);

	mpf_t que;
	mpf_init (que);
			  
	cpx_set_d (ess, 0.5, sim);
	cpx_set_d (ess, 1.5, 14.134725);

	double zmag = sim;

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

	zero = "50.0";
	mpf_set_str (ess[0].im, zero, 10);
#endif

	// printf ("#\n# graph of Hurwitz zeta as function of q, at \n#\n");
	printf ("#\n# graph of periodic zeta as function of q, at \n#\n");
	fp_prt ("# at s=0.5+i ", ess[0].im);
	printf ("\n#\n# prec=%d nbits=%d\n#\n", prec, nbits);
	fflush (stdout);
	for (q=0.02; q<0.991; q+=0.008)
	{
		mpf_set_d (que, q);
		// cpx_hurwitz_zeta (zeta, ess, que, prec);
		// cpx_periodic_beta (zeta, ess, que, prec);
		//
	cpx_set_d (ess, 1.5, 14.134725);
		cpx_periodic_zeta (zeta, ess, que, prec);

		double zre = zmag * cos (2.0*M_PI * q);
		double zim = zmag * sin (2.0*M_PI * q);
		cpx_set_d (zee, zre, zim);
	cpx_set_d (ess, 0.5, 14.134725);
		cpx_polylog_sum (plog, ess, zee, prec);

		printf ("%g",q);
		fp_prt ("\t", zeta[0].re);
		// fp_prt ("\t", zeta[0].im);
		fp_prt ("\t", plog[0].re);
		// fp_prt ("\t", plog[0].im);
		printf ("\n");
		fflush (stdout);
	}

	return 0;
}
