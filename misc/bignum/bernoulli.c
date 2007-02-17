/**
 * bernoulli.c
 *
 * Graphs of the eigenvalues of the Bernoulli map.
 *
 * February 2007
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

/* ============================================================================ */

int
main (int argc, char * argv[])
{
	int prec = 40;
	prec = 20;
	int i;
	double q;

	/* Set the precision (number of binary bits) */
	int nbits = 3.3*prec+100;
	mpf_set_default_prec (nbits);

	if (argc != 3)
	{
		fprintf (stderr, "Usage: %s <sre> <sim>\n", argv[0]);
		exit (1);
	}
	double sre = atof (argv[1]);
	double sim = atof (argv[2]);
	
	cpx_t ess, zeta, zee, plog;
	cpx_init (ess);
	cpx_init (zeta);
	cpx_init (zee);
	cpx_init (plog);

	mpf_t que, l2;
	mpf_init (que);
	mpf_init (l2);
	fp_log2 (l2, prec);
			  
	mpf_t twopi;
	mpf_init (twopi);
	fp_two_pi (twopi, prec);
	mpf_div (twopi, twopi, l2);

	/* compute 1/2 Gamma (s+1+2\pi n) */
	cpx_t gm[10];
	mpf_t agm[10];
	cpx_set_d (ess, sre+1.0, sim);
	for (i=0; i<8; i++)
	{
		cpx_init (gm[i]);
		mpf_init (agm[i]);

		cpx_gamma (gm[i], ess, prec);
		cpx_recip (gm[i], gm[i]);
		cpx_abs (agm[i], gm[i]);

		mpf_add (ess[0].im, ess[0].im, twopi);		
	}


	printf ("#\n# graph of periodic zeta as function of q, at \n#\n");
	printf ("# at s=%g +i %g\n", sre, sim);
	printf ("#\n# prec=%d nbits=%d\n#\n", prec, nbits);
	fflush (stdout);
	for (q=0.001; q<0.999; q+=0.003)
	{
		printf ("%g\t", q);

		// cpx_hurwitz_zeta (zeta, ess, que, prec);
		// cpx_periodic_zeta (zeta, ess, que, prec);

		cpx_set_d (ess, sre, sim);
		for (i=0; i<6;i++)
		{
			mpf_set_d (que, q);
			cpx_periodic_beta (zeta, ess, que, prec);

			// mpf_set_d (que, 1.0-q);
			// cpx_periodic_beta (zee, ess, que, prec);
			// cpx_add (zeta, zeta, zee);

			// cpx_times_mpf (zeta, zeta, agm[i]);
			cpx_mul (zeta, zeta, gm[i]);
			cpx_div_ui (zeta, zeta, 2);

			double zetare = mpf_get_d(zeta[0].re);
			double zetaim = mpf_get_d(zeta[0].im);
			double mag = sqrt (zetare*zetare + zetaim*zetaim);
			double arg = atan2 (zetaim, zetare);
			// printf ("%g\t%g\t", zetare, zetaim);
			// printf ("%g\t", zetare);
			// printf ("%g\t", zetaim);
			// printf ("%g\t", mag);
			printf ("%g\t", arg);

			mpf_add (ess[0].im, ess[0].im, twopi);		
		}

		printf ("\n");
		fflush (stdout);
	}

	return 0;
}
