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
	
	cpx_t ess, zeta, zee, plog, gm;
	cpx_init (ess);
	cpx_init (zeta);
	cpx_init (zee);
	cpx_init (plog);
	cpx_init (gm);

	mpf_t que, tp;
	mpf_init (que);
	mpf_init (tp);
			  
	cpx_set_d (ess, sre, sim);

	cpx_add_ui (ess, ess, 1, 0);
	cpx_gamma (gm, ess, prec);
	cpx_times_ui (gm, gm, 2);
	cpx_recip (gm, gm);
	cpx_sub_ui (ess, ess, 1, 0);

	mpf_t twopi;
	mpf_init (twopi);
	fp_two_pi (twopi, prec);

	printf ("#\n# graph of periodic zeta as function of q, at \n#\n");
	printf ("# at s=%g +i %g\n", sre, sim);
	printf ("#\n# prec=%d nbits=%d\n#\n", prec, nbits);
	fflush (stdout);
	for (q=0.001; q<0.999; q+=0.008)
	{
		printf ("%g\t", q);
		mpf_set_d (que, q);

		// cpx_hurwitz_zeta (zeta, ess, que, prec);
		// cpx_periodic_zeta (zeta, ess, que, prec);

		cpx_periodic_beta (zeta, ess, que, prec);

		mpf_set_d (que, 1.0-q);
		cpx_periodic_beta (zee, ess, que, prec);
		cpx_add (zeta, zeta, zee);

		double zetare = mpf_get_d(zeta[0].re);
		// double zetaim = mpf_get_d(zeta[0].im);
		// printf ("%g\t%g\t", zetare, zetaim);
		printf ("%g\t", zetare);

		printf ("\n");
		fflush (stdout);
	}

	return 0;
}
