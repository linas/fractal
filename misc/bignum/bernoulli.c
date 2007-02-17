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

	if (argc != 4)
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

	mpf_t que, tp;
	mpf_init (que);
	mpf_init (tp);
			  
	cpx_set_d (ess, sre, sim);

	mpf_t twopi;
	mpf_init (twopi);
	fp_two_pi (twopi, prec);

	printf ("#\n# graph of periodic zeta as function of q, at \n#\n");
	printf ("# at s=%g +i %g\n", sre, sim);
	printf ("#\n# prec=%d nbits=%d\n#\n", prec, nbits);
	fflush (stdout);
	for (q=0.001; q<0.999; q+=0.008)
	{
		mpf_set_d (que, q);

		// cpx_hurwitz_zeta (zeta, ess, que, prec);
		// cpx_periodic_zeta (zeta, ess, que, prec);

		cpx_periodic_beta (zeta, ess, que, prec);
		double zetare = mpf_get_d(zeta[0].re);
		double zetaim = mpf_get_d(zeta[0].im);
		printf ("%g\t%g\t", zetare, zetaim);

		printf ("\n");
		fflush (stdout);
	}

	return 0;
}
