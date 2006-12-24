/**
 * exact.c
 *
 * Look at hurtwitz zeta at the location of the Riemann zeta zeros
 *
 * December 2006
 */

#include <stdio.h>
#include <stdlib.h>

#include "mp-misc.h"
#include "mp-polylog.h"

int
main ()
{
	int prec = 30;
	double q;

	/* Set the precision (number of binary bits) */
	int nbits = 3.3*prec+100;
	mpf_set_default_prec (nbits);

	cpx_t ess, zeta;
	cpx_init (ess);
	cpx_init (zeta);

	mpf_t que;
	mpf_init (que);
			  
	cpx_set_d (ess, 0.5, 14.1);
	char * zero = "14.134725141734693790457251983562470270784257115699243175685567460149 \
	               9634298092567649490103931715610127792029715487974367661426914698822545 \
	               8250536323944713778041338123720597054962195586586020055556672583601077";
	mpf_set_str (ess[0].im, zero, 10);

	for (q=0.1; q<1.0; q++)
	{
		mpf_set_d (que, q);
		cpx_hurwitz_zeta (zeta, ess, que, prec);

		printf ("%g",q);
		fp_prt ("\t", zeta[0].re);
		fp_prt ("\t", zeta[0].im);
		printf ("\n");
		fflush (stdout);
	}

	return 0;
}
