/* 
 * hurwitz-euler.c
 *
 * Compute the Hurwitz zeta function for arbitrary complex argument
 * Use the Euler-Maclaurin formla to do so.
 *
 * Linas Vepstas Feb 2007
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/times.h>
#include <stdint.h>
#include <unistd.h>


#include <gmp.h>
#include "mp-binomial.h"
#include "mp-cache.h"
#include "mp-complex.h"
#include "mp-misc.h"
#include "mp-polylog.h"
#include "mp-trig.h"
#include "mp-zeta.h"

/* =========================================================== */
	
void do_perf (int prec)
{
	int i;

	printf ("#\n# graph of periodic zeta as function of precision \n#\n");
	int hz = sysconf (_SC_CLK_TCK);
	printf ("# clock ticks=%d\n#\n", hz);
	fflush (stdout);

	for (prec=20; prec <123123; prec *= 1.41)
	{
		/* Set the precision (number of binary bits) */
		int nbits = 3.3*prec+100;
		mpf_set_default_prec (nbits);

		printf ("%d\t%d\t", prec, nbits);

		cpx_t ess, zeta, zee, plog;
		cpx_init (ess);
		cpx_init (zeta);
		cpx_init (zee);
		cpx_init (plog);
	
		mpf_t que;
		mpf_init (que);

		struct tms start, end;

		cpx_set_d (ess, 0.5, 14.134725);
		cpx_set_d (zee, 0.2, 0.0);
		mpf_set_d (que, 0.2);

#if 1
		/* First we warm the cache */
		times (&start);
		cpx_hurwitz_zeta (zeta, ess, que, prec);
		times (&end);
		printf ("%jd\t", (intmax_t) (end.tms_utime - start.tms_utime));

		/* Then with a hot cache */
		times (&start);
		for (i=0; i<1000; i++)
			cpx_hurwitz_zeta (zeta, ess, que, prec);
		times (&end);
		printf ("%jd\t", (intmax_t) (end.tms_utime - start.tms_utime));
#endif

#if 1
		/* First we warm the cache */
		times (&start);
		cpx_hurwitz_taylor (zeta, ess, zee, prec);
		times (&end);
		printf ("%jd\t", (intmax_t) (end.tms_utime - start.tms_utime));

		/* Then with a hot cache */
		times (&start);
		for (i=0; i<1000; i++)
			cpx_hurwitz_taylor (zeta, ess, zee, prec);
		times (&end);
		printf ("%jd\t", (intmax_t) (end.tms_utime - start.tms_utime));
#endif
		
#if 1
		/* First we warm the cache */
		times (&start);
		cpx_hurwitz_euler (zeta, ess, que, prec);
		times (&end);
		printf ("%jd\t", (intmax_t) (end.tms_utime - start.tms_utime));

		/* Then with a hot cache */
		times (&start);
		for (i=0; i<1000; i++)
			cpx_hurwitz_euler (zeta, ess, que, prec);
		times (&end);
		printf ("%jd\t", (intmax_t) (end.tms_utime - start.tms_utime));
#endif
		printf ("\n");

		cpx_clear (ess);
		cpx_clear (zeta);
		cpx_clear (zee);
		cpx_clear (plog);
		mpf_clear (que);

		fflush (stdout);
	}
}

int main ()
{
	int prec = 40;

	/* Set the precision (number of binary bits) */
	mpf_set_default_prec (3.3*prec+140);

	cpx_t ess, zeta, z2, z3;
	cpx_init (ess);
	cpx_init (zeta);
	cpx_init (z2);
	cpx_init (z3);

	cpx_set_d (ess, 0.5, 4.0);
	cpx_set_d (ess, 2.0, 0.1);
	
	mpf_t que;
	mpf_init (que);
	mpf_set_d (que, 0.555);
	
	cpx_t cq;
	cpx_init (cq);
	mpf_set (cq[0].re, que);
	mpf_set_ui (cq[0].im, 0);
	
	cpx_hurwitz_zeta (z2, ess, que, prec);
	fp_prt ("poly   ", z2[0].re);
	printf ("\n");

	cpx_hurwitz_taylor (z3, ess, cq, prec);
	fp_prt ("taylor ", z3[0].re);
	printf ("\n");

	cpx_hurwitz_euler (zeta, ess, que, prec);
	fp_prt ("its    ", zeta[0].re);
	printf ("\n");

#if 0
	int pee, em;
	for (pee=20; pee<900; pee+=40)
	{
		em = 2*pee +12;
		zeta_euler (zeta, ess, que, prec, em, pee);
		printf ("p=%d m=%d ", pee, em);
		fp_prt ("its ", zeta[0].re);
		printf ("\n");
	}
#endif
	do_perf (prec);
	return 0;
}

