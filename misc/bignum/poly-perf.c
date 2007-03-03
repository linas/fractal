/**
 * poly-perf.c
 *
 * Make comparative performance measurements for
 * polylog algorithms.
 *
 * January 2007
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/times.h>
#include <stdint.h>
#include <unistd.h>

#include "mp-complex.h"
#include "mp-consts.h"
#include "mp-gamma.h"
#include "mp-misc.h"
#include "mp-polylog.h"
#include "mp-trig.h"

/* ================================================================= */

int
main (int argc, char * argv[])
{
	int i;
	int prec = 40;
	// prec = 20;

#if 0
	if (argc != 2)
	{
		fprintf (stderr, "Usage: %s <sim>\n", argv[0]);
		exit (1);
	}
	double sim = atof (argv[1]);
#endif
	
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

#define MEASURE_POLYLOG_PERFORMANCE 1
#ifdef MEASURE_POLYLOG_PERFORMANCE
		cpx_set_d (ess, 0.5, 14.134725);
		cpx_set_d (zee, 0.4, 0.3);

#if 0
		/* First we warm the cache */
		times (&start);
		cpx_polylog (zeta, ess, zee, prec);
		times (&end);
		printf ("%jd\t", (intmax_t) (end.tms_utime - start.tms_utime));
		
		/* Then with a hot cache */
		times (&start);
		for (i=0; i<1000; i++)
			cpx_polylog (zeta, ess, zee, prec);
		times (&end);
		printf ("%jd\t", (intmax_t) (end.tms_utime - start.tms_utime));
#endif

#if 1
		/* First we warm the cache */
		times (&start);
		cpx_polylog_euler (zeta, ess, zee, prec);
		times (&end);
		printf ("%jd\t", (intmax_t) (end.tms_utime - start.tms_utime));
		
		/* Then with a hot cache */
		times (&start);
		for (i=0; i<1000; i++)
			cpx_polylog_euler (zeta, ess, zee, prec);
		times (&end);
		printf ("%jd\t", (intmax_t) (end.tms_utime - start.tms_utime));
#endif

#if 0
		/* First we warm the cache */
		times (&start);
		cpx_polylog_sum (plog, ess, zee, prec);
		times (&end);
		printf ("%jd\t", (intmax_t) (end.tms_utime - start.tms_utime));

		/* Then with a hot cache */
		times (&start);
		for (i=0; i<100; i++)
			cpx_polylog_sum (plog, ess, zee, prec);
		times (&end);
		printf ("%jd\t", (intmax_t) (end.tms_utime - start.tms_utime));
#endif

		printf ("\n");
#endif

// #define MEASURE_HURWITZ_PERFORMANCE
#ifdef MEASURE_HURWITZ_PERFORMANCE
		cpx_set_d (ess, 0.5, 14.134725);
		cpx_set_d (zee, 0.2, 0.0);
		mpf_set_d (que, 0.2);

#if 0
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
		printf ("\n");
#endif

		cpx_clear (ess);
		cpx_clear (zeta);
		cpx_clear (zee);
		cpx_clear (plog);
		mpf_clear (que);

		fflush (stdout);
	}
	return 0;
}
