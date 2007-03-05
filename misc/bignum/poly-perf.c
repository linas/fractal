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
char * restr = "0.1";
char * imstr = "0.1";
/* ================================================================= */

void do_perf(void)
{
	int i;
	int prec = 40;

	printf ("#\n# graph of periodic zeta as function of precision \n#\n");
	int hz = sysconf (_SC_CLK_TCK);
	printf ("# clock ticks=%d\n#\n", hz);
	fflush (stdout);

	for (prec=10; prec <123123; prec *= 1.41)
	{
		/* Set the precision (number of binary bits) */
		int nbits = 3.3*prec+100;
		mpf_set_default_prec (nbits);

		printf ("%d\t%d\t", prec, nbits);

		cpx_t ess, zeta, zee, plog, expected;
		cpx_init (ess);
		cpx_init (zeta);
		cpx_init (zee);
		cpx_init (plog);
		cpx_init (expected);
	
		mpf_t que;
		mpf_init (que);

		mpf_set_str (expected[0].re, restr, 10);
		mpf_set_str (expected[0].im, imstr, 10);
	
		struct tms start, end;

#define MEASURE_POLYLOG_PERFORMANCE 1
#ifdef MEASURE_POLYLOG_PERFORMANCE
		cpx_set_d (ess, 0.5, 14.134725);
		cpx_set_d (zee, 0.4, 0.3);

#if 1
		/* First we warm the cache */
		times (&start);
		cpx_polylog (zeta, ess, zee, prec);
		times (&end);

		cpx_sub (zeta, zeta, expected);
		int rex = get_prec (zeta, prec);
		printf ("%d\t", rex);

		printf ("%jd\t", (intmax_t) (end.tms_utime - start.tms_utime));
		
		/* Then with a hot cache */
		times (&start);
		for (i=0; i<1000; i++)
			cpx_polylog (zeta, ess, zee, prec);
		times (&end);
		printf ("%jd\t", (intmax_t) (end.tms_utime - start.tms_utime));

		cpx_sub (zeta, zeta, expected);
		rex = get_prec (zeta, prec);
		printf ("%d\t", rex);
#endif

#if 0
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
		fflush (stdout);
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
		cpx_clear (expected);
		mpf_clear (que);

		fflush (stdout);
	}
}

/* ================================================================= */


/* Compute to high precision, monitor the precision */
void hiprec(cpx_t zeta, int prec)
{
	/* Set the precision (number of binary bits) */
	int nbits = 3.322*prec+2*3.14*3.32*prec+140;
	// nbits = 3.322*prec+ 7*prec + 140;
	nbits = 3.322*prec+ 140;
	mpf_set_default_prec (nbits);
	cpx_set_prec (zeta, nbits);

	cpx_t cq, ess, prevzeta;
	cpx_init (ess);
	cpx_init (cq);
	cpx_init (prevzeta);
	cpx_set (prevzeta, zeta);

	cpx_set_d (ess, 0.5, 14.134725);
	cpx_set_d (cq, 0.4, 0.3);
	
	cpx_polylog_euler (zeta, ess, cq, prec);
	cpx_sub (ess, zeta, prevzeta);

	printf ("prec=%d ", prec);

	long rex = get_prec (ess, prec);
	printf ("prev=%ld ", rex);

	double vre = mpf_get_d (ess[0].re);
	printf ("vre = %g ", vre);

#if 0
	double vim = mpf_get_d (ess[0].im);
	printf ("vim = %g ", vim);
#endif

#if 1
	// gmp_printf ("re= %.220Ff ", zeta[0].re);
	printf ("re= ");
	mpf_out_str (stdout, 10, prec, zeta[0].re);
	
	printf (" im= ");
	mpf_out_str (stdout, 10, prec, zeta[0].im);
#endif
	printf ("\n");
	fflush (stdout);

	cpx_clear (ess);
	cpx_clear (cq);
	cpx_clear (prevzeta);
}

/* ================================================================= */

int
main (int argc, char * argv[])
{
	int prec = 40;

#if 0
	if (argc != 2)
	{
		fprintf (stderr, "Usage: %s <sim>\n", argv[0]);
		exit (1);
	}
	double sim = atof (argv[1]);
#endif

	cpx_t zeta;
	cpx_init (zeta);
	for (prec=10; prec<15001; prec *= 1.4)
	{
		hiprec(zeta, prec);
	}
	
	// do_perf();
	
	return 0;
}
