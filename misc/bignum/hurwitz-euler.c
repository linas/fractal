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

void zeta_euler(cpx_t zeta, cpx_t ess, mpf_t q, int prec, int em, int pee);

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
	
	mpf_t que;
	mpf_init (que);
	mpf_set_d (que, 0.2);
	
	mpf_set (cq[0].re, que);
	mpf_set_ui (cq[0].im, 0);
	
	// cpx_hurwitz_zeta (zeta, ess, que, prec);
	cpx_hurwitz_taylor (zeta, ess, cq, prec);
	cpx_sub (ess, zeta, prevzeta);

	printf ("prec=%d ", prec);

	long rex, imx;
	mpf_get_d_2exp (&rex, ess[0].re);
	mpf_get_d_2exp (&imx, ess[0].im);
	rex = -0.30103 *rex;
	imx = -0.30103 *imx;
	if (imx < rex) rex = imx;
	printf ("prev=%ld ", rex);

	double vre = mpf_get_d (ess[0].re);
	printf ("vre = %g ", vre);

#if 0
	double vim = mpf_get_d (ess[0].im);
	printf ("vim = %g ", vim);
#endif

#if 0
	// gmp_printf ("re= %.220Ff ", zeta[0].re);
	printf ("re= ");
	mpf_out_str (stdout, 10, prec, zeta[0].re);
	
	printf (" im= ");
	mpf_out_str (stdout, 10, 120, zeta[0].im);
#endif
	printf ("\n");
	fflush (stdout);

	cpx_clear (ess);
	cpx_clear (cq);
	cpx_clear (prevzeta);
}

double err_est (cpx_t ess, int em, int pee)
{
	/* (2pi)^{-2p} ((s+2p)! / (s-1)!)  (M+a)^{-2p}*/
	cpx_t poch;
	cpx_init (poch);

	cpx_poch_rising (poch, ess, 2*pee);

	mpf_t pabs, den;
	mpf_init (pabs);
	mpf_init (den);

	cpx_abs (pabs, poch);

	double fden = 2*3.1416*em;
	mpf_set_d (den, fden);
	fp_log (den, den, 40);
	mpf_mul_ui (den,den, 2*pee);
	mpf_neg (den,den);
	fp_exp(den,den,40);
	
	mpf_mul (pabs,pabs, den);

	double po = mpf_get_d (pabs);

	cpx_clear (poch);
	mpf_clear (pabs);
	mpf_clear (den);
	return po;
}

int main ()
{
	int prec = 40;
	int i;

	/* Set the precision (number of binary bits) */
	mpf_set_default_prec (3.3*prec+140);
	mpf_set_default_prec (1000);

	cpx_t ess, zeta, z2, z3;
	cpx_init (ess);
	cpx_init (zeta);
	cpx_init (z2);
	cpx_init (z3);

	cpx_set_d (ess, 0.5, 4.0);
	cpx_set_d (ess, 2.0, 0.1);
	cpx_set_d (ess, 0.5, 14.134725);
	cpx_set_d (ess, 3.0, 0.1);
	cpx_set_d (ess, 0.5, 4.0);
	
	mpf_t que;
	mpf_init (que);
	mpf_set_d (que, 0.2);
	
	cpx_t cq;
	cpx_init (cq);
	mpf_set (cq[0].re, que);
	mpf_set_ui (cq[0].im, 0);
	
#if 0
	for (i=20; i < 10000; i*=1.5)
	// for (i=100; i > 10; i/=1.3)
	{
		hiprec (zeta, i);
	}
#endif
	
#if 1
	cpx_hurwitz_zeta (z2, ess, que, prec);
	fp_prt ("poly   ", z2[0].re);
	printf ("\n");

	cpx_hurwitz_taylor (z3, ess, cq, prec);
	fp_prt ("taylor ", z3[0].re);
	printf ("\n");

	cpx_hurwitz_euler (zeta, ess, que, prec);
	fp_prt ("euler  ", zeta[0].re);
	printf ("\n");
#endif

#if 1
	int pee, em;
	for (i=20; i<900; i+=40)
	{
		pee = i;
		em = 0.33*pee +12;
		double err = err_est (ess, em, pee);
		printf ("err=%g\n", err);
		zeta_euler (zeta, ess, que, prec, em, pee);
		printf ("p=%d m=%d ", pee, em);
		fp_prt ("its ", zeta[0].re);
		printf ("\n");
	}
#endif
	// do_perf (prec);
	return 0;
}

