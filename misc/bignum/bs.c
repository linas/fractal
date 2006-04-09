/*
 * bs.c
 *
 * High-precison asub_n  using the 
 * Gnu Multiple-precision library.
 *
 * Also, high-precision values of the series a_n 
 * 
 * Linas Vepstas July 2005
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mp_zeta.h"

/* return zero by binary subdivision */

double find_zero (double low, double hi,  double f(double), double prec)
{
	double flow = f(low);
	double fhi = f(hi);
	while (hi-low > prec)
	{
		double mid = low - flow * (hi-low)/ (fhi-flow);
		double fmid = f(mid);
		if (flo*fmid < 0.0)
		{
			hi = mid;
			fhi = fmid;
		}
		else
		{
			lo = mid;
			flo = mid;
		}
	}
	double mid = low - flow * (hi-low)/ (fhi-flow);
	return mid;
}

double eff(double x)
{
	mpf_t re_b, im_b;
	mpf_init (re_b);
	mpf_init (im_b);

	b_sub_s (re_b, im_b, x, 0.0, 100);
	double y = mpf_get_d (re_b);
	double z = mpf_get_d (im_b);
printf ("duuude its %g %g\n", y,z);

	mpf_clear (re_b);
	mpf_clear (im_b);

	return y;
}
/* ==================================================================== */

main (int argc, char * argv[])
{
	char str[4000];

	if (argc < 3)
	{
		fprintf (stderr, "Usage: %s [ndigits] [nterms]\n", argv[0]); 
		exit (1);
	}

	/* the decimal precison (number of decimal places) */
	int prec = atoi (argv[1]);

	/* compute number of binary bits this corresponds to. */
	double v = ((double) prec) *log(10.0) / log(2.0);

	/* the variable-precision calculations are touchy about this */
	/* XXX this should be stirling's approx for binomial */ 
	int bits = (int) (v + 300 + 3*nterms);
	
	/* set the precision (number of binary bits) */
	mpf_set_default_prec (bits);
	
	mpf_t re_a, im_a;
	mpf_init (re_a);
	mpf_init (im_a);


	find_zero (1.0, 2.0, eff, 1.0e-6);
	fflush (stdout);

}

