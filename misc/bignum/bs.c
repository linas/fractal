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

double find_zero (double lo, double hi,  double f(double), double prec)
{
	double flo = f(lo);
	double fhi = f(hi);
	if (flo*fhi > 0.0)
	{
		printf ("error flo=%g  fhi=%g\n", flo, fhi);
		return 0.0;
	}
	while (hi-lo > prec)
	{
		double mid = lo - flo * (hi-lo)/ (fhi-flo);

		if (mid <= lo | mid >= hi)
		{
			mid = 0.5*(lo+hi);
		}
		double fmid = f(mid);
printf ("duude %14.12g %14.12g %14.12g hav %14.12g  del %14.12g\n", lo, mid, hi, fmid,
hi-lo);
		if (flo*fmid < 0.0)
		{
			hi = mid;
			fhi = fmid;
		}
		else
		{
			lo = mid;
			flo = fmid;
		}
	}
	double mid = lo - flo * (hi-lo)/ (fhi-flo);
	return mid;
}

double eff(double x)
{
	mpf_t re_b, im_b;
	mpf_init (re_b);
	mpf_init (im_b);

	b_sub_s (re_b, im_b, x, 0.0, 100);
	double y = mpf_get_d (re_b);

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
	int nterms = 50;
	int bits = (int) (v + 300 + 3*nterms);
	
	/* set the precision (number of binary bits) */
	mpf_set_default_prec (bits);
	
	mpf_t re_a, im_a;
	mpf_init (re_a);
	mpf_init (im_a);

	int i;
	for (i=1;i<10; i++)
	{
		double bn = eff (i);
		printf ("duude %d  is %g\n", i, bn);
	}

	// find_zero (0.0, 5.0, eff, 1.0e-6);
	fflush (stdout);

}

