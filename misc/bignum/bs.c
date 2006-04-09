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

		double b = 0.8*(hi-lo);
		if (mid <= lo+b | mid >= hi-b)
		{
			mid = 0.5*(lo+hi);
		}
		double fmid = f(mid);
// printf ("duude %14.12g %14.12g %14.12g hav %14.12g  del %14.12g\n", lo, mid, hi, fmid, hi-lo);
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

	if (argc < 1)
	{
		fprintf (stderr, "Usage: %s \n", argv[0]); 
		exit (1);
	}

	/* the decimal precison (number of decimal places) */
	int prec = 100;

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

	double z[50];
	z[1] = 2.119811771;
	z[2] = 6.160191799;
	z[3] = 11.77028675;
	z[4] = 19.0188987;
	z[5] = 27.89007811;
	z[6] = 38.28104304;
	z[7] = 50.26296985;
	z[8] = 63.82158993;
	z[9] = 78.94555959;
	z[10] = 95.63052706;
	z[11] = 113.8962431;
	z[12] = 133.72621;
	z[13] = 155.1275029;
	z[14] = 178.1015976;
	z[15] = 202.6474945;
	z[16] = 228.7634104;
	z[17] = 256.4484098;
	z[18] = 285.7058244;
	z[19] = 316.5328576;
	z[20] = 348.9326376;
	z[21] = 382.9013028;
	z[22] = 418.4397673;
	z[23] = 455.5503812;
	z[24] = 494.2314509;
	z[25] = 534.483618;
	z[26] = 576.3063385;
	z[27] = 619.7002344;
	z[28] = 664.664512;
	z[29] = 711.1993656;
	z[30] = 759.3053082;
	z[31] = 808.9828791;
	z[32] = 860.2294413;
	z[33] = 913.0476106;
	z[34] = 967.4367751;
	z[35] = 1023.39658;
	z[36] = 1080.927689;
	z[37] = 1140.028462;
	z[38] = 1200.700888;
	z[39] = 1262.944173;
	z[40] = 1326.75764;
	z[41] = 1392.141956;
	z[42] = 1459.097276;

	int i;
	for (i=1; i<43; i++)
	{
		double zer = find_zero (z[i]-1.0, z[i]+1.0, eff, 1.0e-6);

		printf ("duude %d gu=%g act=%g\n", i, z[i], zer);
		fflush (stdout);
	}

}

