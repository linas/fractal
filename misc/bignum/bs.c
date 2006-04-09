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

	b_sub_s (re_b, im_b, x, 0.0, 400, x*4+150);
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
	int nterms = 1500;
	int bits = (int) (v + 400 + 3*nterms);
	
	/* set the precision (number of binary bits) */
	mpf_set_default_prec (bits);
	
	double z[50];

	z[1] = 2.756189414;
	z[2] = 6.968997419;
	z[3] = 12.72989432;
	z[4] = 20.01889832;
	z[5] = 28.86402965;
	z[6] = 39.28042409;
	z[7] = 51.26259237;
	z[8] = 64.812872;
	z[9] = 79.93504911;
	z[10] = 96.62793513;
	z[11] = 114.8901331;
	z[12] = 134.7234457;
	z[13] = 156.1274901;
	z[14] = 179.101592;
	z[15] = 203.6462289;
	z[16] = 229.7615883;
	z[17] = 257.448081;
	z[18] = 286.7046825;
	z[19] = 317.5324154;
	z[20] = 349.9305071;
	z[21] = 383.8995559;
	z[22] = 419.4395815;
	z[23] = 456.550048;
	z[24] = 495.2314281;
	z[25] = 535.4834264;
	z[26] = 577.3062935;
	z[27] = 620.699736;
	z[28] = 665.6641158;
	z[29] = 712.1993556;
	z[30] = 760.3052747;
	z[31] = 809.9818344;
	z[32] = 861.2294288;
	z[33] = 914.0476105;
	z[34] = 968.4366987;
	z[35] = 1024.396526;
	z[36] = 1081.927037;
	z[37] = 1141.028462;
	z[38] = 1201.700635;
	z[39] = 1263.943587;
	z[40] = 1327.757352;
	z[41] = 1393.141954;
	z[42] = 1460.097275;
	int i;
	for (i=38; i<43; i++)
	{
		double zer = find_zero (z[i]-0.002, z[i]+0.002, eff, 1.0e-9);

		printf ("%d\t%20.16g\n", i, zer);
		fflush (stdout);
	}

}

