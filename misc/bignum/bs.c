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

	/* precision that we should use is about 
	 * x*log(10/log(2) = 3.3*x since b(x) 
	 * rquires 2^x digits to calculate accurately.
	 *
	 * For real x, we need ... ?
	 */
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
	int prec = 1610;

	/* compute number of binary bits this corresponds to. */
	double v = ((double) prec) *log(10.0) / log(2.0);

	/* The largest that a binomial (n,k) will get is 2^n
	 * so need an extra norder bits if going to order norder. 
	 * And pad a bit, just to be safe... */
	int norder = 5000;
	int bits = (int) (v + 100 + norder);
	
	/* set the precision (number of binary bits) */
	mpf_set_default_prec (bits);
	
	printf ("# looking for precise zero locations\n");
	printf ("# computed to precision of %d decimal places\n", prec);
	printf ("# computed up to order of %d \n", norder);
	printf ("# computed with %d bits of default mpf \n", bits);
	fflush (stdout);
	
	double z[150];

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
	z[43] = 1528.62338;
	z[44] = 1598.720279;
	z[45] = 1670.388035;
	z[46] = 1743.626495;
	z[47] = 1818.43581;
	z[48] = 1894.81584;
	z[49] = 1972.766716;
	z[50] = 2052.288429;
	z[51] = 2133.38088;
	z[52] = 2216.044101;
	z[53] = 2300.278171;
	z[54] = 2386.082985;
	z[55] = 2473.458622;
	z[56] = 2562.405047;
	z[57] = 2652.922223;
	z[58] = 2745.010239;
	z[59] = 2838.669045;
	z[60] = 2933.898631;
	z[61] = 3030.69903;
	z[62] = 3129.070226;
	z[63] = 3229.012193;
	z[64] = 3330.52498;
	z[65] = 3433.608541;
	z[66] = 3538.262921;
	z[67] = 3644.488069;
	z[68] = 3752.284028;
	z[69] = 3861.650754;
	z[70] = 3972.5883;
	z[71] = 4085.09664;
	z[72] = 4199.175777;
	z[73] = 4314.82568;
	z[74] = 4432.046412;
	z[75] = 4550.837918;
	z[76] = 4671.200252;
	z[77] = 4793.133354;
							 
																									
	int i;
	for (i=36; i<78; i++)
	{
		double zer = find_zero (z[i]-0.02, z[i]+0.02, eff, 1.0e-9);

		printf ("%d\t%22.18g\n", i, zer);
		fflush (stdout);
	}

}

