/*
 * bs.c
 *
 * Find the zeros of b_sub_n at high preicision, by using 
 * set of initial guesses and a simple root-finding algo.
 * 
 * Linas Vepstas April 2006
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mp_zeta.h"

#ifdef USE_FLT_PT_BINARY_SUBDIVISION
/* Find zero by binary subdivision */
double binary_find_zero (double lo, double hi,  
					 double f(double, void *), double prec)
{
	double flo = f(lo, NULL);
	double fhi = f(hi, NULL);
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
		double fmid = f(mid, NULL);
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
#endif /* USE_FLT_PT_BINARY_SUBDIVISION */

int prec;
int norder;

double eff(double x, void * params)
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
	b_sub_s (re_b, im_b, x, 0.0, prec, norder, -1.0e-30);
	double y = mpf_get_d (re_b);

	mpf_clear (re_b);
	mpf_clear (im_b);

	return y;
}

#define USE_GSL_SLOVER
#ifdef USE_GSL_SLOVER

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

double brent_solver (double x_lo_bound, double x_hi_bound, double dprec)
{
	const gsl_root_fsolver_type *T;
	T = gsl_root_fsolver_brent;

	gsl_root_fsolver *s;
	s = gsl_root_fsolver_alloc (T);

	gsl_function F;
	F.function = eff;
	F.params = NULL;

	gsl_root_fsolver_set (s, &F, x_lo_bound, x_hi_bound);
	// printf ("# using %s method\n", gsl_root_fsolver_name (s));

	double root = 0.0;
	int status;
	int iter=0;	
	do
	{
		iter ++;
		status = gsl_root_fsolver_iterate (s);
		root = gsl_root_fsolver_root (s);
		x_lo_bound = gsl_root_fsolver_x_lower (s);
		x_hi_bound = gsl_root_fsolver_x_upper (s);
		status = gsl_root_test_interval (x_lo_bound, x_hi_bound, 0, dprec);
	}
	while (status == GSL_CONTINUE && iter < 30);

	return root;
}
#endif

// #define BIGNUM_ROOT_FINDER
#ifdef BIGNUM_ROOT_FINDER

void func (mpf_t y, mpf_t x, void * params)
{
	mpf_t im_b;
	mpf_init (im_b);

	b_sub_s (y, im_b, x, 0.0, prec, norder, -1.0e-30);
	double y = mpf_get_d (re_b);

	mpf_clear (im_b);
}

/*
 * @prec: find root to this many decimal places
 */
void find_zero (mpf_t root, double root_bound_lo, double root_bound_hi,  
					 void f(mpf_t, mpf_t, void *), void *params, int prec)
{
	mpf_t ra, rb, rc, rd, re;
	mpf_init (ra);
	mpf_init (rb);
	mpf_init (rc);
	mpf_init (rd);
	mpf_init (re);

	/* Initialize lower, upper bounds */
	mpf_set_d (ra, root_bound_lo);
	mpf_set_d (rb, root_bound_hi);
	
	mpf_t fa, fb, fc, fd;
	mpf_init (fa);
	mpf_init (fb);
	mpf_init (fc);
	mpf_init (fd);

	/* Evaluate the function at the bounds */
	f (fa, ra, params);	
	f (fb, rb, params);	

	/* check that the interval bounds a zero */
	int sig_b = mpf_sign (fb);
	int sig_a = mpf_sign (fa);
	if (sig_b * sig_a > 0)
	{
		fprintf (stderr, "Error duude: endpoints don't bracket a zero\n");
		return;
	}

	/* utility values */
	mpf_t tmp, fba, fcb, fca;
	mpf_init (tmp);
	mpf_init (fba);
	mpf_init (fcb);
	mpf_init (fca);

	/* binary preecision */
	int bin = (int) (3.32192809488 * prec) +1;
	
	int which = 1;
	while (1)
	{
		switch (which)
		{
		case 1: 
			/* linear interpolation */
			mpf_sub(tmp, rb, ra);
			mpf_sub(fba, fb, fa);
			mpf_div(tmp, delt, fba);
			mpf_mul(tmp, fb);
			mpf_sub (rc, rb, tmp);
		
			f (fc, rc, params);	
			int sig_c = mpf_sign (fc);
		
			/* arrange so that c is the worst estimate */
			if (sig_c * sig_a > 0)
			{
				mp_set (tmp, ra);
				mp_set (ra, rc);
				mp_set (rc, tmp);
				mp_set (tmp, fa);
				mp_set (fa, fc);
				mp_set (fc, tmp);
			}
			else
			{
				mp_set (tmp, rb);
				mp_set (rb, rc);
				mp_set (rc, tmp);
				mp_set (tmp, fb);
				mp_set (fb, fc);
				mp_set (fc, tmp);
			}
			// which = 2;
			break;
	
		case 2: 
			/* quadratic interpolation */
			mpf_sub(fba, fb, fa);
			mpf_sub(fcb, fc, fb;
			mpf_sub(fca, fc, fa);
	
			mpf_div (tmp, a, fba);
			mpf_div (tmp, tmp, fca);
			mpf_mul (tmp, tmp, fb);
			mpf_mul (tmp, tmp, fc);
			mpf_set (rd, tmp);
	
			mpf_div (tmp, b, fba);
			mpf_div (tmp, tmp, fcb);
			mpf_mul (tmp, tmp, fa);
			mpf_mul (tmp, tmp, fc);
			mpf_sub (rd, rd, tmp);
	
			mpf_div (tmp, c, fca);
			mpf_div (tmp, tmp, fcb);
			mpf_mul (tmp, tmp, fa);
			mpf_mul (tmp, tmp, fb);
			mpf_add (rd, rd, tmp);
	
			break;
		}

double gb = mpf_get_d(rb);
double gfb = mpf_get_d(fb);
printf ("duude god f(%g)= %g\n", gb, gfb);

		/* estimate bound */
		mpf_sub (tmp, rb, ra);
		mpf_abs (tmp, tmp);
		mpf_mul_2exp (tmp, tmp, bin)
		if (mpf_cmp_ui (tmp, 1) <= 0)
		{
			mpf_set (root, rb);
			break;
		}
	}

	mpf_clear (fba);
	mpf_clear (fcb);
	mpf_clear (fca);
	mpf_clear (tmp);
	
	mpf_clear (ra);
	mpf_clear (rb);
	mpf_clear (rc);
	mpf_clear (rd);
	mpf_clear (re);

	mpf_clear (fa);
	mpf_clear (fb);
	mpf_clear (fc);
	mpf_clear (fd);
}
#endif /* BIGNUM_ROOT_FINDER */

/* ==================================================================== */

int main (int argc, char * argv[])
{
	if (argc < 3)
	{
		fprintf (stderr, "Usage: %s [ndigits] [norder]\n", argv[0]); 
		exit (1);
	}

	/* the decimal precison (number of decimal places) */
	prec = atoi (argv[1]);

	/* Number of binomial terms to sum up to */
	norder = atoi (argv[2]);

	/* compute number of binary bits this corresponds to. */
	double v = ((double) prec) *log(10.0) / log(2.0);

	/* The largest that a binomial (n,k) will get is 2^n
	 * so need an extra norder bits if going to order norder. 
	 * And pad a bit, just to be safe... */
	int bits = (int) (v + 100 + norder);
	
	/* Set the precision (number of binary bits) */
	mpf_set_default_prec (bits);
	
	printf ("# looking for precise zero locations\n");
	printf ("# computed to precision of %d decimal places\n", prec);
	printf ("# computed up to order of %d \n", norder);
	printf ("# computed with %d bits of default mpf \n", bits);
	fflush (stdout);
	
	double z[150];

	z[1] = 2.7245424726404;
	z[2] = 6.9691123091236;
	z[3] = 12.728018112883;
	z[4] = 20.018162314288;
	z[5] = 28.864209826344;
	z[6] = 39.277617847981;
	z[7] = 51.260486636309;
	z[8] = 64.812914781276;
	z[9] = 79.935164294348;
	z[10] = 96.627537143746;
	z[11] = 114.89022965254;
	z[12] = 134.72336032412;
	z[13] = 156.12700944974;
	z[14] = 179.10123558051;
	z[15] = 203.64608219242;
	z[16] = 229.76158198699;
	z[17] = 257.44776001299;
	z[18] = 286.70463579943;
	z[19] = 317.53222480004;
	z[20] = 349.93053940158;
	z[21] = 383.89958964768;
	z[22] = 419.43938376745;
	z[23] = 456.54992856739;
	z[24] = 495.23122972587;
	z[25] = 535.48329201745;
	z[26] = 577.30611948572;
	z[27] = 620.69971557796;
	z[28] = 665.66408325113;
	z[29] = 712.19922505597;
	z[30] = 760.30514320447;
	z[31] = 809.98183962424;
	z[32] = 861.22931600277;
	z[33] = 914.0475738237;
	z[34] = 968.43661439668;
	z[35] = 1024.3964388821;
	
	// z[35] = 1024.396526;
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
	z[78] = 4916.637246;
	z[79] = 5041.711939;
	z[80] = 5168.357449;
	z[81] = 5296.573727;
	z[82] = 5426.360819;
	z[83] = 5557.718682;
	z[84] = 5690.647359;
	z[85] = 5825.146838;

							 
																									
	int i;
	for (i=23; i<86; i++)
	{
		// double zer = bindary_find_zero (z[i]-0.02, z[i]+0.02, eff, 1.0e-12);
		double zer = brent_solver (z[i]-4.0, z[i]+4.0, 1.0e-16);

		printf ("%d\t%22.18g\n", i, zer);
		fflush (stdout);
	}

	return 0;
}

