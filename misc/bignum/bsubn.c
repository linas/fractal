/*
 * bsubn.c
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
#include <gsl/gsl_sf_psi.h>

#include "mp-binomial.h"
#include "mp-consts.h"
#include "mp-misc.h"
#include "mp-trig.h"
#include "mp-zeta.h"

double harmonic (double z, unsigned int prec, unsigned int norder)
{
	mpf_t acc, bin, term;
	mpf_init (acc);
	mpf_init (bin);
	mpf_init (term);

	mpf_set_ui (acc, 0);
	int p;
	for (p=2; p<norder; p++)
	{
		fp_binomial_d (bin, z, p);
		mpf_div_ui (term, bin, p-1);
		if (p%2)
		{
			mpf_sub (acc, acc, term);
		}
		else
		{
			mpf_add (acc, acc, term);
		} 
	}
	double sum = mpf_get_d (acc);

	sum -= 1.0;

	double harm = 0.57721566490153286060-1.0;
	harm += gsl_sf_psi (z);
	sum -= z*harm;
	
// printf ("duude sum=%g\n", sum);

	mpf_clear (term);
	mpf_clear (bin);
	mpf_clear (acc);

	return sum;
}

double b_functional (double z, unsigned int prec, unsigned int norder)
{
	mpf_t acc, bin, term;
	mpf_init (acc);
	mpf_init (bin);
	mpf_init (term);

	z = 1.0 - z;
	mpf_set_ui (acc, 0);
	int p;
	for (p=0; p<norder; p++)
	{
		fp_zeta (term, p+2, prec);
		mpf_sub_ui (term, term, 1);
		fp_binomial_d (bin, z-1.0, p+1);
		mpf_mul (term, term, bin);
		mpf_add (acc, acc, term);
	}
	double sum = mpf_get_d (acc);

	sum -= 1.5;
	sum += pow (2.0, z-1.0);

	double harm = 2.0*(0.57721566490153286060-1.0);
	harm += gsl_sf_psi (z);
	harm += M_PI / tan (M_PI*z);
	sum += (z-1.0)*harm;
	
// printf ("duude sum=%g\n", sum);

	mpf_clear (term);
	mpf_clear (bin);
	mpf_clear (acc);

	return sum;
}

double b_fa (double z, unsigned int prec, unsigned int norder)
{
	mpf_t acc, bin, term;
	mpf_init (acc);
	mpf_init (bin);
	mpf_init (term);

	mpf_set_ui (acc, 0);
	int p;
	for (p=2; p<norder; p++)
	{
		fp_zeta (term, p, prec);
		mpf_set_ui (bin, 1);
		mpf_div_ui (bin, bin, p-1);
		mpf_sub (term, term, bin);
		fp_binomial_d (bin, z, p);
		mpf_mul (term, term, bin);
		if (p%2)
		{
			mpf_sub (acc, acc, term);
		}
		else
		{
			mpf_add (acc, acc, term);
		} 
	}
	double sum = mpf_get_d (acc);

	sum += 0.5;
	sum -= z*0.57721566490153286060;

	mpf_clear (term);
	mpf_clear (bin);
	mpf_clear (acc);

	return sum;
}

double b_fb (double z, unsigned int prec, unsigned int norder)
{
	mpf_t acc, bin, term;
	mpf_init (acc);
	mpf_init (bin);
	mpf_init (term);

	mpf_set_ui (acc, 0);
	int p;
	for (p=2; p<norder; p++)
	{
		fp_zeta (term, p, prec);
		fp_binomial_d (bin, z, p);
		mpf_mul (term, term, bin);
		if (p%2)
		{
			mpf_sub (acc, acc, term);
		}
		else
		{
			mpf_add (acc, acc, term);
		} 
	}
	double sum = mpf_get_d (acc);

	sum += 0.5;
	sum -= z*0.57721566490153286060;

	double harm = 0.57721566490153286060-1.0;
	harm += gsl_sf_psi (z);
	sum -= 1.0+ z*harm;
	
	mpf_clear (term);
	mpf_clear (bin);
	mpf_clear (acc);

	return sum;
}

double b_fc (double z, unsigned int prec, unsigned int norder)
{
	mpf_t acc, bin, term;
	mpf_init (acc);
	mpf_init (bin);
	mpf_init (term);

	mpf_set_ui (acc, 0);
	int p;
	for (p=2; p<norder; p++)
	{
		fp_zeta (term, p, prec);
		mpf_sub_ui (term, term, 1);
		fp_binomial_d (bin, z, p);
		mpf_mul (term, term, bin);
		if (p%2)
		{
			mpf_sub (acc, acc, term);
		}
		else
		{
			mpf_add (acc, acc, term);
		} 
	}
	double sum = mpf_get_d (acc);

#if 0
	sum += 0.5;
	sum -= z*0.57721566490153286060;

	double harm = 0.57721566490153286060-1.0;
	harm += gsl_sf_psi (z);
	sum -= 1.0+ z*harm;

	sum += z-1.0;
#endif

	sum -= 1.5;

	double harm = 2.0* (1.0-0.57721566490153286060);
	harm -= gsl_sf_psi (z);
	sum += z*harm;

	mpf_clear (term);
	mpf_clear (bin);
	mpf_clear (acc);

	return sum;
}

double b_b (double z, unsigned int prec, unsigned int norder)
{
	mpf_t acc, bin, term;
	mpf_init (acc);
	mpf_init (bin);
	mpf_init (term);

	mpf_set_ui (acc, 0);
	int p;
	for (p=0; p<norder; p++)
	{
		fp_binomial_d (term, z, p);
		if (p%2)
		{
			mpf_sub (acc, acc, term);
		}
		else
		{
			mpf_add (acc, acc, term);
		} 
	}
	double sum = mpf_get_d (acc);

	mpf_clear (term);
	mpf_clear (bin);
	mpf_clear (acc);

	return sum;
}

/* ==================================================================== */

int main (int argc, char * argv[])
{
	if (argc < 3)
	{
		fprintf (stderr, "Usage: %s [ndigits] [norder]\n", argv[0]); 
		exit (1);
	}

	/* the decimal precison (number of decimal places) */
	int prec = atoi (argv[1]);

	/* number of an's to compute */
	int norder = atoi (argv[2]);

	/* compute number of binary bits this corresponds to. */
	double v = ((double) prec) *log(10.0) / log(2.0);

	/* The largest that a binomial (n,k) will get is 2^n
	 * so need an extra norder bits if going to order norder. 
	 * And pad a bit, just to be safe... */
	// int bits = (int) (v + 100 + 3.3*norder);
	int bits = (int) (v + 100 + norder);
	
	/* set the precision (number of binary bits) */
	mpf_set_default_prec (bits);
	
	
// #define A_SUB_N
#ifdef A_SUB_N

#ifdef PRECOMPUTE
	/* precompute values */
	mpf_t zeta;
	mpf_init (zeta);
	int i;
	int pr = prec;
	for (i=3; i<norder; i++ ) {
		fp_zeta (zeta, i, pr);
	}
#endif /* PRECOMPUTE */
	
	mpf_t a_n, b_n, en, pi, sq, term, p_n, prod, w;
	mpf_init (a_n);
	mpf_init (b_n);
	mpf_init (pi);
	mpf_init (en);
	mpf_init (sq);
	mpf_init (term);
	mpf_init (p_n);
	mpf_init (prod);
	mpf_init (w);

	// The standard w value is 1 ... 
	mpf_set_d (w, 1.0);
	fp_pi (pi, prec);
	
	int n;
	printf ("#\n# zeta expansion -- using b_n = n*a_{n-1} variant \n#\n");
	printf ("# computed to precision of %d decimal places\n", prec);
	printf ("# computed up to order of %d \n", norder);
	printf ("# computed with %d bits of default mpf \n", bits);
	fflush (stdout);
	// for (n=2590; n<2600; n++)
	for (n=0; n<norder; n++)
	{
		a_sub_n (a_n, w, n, prec);

// #define EXACT_BND
#ifdef EXACT_BND
		/* compute the bound */
		mpf_set_ui (en, n+1);
		mpf_sqrt (sq, en);
		mpf_mul_ui (term, sq, 4);
		mpf_neg (en, term);
		fp_exp (p_n, en, prec);
		mpf_div (prod, a_n, p_n);
#endif
		
#ifdef FLT_BND
		double dbn = 1.0/exp (-4.0*sqrt (n+1));
		mpf_set_d (b_n, dbn);
		mpf_mul(prod, a_n, b_n);
#endif

		// b_n = n a_{n-1}
		mpf_mul_ui (b_n, a_n, n+1);

#define B_N_SCALE
#ifdef B_N_SCALE
		mpf_set_ui (en, n+1);
		mpf_mul_ui (term, en, 4);
		mpf_mul (term, term, pi);
		mpf_sqrt (sq, term);
		mpf_neg (sq, sq);
		fp_exp (p_n, sq, prec);
		mpf_div (prod, b_n, p_n);
		mpf_sqrt (sq, en);
		mpf_sqrt (sq, sq);
		mpf_div (prod, prod, sq);
#endif /* B_N_SCALE */
		
		printf ("%d\t",n+1);
		fp_prt ("", prod);
		// fp_prt ("", a_n);
		// fp_prt ("", b_n);
		printf("\n");
		fflush (stdout);
	}
#endif /* A_SUB_N */
	
#define B_SUB_N
#ifdef B_SUB_N

	mpf_t b_n, en, pi, sq, term, p_n, prod;
	mpf_init (b_n);
	mpf_init (pi);
	mpf_init (en);
	mpf_init (sq);
	mpf_init (term);
	mpf_init (p_n);
	mpf_init (prod);

	fp_pi (pi, prec);
	
	int n;
	printf ("#\n# zeta expansion terms b_n straight up. \n#\n");
	printf ("# computed to precision of %d decimal places\n", prec);
	printf ("# computed up to order of %d \n", norder);
	printf ("# computed with %d bits of default mpf \n", bits);
	fflush (stdout);
	for (n=1; n<norder; n++)
	{
		b_sub_n (b_n, n, prec);

#define B_N_SCALE
#ifdef B_N_SCALE
		mpf_set_ui (en, n);
		mpf_mul_ui (term, en, 4);
		mpf_mul (term, term, pi);
		mpf_sqrt (sq, term);
		mpf_neg (sq, sq);
		fp_exp (p_n, sq, prec);
		mpf_div (prod, b_n, p_n);
		mpf_sqrt (sq, en);
		mpf_sqrt (sq, sq);
		mpf_div (prod, prod, sq);
#endif
		
		printf ("%d\t",n);
		fp_prt ("", prod);
		// fp_prt ("", a_n);
		// fp_prt ("", b_n);
		printf("\n");
		fflush (stdout);
	}
#endif /* B_SUB_N */
	
// #define A_SUB_S
#ifdef A_SUB_S

	mpf_t re_a, im_a;
	mpf_init (re_a);
	mpf_init (im_a);

	int n;
	printf ("#\n# zeta expansion terms \n#\n");
	printf ("# computed with variable precision of %d decimal places\n", prec);
	printf ("# computed up to order of %d \n", norder);
	printf ("# computed with %d bits of default mpf \n", bits);
	fflush (stdout);
	for (n=0; n<20; n++)
	{
		// double re_s = 0;
		// double im_s = -10.0+n/10.0;
		double re_s = 1.0 - 0.1*n;
		double im_s = 0.0;
		// a_sub_s (re_a, im_a, re_s, im_s, prec);
		b_sub_s (re_a, im_a, re_s, im_s, prec, norder);

		printf ("%d\t%12.9g\t%12.9g\t", n, re_s, im_s);
		mpf_out_str (stdout, 10, 21, re_a);
		printf ("\t");
		mpf_out_str (stdout, 10, 21, im_a);
		printf ("\n");
		fflush (stdout);
	}
#endif

// #define B_FUNC
#ifdef B_FUNC
	mpf_t re_a, im_a;
	mpf_init (re_a);
	mpf_init (im_a);

	int n;
	for (n=1; n<50; n++)
	{
		double re_s = -1.05 - 0.1*n;
		double im_s = 0.0;
		// b_sub_s (re_a, im_a, re_s, im_s, prec, norder);
		double bs = mpf_get_d (re_a);
		double bz = b_b (re_s, prec, norder);
		// double bf = b_fa (re_s, prec, norder);
		// double bf = b_fb (re_s, prec, norder);
		// double bf = bz + b_fc (re_s, prec, norder);
		double bf = b_functional (re_s, prec, norder);
		printf ("z=%8.4g   bs=%g  \tbf=%g\n", re_s, bs, bf);
	}
#endif

// #define TEST
#ifdef TEST
	mpf_t s, bin, term, acc;
	mpf_init (s);
	mpf_init (bin);
	mpf_init (term);
	mpf_init (acc);
	printf ("# computed to precision of %d decimal places\n", prec);
	printf ("# computed with %d bits of default mpf \n", bits);

	mpf_set_ui (s, 1);
	mpf_div_ui (s, s, 2);
	mpf_set_ui (acc, 0);

	double es = -0.5;
	int n;
	for (n=0; n<5000; n++)
	{
		fp_binomial_d (bin, es, n);
		if (n%2)
		{
			mpf_sub (acc, acc, bin);
		}
		else
		{
			mpf_add (acc, acc, bin);
		}

		printf ("n= %d  ", n);
		mpf_out_str (stdout, 10, 30, acc);
		printf ("\n");
		fflush (stdout);
	}
#endif


	return 0;
}

