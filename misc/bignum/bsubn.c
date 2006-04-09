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

#include "mp_zeta.h"

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

	/* number of an's to compute */
	int nterms = atoi (argv[2]);

	/* compute number of binary bits this corresponds to. */
	double v = ((double) prec) *log(10.0) / log(2.0);

	/* the variable-precision calculations are touchy about this */
	/* XXX this should be stirling's approx for binomial */ 
	int bits = (int) (v + 300 + 3*nterms);
	
	/* set the precision (number of binary bits) */
	mpf_set_default_prec (bits);
	
	
#define A_SUB_N
#ifdef A_SUB_N

#ifdef PRECOMPUTE
	/* precompute values */
	mpf_t zeta;
	mpf_init (zeta);
	int i;
	int pr = prec;
	for (i=3; i<nterms; i++ ) {
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
	fp_pi (pi);
	
	int n;
	printf ("#\n# zeta expansion terms \n#\n");
	printf ("# computed with variable precision of %d decimal places\n", prec);
	printf ("# computed with %d bits of default mpf \n", bits);
	for (n=0; n<nterms; n++)
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
#endif
		
		printf ("%d\t",n+1);
		fp_prt ("", prod);
		// fp_prt ("", a_n);
		// fp_prt ("", b_n);
		fflush (stdout);
	}
#endif
	
// #define A_SUB_S
#ifdef A_SUB_S

	mpf_t re_a, im_a;
	mpf_init (re_a);
	mpf_init (im_a);

	int n;
	printf ("#\n# zeta expansion terms \n#\n");
	printf ("# computed with variable precision of %d decimal places\n", prec);
	printf ("# computed with %d bits of default mpf \n", bits);
	for (n=0; n<nterms; n++)
	{
		double re_s = 0;
		double im_s = -10+n/10.0;
		a_sub_s (re_a, im_a, re_s, im_s, prec);

		printf ("%d\t%12.9g\t%12.9g\t", n, re_s, im_s);
		mpf_out_str (stdout, 10, 21, re_a);
		printf ("\t");
		mpf_out_str (stdout, 10, 21, im_a);
		printf ("\n");
		fflush (stdout);
	}
#endif

}

