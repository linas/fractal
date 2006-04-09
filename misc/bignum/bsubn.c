/*
 * mp_zeta.c
 *
 * High-precison Riemann zeta function, using the 
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


/* ============================================================================= */

main (int argc, char * argv[])
{
	char str[4000];

#ifdef FACT_TEST
	mpz_t fact;
	mpz_init (fact);

	i_factorial (fact, 5);
	mpz_get_str (str, 10, fact);
	printf ("fact = %s\n", str);
#endif

#ifdef I_BINOMIAL_TEST
	int n, k;
	mpz_t bin;
	mpz_init (bin);

	for (n=1; n<7; n++)
	{
		for (k=0; k<=n; k++)
		{
			i_binomial (bin, n ,k);
			mpz_get_str (str, 10, bin);
			printf ("bin (%d %d) = %s\n", n, k, str);
		}
		printf ("---\n");
	}
#endif

// #define F_BINOMIAL_TEST
#ifdef F_BINOMIAL_TEST
	int n, k;
	mpf_t bin;
	mpf_init (bin);

	for (n=1; n<7; n++)
	{
		for (k=0; k<=n; k++)
		{
			fp_binomial (bin, (double)n ,k);
			printf ("bin (%d %d) = ", n, k);
			mpf_out_str (stdout, 10, 60, bin);
			printf ("\n");
		}
		printf ("---\n");
	}
#endif
	
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
	
#ifdef ZETA_STUFF
	mpf_t zeta;
	mpf_init (zeta);
	
	printf ("           000000000011111111112222222222333333333344444444445555555555666666666677777777778\n");
	printf ("           012345678901234567890123456789012345678901234567890123456789012345678901234567890\n");
	fp_zeta (zeta, 2, 13);
	fp_prt ("13 digs= ", zeta);
	fp_zeta (zeta, 8, 30);
	fp_prt ("30 digs= ", zeta);
	fp_zeta (zeta, 8, 40);
	fp_prt ("40 digs= ", zeta);
	fp_zeta (zeta, 8, 50);
	fp_prt ("50 digs= ", zeta);
	fp_zeta (zeta, 8, 60);
	fp_prt ("60 digs= ", zeta);
	fp_zeta (zeta, 8, 70);
	fp_prt ("70 digs= ", zeta);
	fp_zeta (zeta, 8, 80);
	fp_prt ("0 digs= ", zeta);
#endif
	
// #define TEST_ZETA
#ifdef TEST_ZETA
	mpf_t zeta;
	mpf_init (zeta);
	// fp_zeta_odd (zeta, 3, 180, 7, 360, 0, 60); 
	// fp_prt ("duude zeta3= ", zeta);
	int i;
	int pr = prec;
	for (i=3; i<nterms; i++ ) {
		fp_zeta (zeta, i, pr);
		printf ("char * zeta_%d_%d = \"", i, pr);
		mpf_out_str (stdout, 10, pr, zeta);
		printf ("\";\n");
		fflush (stdout);
	}
#endif

#ifdef TEST_EXP
	mpf_t ex, one;
	mpf_init (ex);
	mpf_init (one);
	mpf_set_ui(one, 1);
	fp_exp (ex, one, 50);
	fp_prt ("e= ", ex);
	mpf_clear (ex);
	mpf_clear(one);
#endif
	
#ifdef TEST_DIGIT_COUNT
	mpz_t ival, tmpa, tmpb;
	mpz_init (ival);
	mpz_init (tmpa);
	mpz_init (tmpb);
	mpz_set_ui (ival, 3000000);
	int nd = num_digits (ival, tmpa, tmpb);
	printf ("found %d digits\n", nd);
#endif

// #define TEST_BERNOULLI
#ifdef TEST_BERNOULLI
	mpq_t bern;
	mpq_init (bern);
	int n = 4;
	for (n=8; n<30; n++) 
	{
		q_bernoulli (bern, n);
		printf ("bernoulli (%d)= ", n);
		mpq_out_str (stdout, 10, bern);
		printf ("\n");
	}
#endif /* TEST_BERNOULLI */
	
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

