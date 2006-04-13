/*
 * mp_zeta_test.c
 *
 * Small test suite for 
 * High-precison Riemann zeta function, using the 
 * Gnu Multiple-precision library.
 *
 * Currently not automated, done purely by visual inspection
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

// #define I_STIRLING_TEST
#ifdef I_STIRLING_TEST
	int n, k;
	mpz_t sitrly;
	mpz_init (sitrly);

	for (n=0; n<21; n++)
	{
		for (k=0; k<=n; k++)
		{
			i_stirling_first (sitrly, n ,k);
			mpz_get_str (str, 10, sitrly);
			printf ("sitrly (%d %d) = %s\n", n, k, str);
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

// #define TEST_STIELTJES
#ifdef TEST_STIELTJES
	mpf_t stie;
	mpf_init (stie);
	int i;
	for (i=0; i<40; i++ ) {
		stieltjes_gamma (stie, i);
		printf ("gamma[%d] = ", i);
		mpf_out_str (stdout, 10, 60, stie);
		printf (";\n");
		fflush (stdout);
	}
#endif

}

/* =============================== END OF FILE =========================== */

