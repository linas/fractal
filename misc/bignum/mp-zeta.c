/*
 * mp_zeta.c
 *
 * High-precison Riemann zeta function, using the 
 * Gnu Multiple-precision library.
 * 
 * Linas Vepstas July 2005
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>

void fp_prt (char * str, mpf_t val)
{
	printf (str);
	mpf_out_str (stdout, 10, 0, val);
	printf ("\n");
}

/* ============================================================================= */
/* i_poch_rising
 * rising pochhammer symbol, for integer values.
 *
 * Brute force, simple.
 */

void i_poch_rising (mpz_t poch, unsigned int k, unsigned int n)
{
	mpz_t acc;
	
	mpz_init (acc);

	mpz_set_ui (poch, 1);
	unsigned int i;
	for (i=0; i<n; i++)
	{
		mpz_mul_ui (acc, poch, i+k);
		mpz_set (poch, acc);
	}

	mpz_clear (acc);
}

/* i_factorial -- the factorial
 */
void i_factorial (mpz_t fact, unsigned int n)
{
	i_poch_rising (fact, 1, n);
}

/* ============================================================================= */
/* fp_binomial
 * Binomial coefficient
 */

void i_binomial (mpz_t bin, unsigned int n, unsigned int k)
{
	mpz_t top, bot;

	if (2*k < n) k = n-k;

	mpz_init (top);
	mpz_init (bot);
	i_poch_rising (top, k+1, n-k);
	i_factorial (bot, n-k); 

	mpz_divexact (bin, top, bot);
	
	mpz_clear (top);
	mpz_clear (bot);
}

/* ============================================================================= */
/* fp_euler
 * return Euler-Mascheroni const
 */
void fp_euler_mascheroni (mpf_t gam)
{
	char * g="0.5772156649015328606065120900824024310421593359399235988057672348848677267776646709369470632917467495e0";
	
	mpf_set_str (gam, g, 10);
}

/* ============================================================================= */
/* fp_zeta
 * Floating-point-valued Riemann zeta for positive integer arguments 
 * return value placed in the arg "zeta".
 *
 * Simple-minded algo, carries out math to prec decimal digits
 */
void fp_zeta (mpf_t zeta, unsigned int s, int prec)
{
	unsigned long int us = s;
	mpf_t acc;
	mpf_t term;
	mpf_t en;
	mpf_t inv;
	
	mpf_init (acc);
	mpf_init (term);
	mpf_init (en);
	mpf_init (inv);
	
	mpf_set_ui (zeta, 1);

	/* Compute number of terms to be carried out 
	 * However, this estimate is wrong; it stops 
	 * when next term is "smaller than" rather 
	 * than when its converged.
	 */
	double fprec = prec;
	fprec /= (double) s;
	double dig = pow (10.0, fprec);
	if (1.0e10 < dig)
	{
		fprintf (stderr, "Sorry bucko, can't do it\n");
		return;
	}
	int nmax = dig+1.0;
	printf ("zeta will be computed with %d terms\n", nmax);
	
	int n;
	for (n=2; n<= nmax; n++)
	{
		mpf_set_ui (en, n);
		mpf_ui_div (inv, 1, en);  /* inv = 1/n */
		mpf_pow_ui (term, inv, us); /* term = 1/n^s */
		mpf_add (acc, zeta, term);
		mpf_set (zeta, acc);
	}

	mpf_clear (acc);
	mpf_clear (term);
	mpf_clear (en);
	mpf_clear (inv);
}

/* ============================================================================= */
/* compute a_sub_n
 */
void a_sub_n (mpf_t a_n, unsigned int n, unsigned int prec)
{
	int k;
	mpf_t fbin, term, zt, ok, one, acc, zeta;
	mpf_t gam;

	mpf_init (term);
	mpf_init (acc);
	mpf_init (zeta);
	mpf_init (zt);
	mpf_init (ok);
	mpf_init (one);
	mpf_init (fbin);
	mpf_init (gam);
	
	mpf_set_ui (one, 1);

	mpz_t ibin;
	mpz_init (ibin);

	for (k=1; k<= n; k++)
	{
		fp_zeta (zeta, k+1, prec);
		mpf_div_ui (zt, zeta, k+1);
		mpf_div_ui (ok, one, k);
		mpf_add (term, zt, ok);

		i_binomial (ibin, n, k);
		mpf_set_z (fbin, ibin);

		mpf_mul (zeta, term, fbin);

		if (k%2) mpf_neg (term, zeta);
		
		mpf_add (acc, a_n, term);
		mpf_set (a_n, acc);
	}

	/* add const terms */
	mpf_add_ui (term, a_n, 1);
	fp_euler_mascheroni (gam);
	mpf_sub (a_n, term, gam);

	/* subtract 1/2(n+1) */
	mpf_div_ui (ok, one, 2*(n+1));
	mpf_sub (term, a_n, ok);
	mpf_set (a_n, term);
	
	mpf_clear (term);
	mpf_clear (acc);
	mpf_clear (zeta);
	mpf_clear (zt);
	mpf_clear (ok);
	mpf_clear (one);
	mpf_clear (fbin);
	mpf_clear (gam);

	mpz_clear (ibin);
}

/* ============================================================================= */
main ()
{
	char str[4000];

#ifdef FACT_TEST
	mpz_t fact;
	mpz_init (fact);

	i_factorial (fact, 5);
	mpz_get_str (str, 10, fact);
	printf ("fact = %s\n", str);
#endif

#ifdef BINOMIAL_TEST
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
	
	/* set the precision */
	mpf_set_default_prec (400);

	mpf_t a_n;
	mpf_init (a_n);

	int prec = 10;
	a_sub_n (a_n, 0, prec);
	fp_prt ("a_0= ", a_n);
	a_sub_n (a_n, 1, prec);
	fp_prt ("a_1= ", a_n);
	a_sub_n (a_n, 2, prec);
	fp_prt ("a_2= ", a_n);
	
#ifdef ZETA_STUFF
	mpf_t zeta;
	mpf_init (zeta);
	
	printf ("           000000000011111111112222222222333333333344444444445555555555666666666677777777778\n");
	printf ("           012345678901234567890123456789012345678901234567890123456789012345678901234567890\n");
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
}

