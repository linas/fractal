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
void a_sub_n (mpf_t a_n)
{
}

/* ============================================================================= */
main ()
{
	char str[4000];
	mpz_t fact;
	mpz_init (fact);

	i_factorial (fact, 5);
	mpz_get_str (str, 10, fact);
	printf ("fact = %s\n", str);
	
	/* set the precision */
	mpf_set_default_prec (400);

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
}

