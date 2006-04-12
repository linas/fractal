/*
 * stieltjes.c
 *
 * High-precison Steiltjes constants
 * Gnu Multiple-precision library.
 *
 * Linas Vepstas April 2006
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mp_zeta.h"

/* ==================================================================== */
/* Return the Steiltjes constants */

void stieltjes_gamma (mpf_t gam, int n, int prec, int nterms)
{
	int k;

	mpz_t isb;
	mpz_init (isb);

	mpf_t term, sb;
	mpf_init (term);
	mpf_init (sb);

	mpf_set_ui (gam, 0);
	// XXXX precision violation !!
	for (k=n; k<n+nterms; k++)
	{
		b_sub_n (term, k, prec);
		i_stirbin_sum (isb, k,n);
		mpf_set_z (sb, isb);
		mpf_mul (term, term, sb);

		i_factorial (isb, k);
		mpf_set_z (sb, isb);
		mpf_div (term, term, sb);
		mpf_add (gam, gam, term);
	}
	i_factorial (isb, n);
	mpf_set_z (sb, isb);
	mpf_mul (gam, gam, sb);
	if (n%2) mpf_neg (gam, gam);

	mpf_clear (term);
	mpf_clear (sb);
	mpz_clear (isb);
}

/* ==================================================================== */

main (int argc, char * argv[])
{
	char str[4000];

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
	int bits = (int) (v + 100 + norder);
	
	/* set the precision (number of binary bits) */
	mpf_set_default_prec (bits);

	printf ("# Computing Stieltjes constants\n");
	printf ("# computed to precision of %d decimal places\n", prec);
	printf ("# computed up to order of %d \n", norder);
	printf ("# computed with %d bits of default mpf \n", bits);

	mpf_t stie;
	mpf_init (stie);
	int i;
	for (i=0; i<80; i++ ) {
		stieltjes_gamma (stie, i, prec, norder);
		printf ("gamma[%d] = ", i);
		mpf_out_str (stdout, 10, 60, stie);
		printf (";\n");
		fflush (stdout);
	}

}

/* =============================== END OF FILE =========================== */

