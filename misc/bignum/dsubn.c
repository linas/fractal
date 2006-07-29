/*
 * dsubn.c
 *
 * High-precison dsub_n  using the 
 * Gnu Multiple-precision library.
 *
 * Here, d_sub_n is sum 1/zeta * sunomial
 *  
 * Linas Vepstas July 2005, July 2006
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sf_psi.h>

#include "mp_zeta.h"

void d_sub_n (mpf_t acc, int en, unsigned int prec)
{
	mpz_t ibin;
	mpz_init (ibin);

	mpf_t bin, term;
	mpf_init (bin);
	mpf_init (term);

	mpf_set_ui (acc, 0);
	int p;
	for (p=2; p<=en; p++)
	{
		i_binomial (ibin, en, p);
		mpf_set_z (bin, ibin);
		fp_zeta (term, p, prec);
		// mpf_set_ui (term, 1);
		mpf_div (term, bin, term);
		if (p%2)
		{
			mpf_sub (acc, acc, term);
		}
		else
		{
			mpf_add (acc, acc, term);
		} 
	}

	mpf_clear (bin);
	mpf_clear (term);
	mpz_clear (ibin);
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

	/* number of d_n's to compute */
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
	
	
#define D_SUB_N
#ifdef D_SUB_N

	mpf_t term, d_n, en, pi;
	mpf_init (term);
	mpf_init (d_n);
	mpf_init (pi);
	mpf_init (en);

	fp_pi (pi, prec);
	
	int n;
	printf ("#\n# zeta expansion terms d_n straight up. \n#\n");
	printf ("# computed to precision of %d decimal places\n", prec);
	printf ("# computed up to order of %d \n", norder);
	printf ("# computed with %d bits of default mpf \n", bits);
	fflush (stdout);
	for (n=2; n<=norder; n++)
	{
		d_sub_n (d_n, n, prec);

#define D_N_SCALE
#ifdef D_N_SCALE
		mpf_ui_sub (d_n, 2, d_n);
		mpf_ui_div (term, 1, d_n);
		mpf_div_ui (term, term, n*n);
		mpf_ui_div (term, 1, term);
#endif
		
		printf ("%d\t",n);
		fp_prt ("", term);
		fflush (stdout);
	}
#endif
	

	return 0;
}

