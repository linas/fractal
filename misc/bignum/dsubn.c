/*
 * dsubn.c
 *
 * High-precison dsub_n  using the 
 * Gnu Multiple-precision library.
 *
 * Here, d_sub_n is sum 1/zeta * binomial
 *  
 * Linas Vepstas July 2005, July 2006
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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
	// for (p=3; p<=en; p++)
	{
		i_binomial (ibin, en, p);
		mpf_set_z (bin, ibin);
		fp_zeta (term, p, prec);
		// mpf_set_ui (term, 1);
		mpf_div (term, bin, term);

// #define LIOUVILLE
#ifdef LIOUVILLE
		fp_zeta (bin, 2*p, prec);
		mpf_mul (term, term, bin);
#endif
// #define TOTIENT
#ifdef TOTIENT
		fp_zeta (bin, p-1, prec);
		mpf_mul (term, term, bin);
#endif
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

	mpf_t term, d_n;
	mpf_init (term);
	mpf_init (d_n);

	int n;
	printf ("#\n# zeta expansion terms n^2 * (2-d_n)   \n#\n");
	printf ("# computed to precision of %d decimal places\n", prec);
	printf ("# computed up to order of %d \n", norder);
	printf ("# computed with %d bits of default mpf \n", bits);
	fflush (stdout);
	time_t then = time(0);
	for (n=2; n<=norder; n+=10)
	{
		d_sub_n (d_n, n, prec);
		mpf_set (term, d_n);

#define D_N_SCALE
#ifdef D_N_SCALE
		mpf_ui_sub (d_n, 2, d_n);
		mpf_mul_ui (term, d_n, n*n);
#endif
		time_t now = time(0);
		int elapsed = now-then;
		then = now;
		
		printf ("%d\t",n);
		fp_prt ("", term);
		printf ("\n");
		fprintf (stderr, "that took %d secs\n",elapsed);
		fflush (stdout);
	}
#endif
	

	return 0;
}

