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

double d_sub_n (double z, unsigned int prec, unsigned int norder)
{
	mpf_t acc, bin, term;
	mpf_init (acc);
	mpf_init (bin);
	mpf_init (term);

	mpf_set_ui (acc, 0);
	int p;
	for (p=2; p<norder; p++)
	{
		fp_binomial (bin, z, p);
		fp_zeta (term, p+2, prec);
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

// printf ("duude sum=%g\n", sum);

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

// #define B_N_SCALE
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
		fflush (stdout);
	}
#endif
	

	return 0;
}

