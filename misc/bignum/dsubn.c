/*
 * dsubn.c
 *
 * High-precison dsub_n  using the 
 * Gnu Multiple-precision library.
 *
 * Here, d_sub_n is d_n = sum_{k=2}^n (-1)^k {n \choose k} /zeta(k)
 * which is given by the Mobious function.
 *
 * Similar sums for the totient ad Liouville function are also
 * implemented.
 *  
 * Linas Vepstas July 2005, July 2006
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mp-binomial.h"
#include "mp-consts.h"
#include "mp-misc.h"
#include "mp_zeta.h"

#ifdef TESTING_123
void bino (void)
{
	mpz_t ibin;
	mpz_init (ibin);

	int n;
	for (n=0;n<2000; n++)
	{
		int k;
		for (k=0; k<=n; k++)
		{
			// printf ("%d %d  ", n, k);
			i_binomial_sequence (ibin, n, k);
			// i_binomial (ibin, n, k);
			// i_prt ("", ibin);
			// printf ("\n");
		}
	}
	mpz_clear (ibin);
}
#endif

/* ==================================================================== */
/* Provide summation of 
 * sum_{p=2}^n (-1)^p {n \choose p} \frac{1}{p}
 */

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
		i_binomial_sequence (ibin, en, p);
		mpf_set_z (bin, ibin);
		fp_zeta (term, p, prec);
		// mpf_set_ui (term, 1);
		mpf_div (term, bin, term);

		fp_zeta (bin, 2*p, prec);
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

#define D_N_SCALE
#ifdef D_N_SCALE
	mpf_ui_sub (acc, 2, acc);
	mpf_mul_ui (acc, acc, en*en);
#endif
	mpf_clear (bin);
	mpf_clear (term);
	mpz_clear (ibin);
}

/* ==================================================================== */
/* Same as above, but for totient function, instead of mobius
 * That is, its for zeta(n-1)/zeta(n)
 */

void d_totient_n (mpf_t acc, int en, unsigned int prec)
{
	mpz_t ibin;
	mpz_init (ibin);

	mpf_t bin, term;
	mpf_init (bin);
	mpf_init (term);

	mpf_set_ui (acc, 0);
	
	int p;
	for (p=3; p<=en; p++)
	{
		i_binomial_sequence (ibin, en, p);
		mpf_set_z (bin, ibin);
		fp_zeta (term, p, prec);
		// mpf_set_ui (term, 1);
		mpf_div (term, bin, term);

		fp_zeta (bin, p-1, prec);
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
#define T_N_SCALE
#ifdef T_N_SCALE
	mpf_div_ui (acc, acc, en*en);
	// mpf_set_ui (bin, en);
	// fp_log (term, bin, prec);
	// mpf_div (acc, acc, term);
#endif

	mpf_clear (bin);
	mpf_clear (term);
	mpz_clear (ibin);
}

/* ==================================================================== */
/* Same as above, but for the Liouville function, instead of
 * the totiet function. That is, its for zeta(2n)/zeta(n)
 */

void d_liouville_n (mpf_t acc, int en, unsigned int prec)
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
		i_binomial_sequence (ibin, en, p);
		mpf_set_z (bin, ibin);
		fp_zeta (term, p, prec);
		// mpf_set_ui (term, 1);
		mpf_div (term, bin, term);

		fp_zeta (bin, 2*p, prec);
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

#define L_N_SCALE
#ifdef L_N_SCALE
	mpf_add_ui (acc, acc, 1);
	
	mpz_fac_ui (ibin, en);
	mpf_set_z (term, ibin);
	
	fp_poch_rising (bin, 0.5, en);
	mpf_div (term, term, bin);
	
	fp_zeta_half (bin, 60);
	mpf_div (term, term, bin);

	mpf_add (acc, acc, term);
#endif

	mpf_clear (bin);
	mpf_clear (term);
	mpz_clear (ibin);
}

/* ==================================================================== */

int main (int argc, char * argv[])
{
	if (argc < 4)
	{
		fprintf (stderr, "Usage: %s <ndigits> <nlow> <nhigh>\n", argv[0]); 
		exit (1);
	}

	/* the decimal precison (number of decimal places) */
	int prec = atoi (argv[1]);

	/* number of d_n's to compute */
	int nlow = atoi (argv[2]);
	int nhigh = atoi (argv[3]);
	int norder = nhigh;

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
	// printf ("#\n# zeta expansion terms n^2 * (2-d_n)   \n#\n");
	// printf ("#\n# totient zeta expansion terms d_n/(n*n)   \n#\n");
	// printf ("#\n# liouville zeta expansion terms d_n   \n#\n");
	printf ("#\n# liouville zeta expansion terms d_n + 1 + gamma(n+0.5)/gamma(n+1)/zeta(1/2)\n#\n");
	printf ("# computed to precision of %d decimal places\n", prec);
	printf ("# computed up to order of %d \n", norder);
	printf ("# computed with %d bits of default mpf \n", bits);
	fflush (stdout);
	time_t then = time(0);
	
	mpz_t ibin;
	mpz_init (ibin);
	i_binomial_sequence (ibin, 0, 0);
	mpz_clear (ibin);

	for (n=nlow; n<=nhigh; )
	{
		// d_sub_n (d_n, n, prec);
		// d_totient_n (d_n, n, prec);
		d_liouville_n (d_n, n, prec);

		time_t now = time(0);
		int elapsed = now-then;
		then = now;
		
		printf ("%d\t",n);
		fp_prt ("", d_n);
		printf ("\t%d\n", elapsed);
		fprintf (stderr, "n=%d took %d secs\n", n, elapsed);
		fflush (stdout);

		// int step = 1+n/10;
		// int step = 1+n/50;
		int step = 1+n/100;
		n += step;
	}
#endif

	return 0;
}

