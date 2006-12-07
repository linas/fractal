/*
 * mp-misc.c
 *
 * High-precison misc functions, using the 
 * Gnu Multiple-precision library.
 *
 * Linas Vepstas July 2005
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>
#include "mp-complex.h"
#include "mp-misc.h"

void i_prt (char * str, mpz_t val)
{
	printf (str);
	mpz_out_str (stdout, 10, val);
}

void fp_prt (char * str, mpf_t val)
{
	printf (str);
	mpf_out_str (stdout, 10, 60, val);
}

void cpx_prt (char * str, const cpx_t const val)
{
	printf (str);
	mpf_out_str (stdout, 10, 30, val[0].re);
	printf (" + i ");
	mpf_out_str (stdout, 10, 30, val[0].im);
}

/**
 * fp_epsilon - return 10^{-prec} 
 */
void fp_epsilon (mpf_t eps, int prec)
{
	static int cache_prec = -1;
	mpf_t cache_eps;

	if (-1 == cache_prec)
	{
		mpf_init (cache_eps);
	}

	if (prec == cache_prec)
	{
		mpf_set (eps, cache_eps);
		return;
	}

	/* double mex = ((double) prec) * log (10.0) / log(2.0); */
	double mex = ((double) prec) * 3.321928095;
	unsigned int imax = (unsigned int) (mex +1.0);
	mpf_t one;
	mpf_init (one);
	mpf_set_ui (one, 1);
	mpf_div_2exp (cache_eps, one, imax);

	mpf_set (eps, cache_eps);
	cache_prec = prec;
	mpf_clear (one);
}

/* prec is the decimal precison (number of decimal places) */
/* nterms is the number of an's to compute */
void set_bits (int prec, int nterms)
{
	/* Compute number of binary bits this corresponds to. */
	double v = ((double) prec) *log(10.0) / log(2.0);

	/* The variable-precision calculations are touchy about this */
	/* XXX this should be stirling's approx for binomial */
	int bits = (int) (v + 300 + 3*nterms);

	/* Set the precision (number of binary bits) */
	mpf_set_default_prec (bits);
}

/* =============================== END OF FILE =========================== */

