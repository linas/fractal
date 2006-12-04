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

