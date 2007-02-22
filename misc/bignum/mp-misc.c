/*
 * mp-misc.c
 *
 * High-precison misc functions, using the 
 * Gnu Multiple-precision library.
 *
 * Copyright (C) 2005 Linas Vepstas
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301  USA
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

void cpx_prt (char * str, const cpx_t val)
{
	printf (str);
	mpf_out_str (stdout, 10, 30, val[0].re);
	printf (" + i ");
	mpf_out_str (stdout, 10, 30, val[0].im);
}

void ecpx_prt (char * str, const cpx_t val)
{
	fprintf (stderr, str);
	mpf_out_str (stderr, 10, 30, val[0].re);
	fprintf (stderr, " + i ");
	mpf_out_str (stderr, 10, 30, val[0].im);
}

/* ===================================================== */
/**
 * fp_epsilon - return 10^{-prec} 
 */
void fp_epsilon (mpf_t eps, int prec)
{
	static int cache_prec = -1;
	static mpf_t cache_eps;

	if (-1 == cache_prec)
	{
		mpf_init (cache_eps);
	}

	if (prec == cache_prec)
	{
		mpf_set (eps, cache_eps);
		return;
	}
	mpf_set_prec (cache_eps, 3.322*prec+50);

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

/* ===================================================== */

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

