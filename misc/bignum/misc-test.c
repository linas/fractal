/*
 * misc-test.c
 *
 * Spot checks.
 *
 * Copyright (C) 2010 Linas Vepstas
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
#include "mp-zeta.h"


void test_complex_riemann_zeta (int nterms, int prec)
{
	/* This test was passing in Nov 2006 */
	mpf_t nzeta;
	mpf_init (nzeta);
	
	cpx_t zeta, ess;
	cpx_init (zeta);
	cpx_init (ess);

	int i = nterms;

	double z_expected, z_got;
	fp_zeta (nzeta, i, prec);
	z_expected = mpf_get_d (nzeta);
		
	cpx_set_ui (ess, i, 0);
	cpx_borwein_zeta (zeta, ess, prec);
	z_got = mpf_get_d(zeta[0].re);

	printf("duuude s=%d expcected=%f got=%f\n", i, z_expected, z_got);
}


int main (int argc, char * argv[])
{
	if (argc < 3)
	{
		fprintf (stderr, "Usage: %s ndigits nterms\n", argv[0]); 
		exit (1);
	}

	/* the decimal precison (number of decimal places) */
	int prec = atoi (argv[1]);
	int nterms = atoi (argv[2]);

	/* compute number of binary bits this corresponds to. */
	double v = ((double) prec) *log(10.0) / log(2.0);

	/* the variable-precision calculations are touchy about this */
	/* XXX this should be stirling's approx for binomial */ 
	int bits = (int) (v + 300);
	
	/* set the precision (number of binary bits) */
	mpf_set_default_prec (bits);
	
	test_complex_riemann_zeta (nterms, prec);
	
	return 0;
}

/* =============================== END OF FILE =========================== */

