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

/* =============================== END OF FILE =========================== */

