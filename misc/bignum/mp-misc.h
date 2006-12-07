/*
 * mp-misc.h
 *
 * High-precison misc functions, using the 
 * Gnu Multiple-precision library.
 *
 * Linas Vepstas July 2005
 */

#include <gmp.h>
#include "mp-complex.h"

void i_prt (char * str, mpz_t val);
void fp_prt (char * str, mpf_t val);
void cpx_prt (char * str, const cpx_t const val);

/**
 * get_epsilon - return 10^{-prec} 
 */
void get_epsilon (mpf_t eps, int prec);

/* prec is the decimal precison (number of decimal places) */
/* nterms is the number of an's to compute */
void set_bits (int prec, int nterms);
