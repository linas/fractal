/*
 * mp-hyper.h
 *
 * High-precison Hypergeometric functions, using the 
 * Gnu Multiple-precision library.
 *
 * Linas Vepstas November 2006
 */

#include <gmp.h>
#include "mp-complex.h"

/**
 * cpx_confluent -- Confluent hypergeometric function
 * Implemented using a brute-force, very simple algo, with 
 * no attempts at optimization. Also, does not assume any 
 * precomputed constants.
 */

void cpx_confluent (cpx_t em, cpx_t a, cpx_t b, cpx_t z, unsigned int prec);
