/* 
 * db-cache.h
 *
 * File cache for pre-computed bignum values.
 *
 * Linas Vepstas July 2006
 */

#include <gmp.h>

void fp_cache_put (const char * dbname, mpf_t val, int idx, int nprec);
