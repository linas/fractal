/* 
 * db-cache.h
 *
 * File cache for pre-computed bignum values. Stores values,
 * with indicated precision, so we don't waste too much compute 
 * time.
 *
 * You may redistribute this, and/or modify this, under the terms 
 * of the GNU Lesser General Public License as  published by the
 * Free Software Foundation; either version 2.1 of the License,
 * or (at your option) any later version.
 * 
 * Copyright Linas Vepstas July 2006
 */

#include <gmp.h>

/**
 * fp_cache_put -- put mpf_t value into the database.
 * @prec: number of decimal places of precision to store.
 * @idx:  index under which to store the value
 */
void fp_cache_put (const char * dbname, const mpf_t val, int idx, int nprec);

/**
 * fp_cache_get -- get mpf_t from database
 * Returns 0 if no value in the database, or if the value in the
 * database has fewer than nprec digits. Thus, nprec is a minimum
 * requirement.
 */ 
int fp_cache_get (const char * dbname, mpf_t val, int idx, int nprec);
