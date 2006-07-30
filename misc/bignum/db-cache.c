
/* db-cache.c
 *
 * File cache for pre-computed bignum values.
 *
 * Linas Vepstas July 2006
 */

#include <db.h>
#include <fcntl.h>
#include <limits.h>
#include <stdio.h>
#include <sys/types.h>

#include <gmp.h>

void fp_cache_put (const char * dbname, mpf_t val, int idx, int nprec)
{
	DB *db;

	db = dbopen (dbname, O_RDWR, O_CREAT|O_EXLOCK, DB_RECNO, NULL);

	if (!db)
	{
		fprintf (stderr, "Error: cannot open the cache file\n");
		return;
	}

	db->close (db);
}

