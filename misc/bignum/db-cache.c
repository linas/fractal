
/* db-cache.c
 *
 * File cache for pre-computed bignum values.
 *
 * Linas Vepstas July 2006
 */

#include <db_185.h>
#include <fcntl.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>

#include <gmp.h>
#include "db-cache.h"

void fp_cache_put (const char * dbname, mpf_t val, int idx, int nprec)
{
	DB *db;

	db = dbopen (dbname, O_RDWR, O_CREAT, DB_RECNO, NULL);

	if (!db)
	{
		fprintf (stderr, "Error: cannot open the cache file\n");
		return;
	}

	DBT pkey, vkey;
	char buf[50];
	sprintf (buf, "val[%d]", idx);
	vkey.data = &buf;
	vkey.size = strlen(buf)+1;

	DBT vdat;
	vdat.data = "asdf";
	vdat.size = 5;

	db->put (db, &vkey, &vdat, R_SETCURSOR);

	db->close (db);
}

