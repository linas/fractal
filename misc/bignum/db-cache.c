
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
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <gmp.h>
#include "db-cache.h"

void fp_cache_put (const char * dbname, mpf_t val, int idx, int nprec)
{
	DB *db;

	db = dbopen (dbname, O_RDWR|O_CREAT, 
	             S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH, DB_HASH, NULL);

	if (!db)
	{
		fprintf (stderr, "Error: cannot open the cache file\n");
		return;
	}

	/* Save the value data to the file */
	DBT vkey;
	char buf[50];
	sprintf (buf, "val[%d]", idx);
	vkey.data = &buf;
	vkey.size = strlen(buf)+1;

	/* the printf is floating point */
	size_t allo = nprec*sizeof(char)+50;
	char *vstr= malloc (allo);
	gmp_snprintf (vstr, allo, "%.*Fg", nprec, val);

	DBT vdat;
	vdat.data = vstr;
	vdat.size = strlen(vstr)+1;

	db->put (db, &vkey, &vdat, 0);
	free (vstr);

	/* Save the precision data to the file */
	DBT pkey;
	sprintf (buf, "prec[%d]", idx);
	pkey.data = &buf;
	pkey.size = strlen(buf)+1;

	DBT pdat;
	pdat.data = &nprec;
	pdat.size = sizeof(int);

	db->put (db, &pkey, &pdat, 0);

	db->close (db);
}

int fp_cache_get (const char * dbname, mpf_t val, int idx, int nprec)
{
	DB *db;

	db = dbopen (dbname, O_RDONLY, 0, DB_HASH, NULL);

	if (!db) return 0;

	/* Look for the stored precisio of data in the file */
	DBT pkey;
	char buf[50];
	sprintf (buf, "prec[%d]", idx);
	pkey.data = &buf;
	pkey.size = strlen(buf)+1;

	DBT pdat;

	/* if not key, then nothing */
	int rc = db->get (db, &pkey, &pdat, 0);
	if (rc) return 0;

	/* check to see if we have enough precision */
	int have_prec = *((int *)pdat.data);
	if (nprec > have_prec) return 0;

	/* Get the value data from the file */
	DBT vkey;
	sprintf (buf, "val[%d]", idx);
	vkey.data = &buf;
	vkey.size = strlen(buf)+1;

	rc = db->get (db, &vkey, &pdat, 0);
	if (rc) return 0;

	// printf ("found nprec=%d for %d = %s\n", have_prec, idx, pdat.data);
	mpf_set_str (val, pdat.data, 10);
	
	db->close (db);
	return have_prec;
}

/* ========================== END OF FILE ============= */
