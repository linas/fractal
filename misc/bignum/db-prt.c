/* 
 * db-merge.c
 *
 * Merge two file caches containing pre-computed bignum values. 
 * The value with the highest precision is kept.
 *
 * Linas Vepstas July 2006
 */

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "db-cache.h"
#include "mp_zeta.h"

int
main (int argc, char * argv[])
{
	if (argc<3)
	{
		fprintf (stderr,"Usage: %s <db> <count>", argv[0]);
		exit (1);
	}
	char * db = argv[1];
	int maxidx = atoi (argv[2]);

	printf ("printout of %s up to n=%d\n", db, maxidx);

	mpf_t val;
	mpf_init (val);
	int n;
	for (n=2; n<= maxidx; n++)
	{
		int prec = fp_cache_get (db, val, n, 10);
		printf ("%d to %d\t", n, prec);
		fp_prt("", val);
		printf ("\n");
	}
	return 0;
}
