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

int
main (int argc, char * argv[])
{
	int maxidx=10;
	char * dbina;
	char * dbinb;
	char * dbout;

	if (argc<5)
	{
		fprintf (stderr,"Usage: %s <out-db> <in-dba> <in-dbb> <count>\n", argv[0]);
		exit (1);
	}
	dbout = argv[1];
	dbina = argv[2];
	dbinb = argv[3];
	maxidx = atoi (argv[4]);

	printf ("Merging %s and %s into %s for n=%d\n", dbina, dbinb, dbout, maxidx);

	mpf_t vala, valb;
	mpf_init (vala);
	mpf_init (valb);
	int n;
	for (n=2; n<= maxidx; n++)
	{
		int preca = fp_cache_get (dbina, vala, n, 10);
		int precb = fp_cache_get (dbinb, valb, n, 10);

		if (0 >= preca && 0 >= precb) continue;
		
		if (precb < preca)
		{
			fp_cache_put (dbout, vala, n, preca);
			printf ("%d to %d from %s\n", n, preca, dbina);
		}
		else
		{
			fp_cache_put (dbout, valb, n, precb);
			printf ("%d to %d from %s\n", n, precb, dbinb);
		}
	}
	return 0;
}
