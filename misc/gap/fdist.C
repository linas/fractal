
/*
 * fdist.C
 *
 * Distribution of the Farey Numbers on the unit interval
 *
 * Linas October 2004
 */
#include <stdio.h>
#include <stdlib.h>

#include "FareyTree.h"

main(int argc, char *argv[])
{
	int i;

	if (argc <2)
	{
		fprintf (stderr, "Usage: %s <num>\n", argv[0]);
		exit (1);
	}
	int max = atoi (argv[1]);

	FareyIterator fi;

#define BINSZ 400
	double bin[BINSZ];
	for (i=0; i<BINSZ; i++)
	{
		bin[i] = 0.0;
	}

	for (i=0; i<max; i++)
	{
		int n,d;
		fi.GetNextFarey (&n, &d);

		double x = ((double) n)/ ((double) d);
		x *= BINSZ;
		int ib = (int) x;
		bin [ib] += 1.0;
	}

	for (i=0; i<BINSZ; i++)
	{
		bin[i] /= (double) max;
		printf ("%d	%g\n", i, bin[i]);
	}

}
