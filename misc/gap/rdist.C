
/*
 * rdist.C
 *
 * Distribution of the rationals
 *
 * Linas October 2004
 */
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"

main(int argc, char *argv[])
{
	int i;

	if (argc <2)
	{
		fprintf (stderr, "Usage: %s <num>\n", argv[0]);
		exit (1);
	}
	int max = atoi (argv[1]);

#define BINSZ 400
	double bin[BINSZ];
	for (i=0; i<BINSZ; i++)
	{
		bin[i] = 0.0;
	}

	int n, d;
	for (d=1; d<max; d++)
	{
		for (n=1; n<d; n++)
		{

			double x = ((double) n)/ ((double) d);
			x *= BINSZ;
			int ib = (int) x;
			bin [ib] += 1.0;
		}
	}

	for (i=0; i<BINSZ; i++)
	{
		bin[i] /= (double) max;
		double x = ((double) i) / ((double) BINSZ);
		printf ("%8.6g	%g\n", x, bin[i]);
	}
}
