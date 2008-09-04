/*
 * jdist.C
 *
 * Distribution of the Farey Numbers on the unit interval
 *
 * Linas October 2004
 */
#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"

void bincount(int nbins, int depth)
{
	int i;

   ContinuedFraction f;

	int max = 1 << depth;
	printf ("#\n# nbins=%d   depth=%d\n#\n",nbins,depth);
	printf ("# Legend:\n");
	printf ("# i, x, bin_cnt, bin_cnt_sum, exact_farey\n");
	fflush (stdout);

	int *bin = (int *) malloc (nbins * sizeof (int));
	for (i=0; i<nbins; i++)
	{
		bin[i] = 0;
	}
	bin[0] = 0;
	bin[nbins-1] = 0;

	/* Compute the distribution by bining */
	int cnt = 0;
	for (i=0; i<max; i++)
	{
		int n = rand();

   	f.SetRatio (n, RAND_MAX + 1);
   	double x = f.ToFarey (); 
		x *= nbins;
		int ib = (int) x;
		bin [ib] ++;
		cnt ++;
	}

	/* Compute the integral of the distribution */
	double gral = 0.0;
	double egral = 0.0;
	for (i=0; i<nbins; i++)
	{
		/* gral is the ordinary integral of the bin count */
		double rect = bin[i] / ((double) cnt);
		gral += rect;

		/* Be careful to bin-count the entropy 
		 * (use discrete not continuous formula) 
		 */
		double entropy = - rect * log(rect);
		if (bin[i] == 0) entropy = 0;
		egral += entropy;

		double x = ((double) i) / ((double) nbins);

#if 1
		printf ("%6d	%8.6g	%8.6g	%8.6g	%8.6g %8.6g\n", 
			i, x, rect, gral, entropy, egral);
		fflush (stdout);
#endif
	}

	printf ("#Total entropy for ?' =  %18.16g\n", egral);
}

main(int argc, char *argv[])
{
	int i;

	if (argc <3)
	{
		fprintf (stderr, "Usage: %s <nbins> <tree-depth>\n", argv[0]);
		exit (1);
	}
	int nbins = atoi (argv[1]);
	int depth = atoi (argv[2]);

	bincount (nbins, depth);
}

