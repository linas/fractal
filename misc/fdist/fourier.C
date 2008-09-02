/*
 * fourier.C
 * 
 * Fourier transform of
 * Distribution of the Farey Numbers on the unit interval
 *
 * Linas October 2004
 * Linas Sept 2008
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"
#include "FareyTree.h"

void bincount(int nbins, int depth)
{
	int i;

	int max = 1 << depth;
	printf ("#\n# nbins=%d   tree depth=%d\n#\n",nbins,depth);
	printf ("# Legend:\n");
	printf ("# i, x, bin_cnt, bin_cnt_sum, exact_farey\n");
	fflush (stdout);

	FareyIterator fi;

	bin = (int *) malloc (nbins * sizeof (int));
	for (i=0; i<nbins; i++)
	{
		bin[i] = 0;
	}
	bin[0] = 1;
	bin[nbins-1] = 1;

	/* Compute the distribution by bining */
	int cnt =2;
	for (i=0; i<max; i++)
	{
		int n,d;
		fi.GetNextFarey (&n, &d);
		// GetNextDyadic (&n, &d);

		double x = ((double) n)/ ((double) d);
		x *= nbins;
		int ib = (int) x;
		bin [ib] ++;
		cnt ++;
	}

	/* Compute the integral of the distribution */
   ContinuedFraction f;
	double gral = 0.0;
	double egral = 0.0;
	double dgral = 0.0;
	double fprev = 0.0;
	for (i=0; i<nbins; i++)
	{
		/* gral is the ordinary integral of the bin count */
		double rect = bin[i] / ((double) cnt);
		gral += rect;

		/* Be careful to bin-count the entropy 
		 * (use discrete not continuous formula) */
		double entropy = - rect * log(rect);
		if (bin[i] == 0) entropy = 0;
		egral += entropy;

		double x = ((double) i) / ((double) nbins);

   	f.SetRatio (2*i+1, 2*nbins);
   	double far = f.ToFarey (); 

		/* Integral of the jacobian */
		double delt = (far - fprev)*nbins;

		if (1.0e-8 < delt) dgral += rect / delt;

#if 1
		printf ("%6d	%8.6g	%8.6g	%8.6g	%8.6g	%8.6g	%8.6g	%8.6g\n", 
			i, x, bcnt, gral, far, entropy, dgral, egral);
		fflush (stdout);
#endif
		fprev = far;
	}

	printf ("#Total entropy =  %18.16g\n", egral);
}

main(int argc, char *argv[])
{
	int i;

	if (argc <3)
	{
		fprintf (stderr, "Usage: %s <nbins> <tree-depth> <arg>\n", argv[0]);
		exit (1);
	}
	int nbins = atoi (argv[1]);
	int depth = atoi (argv[2]);
	double misc_arg = atof (argv[3]);

	bincount (nbins, depth);
}

