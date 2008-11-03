/*
 * hardy.C
 *
 * Poison integral of the Minkowski measure (i.e. of the
 * distribution of the Farey Numbers on the unit interval)
 * The Poison ingtegral takes the form of a so-called 
 * "singular inner function" in the theory of Hardy spaces.
 *
 * Linas October 2004
 * Linas October 2008
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"
#include "FareyTree.h"

double *bin, *si, *co;

void bincount(int nbins, int depth)
{
	int i;
	int max = 1 << depth;

	FareyIterator fi;

	bin = (double *) malloc (nbins * sizeof (double));
	for (i=0; i<nbins; i++)
	{
		bin[i] = 0.0;
	}
	bin[0] = 1.0;
	bin[nbins-1] = 1.0;

	/* Compute the distribution by bining */
	int cnt = 2;
	for (i=0; i<max; i++)
	{
		int n,d;
		fi.GetNextFarey (&n, &d);
		// GetNextDyadic (&n, &d);

		double x = ((double) n)/ ((double) d);
		x *= nbins;
		int ib = (int) x;
		bin [ib] += 1.0;
		cnt ++;
	}

	/* renormalize */
	for (i=0; i<nbins; i++)
	{
		bin[i] *= ((double) nbins) / ((double) cnt);
	}
}

void init(int nbins)
{
	si = (double *) malloc (nbins * sizeof (double));
	co = (double *) malloc (nbins * sizeof (double));

	for (i=0; i<nbins; i++)
	{
		/* x is the midpoint of the bin */
		double x = ((double) 2*i+1) / ((double) 2*nbins);

		/* Likewise, the midpoint */
   	f.SetRatio (2*i+1, 2*nbins);
   	double far = f.ToFarey (); 

		si[i] = sin(2.0*M_PI*far);
		co[i] = cos(2.0*M_PI*far);
	}
}

void hardy(double re, double im)
{
	int i, n;

	/* Compute the integral of the distribution */

	for (i=0; i<nbins; i++)
	{
		/* x is the midpoint of the bin */
		double x = ((double) 2*i+1) / ((double) 2*nbins);

	}

	/* renormalize */
	fre[n] /= (double) nbins;
	fim[n] /= (double) nbins;

}

main(int argc, char *argv[])
{
	int i;

	if (argc < 4)
	{
		fprintf (stderr, "Usage: %s <nbins> <tree-depth> <freq>\n", argv[0]);
		exit (1);
	}
	int nbins = atoi (argv[1]);
	int depth = atoi (argv[2]);
	int max_freq = atoi (argv[3]);

	bincount (nbins, depth);
	fourier (nbins, max_freq);
}

