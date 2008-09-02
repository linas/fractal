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

double * bin;

void bincount(int nbins, int depth)
{
	int i;

	int max = 1 << depth;
	printf ("#\n# nbins=%d   tree depth=%d\n#\n",nbins,depth);
	printf ("# Legend:\n");
	printf ("# ........ \n");
	fflush (stdout);

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

void fourier (int nbins, int freq_max)
{
	int i, n;

	double *fre = (double *) malloc (freq_max*sizeof(double));
	double *fim = (double *) malloc (freq_max*sizeof(double));

	for (n=0; n<freq_max; n++)
	{
		fre[n] = 0.0;
		fim[n] = 0.0;
	}

	/* Compute the integral of the distribution */
   ContinuedFraction f;

	for (i=0; i<nbins; i++)
	{
		/* x is the midpoint of the bin */
		double x = ((double) 2*i+1) / ((double) 2*nbins);

		/* Likewise, the midpoint */
   	f.SetRatio (2*i+1, 2*nbins);
   	double far = f.ToFarey (); 

		for (n=0; n<freq_max; n++)
		{
			double re = cos(2.0*M_PI*n*far);
			double im = sin(2.0*M_PI*n*far);
#ifdef JACOB
			fre[n] += re;
			fim[n] += im;
#else
			fre[n] += bin[i] * bin[i] * re;
			fim[n] += bin[i] * bin[i] * im;
#endif
		}
	}

	/* renormalize */
	for (n=0; n<freq_max; n++)
	{
		fre[n] /= (double) nbins;
		fim[n] /= (double) nbins;
	}

#if 0
	for (n=0; n<freq_max; n++)
	{
		printf ("%d	%8.6g	%8.6g\n", n, fre[n], fim[n]);
	}
#endif

#if 1
	/* rebin */
	int npts = 1200;
	for(i=0; i<npts; i++)
	{
		double x = ((double) i) / ((double) npts);
		double fx = 0.0;
		for (n=1; n<freq_max; n++)
		{
			fx += 2.0* fre[n] * cos(2.0*M_PI*n*x);
		}
		printf ("%8.6g	%8.6g\n", x, fx);
	}
#endif

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

