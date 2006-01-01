
/*
 * f-wavelet.C
 *
 * Obtain binaryu coefficieints corresponding to 
 * Distribution of the Farey Numbers on the unit interval
 * This may or may not require wavelet-like analysis.
 *
 * Linas October 2004
 */
#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"
#include "FareyTree.h"

void GetNextDyadic (unsigned int *n, unsigned int *d)
{
	static unsigned int last_d = 1;
	static unsigned int last_n = 1;

	last_n += 2;
	if (last_n > last_d)
	{
		last_d *= 2;
		last_n = 1;
	}

	*n = last_n;
	*d = last_d;
}

/* bincount the farey fractions */
double * 
bincount(int npow, int oversamp)
{
	int i;
	int nbins = 1<<npow;
	int max = nbins * oversamp;

	printf ("#\n# nbins=%d   maxiter=%d\n#\n", nbins, max);

	FareyIterator fi;

	double *bin = (double *) malloc ((nbins+1)*sizeof (double));
	for (i=0; i<nbins; i++)
	{
		bin[i] = 0.0;
	}
	bin[0] = 1.0;
	bin[nbins-1] = 1.0;

	/* Compute the distribution by bining */
	int cnt =2;
	for (i=0; i<max-1; i++)
	{
		int n,d;
		fi.GetNextFarey (&n, &d);
		// GetNextDyadic (&n, &d);
		// printf ("duude i=%d n/d= %d/%d\n", i, n,d);

		double x = ((double) n)/ ((double) d);
		x *= nbins;
		int ib = (int) x;

		// Is d a power of two ? if it is, then split its
		// contribution over two neighboring bins
		int dd = d;
		while (dd%2 == 0 && dd>1)
		{
			dd >>= 1;
		}
		if (dd == 1)
		{
			bin [ib-1] += 0.5;
			bin [ib] += 0.5;
		}
		else
		{
			bin [ib] += 1.0;
		}
		cnt ++;
	}

#if 0
	/* Compute the integral of the distribution */
   ContinuedFraction f;
	double gral = 0.0;
	for (i=0; i<nbins; i++)
	{
		double bcnt = bin[i];
		bcnt /= (double) cnt;
		gral += bcnt;
		bcnt *= nbins;
		double x = ((double) i) / ((double) nbins);

   	f.SetRatio (2*i+1, 2*nbins);
   	double far = f.ToFarey (); 

		printf ("%6d	%8.6g	%8.6g	%8.6g	%8.6g\n", i, x, bcnt, gral, far);
	}
#endif

	return bin;
}

/* perform fourier transform; However, this explcitly omits
 * non-scaling wavelengths: only power of two are considered.
 */

void fourier (double *bins, int npow)
{
	int i;
	int p;
	int nbins = 1<<npow;

	int step = nbins >>2;
	for (p=1; p<npow; p++)
	{
		double aleft = 0.0;
		double aright =  0.0;
		int ncnt = 0;
		while (ncnt < nbins-1)
		{
			for (i=0; i< step; i++)
			{
				aleft += bins[ncnt+i];
			}
			ncnt += step;
			for (i=0; i< step; i++)
			{
				aright += bins[ncnt+i];
			}
			ncnt += step;

			for (i=0; i< step; i++)
			{
				aright += bins[ncnt+i];
			}
			ncnt += step;

			for (i=0; i< step; i++)
			{
				aleft += bins[ncnt+i];
			}
			ncnt += step;
		}
		double lpr = aleft+aright;
		double lmr = aleft-aright;
		// printf ("duude p=%d step=%d \t(l,r)=( %g\t%g )\n", p, step, aleft, aright);
		printf ("duude p=%d step=%d \t%g\t%g\t%g\t%g\n", p, step, lpr, lmr, lmr/lpr, lpr/lmr);
		step >>= 1;
	}
}

main(int argc, char *argv[])
{
	int i;

	if (argc <2)
	{
		fprintf (stderr, "Usage: %s <npow> <oversamp>\n", argv[0]);
		exit (1);
	}
	int npow = atoll (argv[1]);
	int oversamp = atoll (argv[2]);

	double * bins = bincount (npow, oversamp);
	fourier (bins, npow);
}

