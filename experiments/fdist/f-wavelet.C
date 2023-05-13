
/*
 * f-wavelet.C
 *
 * Obtain binary coefficieints corresponding to 
 * Distribution of the Farey Numbers on the unit interval
 * This may or may not require wavelet-like analysis.
 *
 * derived from fdist.C orginally
 *
 * Linas October 2004
 */
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
	int max = nbins * (1<<oversamp);

	printf ("#\n# nbins=%d   maxiter=%d\n#\n", nbins, max);

	FareyIterator fi;

	double *bin = (double *) malloc ((nbins+1)*sizeof (double));
	for (i=0; i<nbins; i++)
	{
		bin[i] = 0.0;
	}
	bin[0] = 1.0;
	bin[nbins-1] = 1.0;

	/* Compute the distribution by histogramming */
	int cnt =2;
	for (i=1; i<max; i++)
	{
		int n,d;
		fi.GetNextFarey (&n, &d);
		// GetNextDyadic (&n, &d);
		// printf ("duude i=%d n/d= %d/%d\n", i, n,d);

		long double x = ((long double) n)/ ((long double) d);
		x *= (long double) nbins;
		int ib = (int) floorl(x);

		// Is d a power of two ? if it is, then split its
		// contribution over two neighboring bins
		// make sure it really is on a bin boundary
		int dd = d;
		int np = npow;
		while (dd%2 == 0 && np >0 && dd>1)
		{
			dd >>= 1;
			np --;
		}
		if (dd == 1)
		{
			bin [ib-1] += 0.5;
			bin [ib] += 0.5;
			// printf ("split i=%d n/d= %d/%d	ibin=%d\n", i, n,d, ib);
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

	int step = nbins >>1;
	for (p=0; p<npow; p++)
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

			double lpr = aleft+aright;
			double lmr = aleft-aright;
			printf ("p=%d step=%d/%d \ttot= %8.1f\tl=%8.1f (%6.4f) r=%8.1f (%6.4f) D=%8.1f (%6.5f)\t%g\n", 
				p, ncnt/(2*step), nbins/(2*step), lpr, aleft, aleft/lpr, aright, aright/lpr, lmr, lmr/lpr, lpr/lmr);
			aleft = 0.0;
			aright = 0.0;
#if 0
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
#endif
		}
		printf ("\n");

		double lpr = aleft+aright;
		double lmr = aleft-aright;
		// printf ("duude p=%d step=%d \t(l,r)=( %g\t%g )\n", p, step, aleft, aright);
		// printf ("duude p=%d step=%d \t%g\t%g\t%g\t%g\n", p, step, lpr, lmr, lmr/lpr, lpr/lmr);
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

