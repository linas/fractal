/*
 * cantor-dist.C
 *
 * Distribution of the cantor function
 *
 * Linas October 2004
 * Linas December 2005
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void
bincount (int nbins, int pmax, long double z)
{
	int i;

	printf ("#\n# bincount of projection of Cantor polynomial\n#\n");
	printf ("#\n# nbins=%d   pow=%d	z=%Lg\n#\n",nbins,pmax, z);

#define BINSZ 45720
	long double bin[BINSZ];
	for (i=0; i<BINSZ; i++)
	{
		bin[i] = 0.0;
	}

	int max = 1<<pmax;

	int n, d;
	long double cnt = 0.0;
	for (d=0; d<max; d+=2)
	{
		int id;
		// work out the binary digit expansion of d / pmax
		int mask = 1<<(pmax-1);
		long double clo = 0.0;
		long double chi = 0.0;
		long double zn = 1.0;
		for (id=0; id < pmax-1; id++)
		{
			if (mask & d)	
			{
				clo += zn;
			}
			mask >>= 1;
			zn *= z;
		}
		chi = clo + zn/(1-z);

		clo *= 1.0L-z;
		chi *= 1.0L-z;

		// clo and chi are the enpoints of the cantor intervals
		// i.e. the stuff that's left after midpoint removal.
		// printf ("Cantor interval=%d/%d lo=%g hi=%g\n", d/2, max/2, clo, chi);

		// now bincount
		clo *= nbins;
		chi *= nbins;
		long double nlo = ceill (clo);
		long double nhi = floorl (chi);

		int ilo = (int) nlo;
		int ihi = (int) nhi;

		// if the interval is spread across multiple bins ... 
		// or if the interval fits entirely within one bin
		if (ilo-1 == ihi)
		{
			// printf ("onebin (%g, %g) i=(%d, %d), \tdelt= %Lg\n", clo, chi, ilo, ihi, chi-clo);
			bin [ihi] += chi - clo;
			cnt += chi - clo;
		}
		else
		{ 
			if (ilo-1 >= 0) {
				bin [ilo-1] += nlo - clo;
			}
			bin [ihi] += chi - nhi;

			// printf ("duude (%g, %g) i=(%d, %d), \tdelt= %Lg %Lg\n", clo, chi, ilo, ihi, nlo-clo, chi-nhi);

			cnt += (nlo-clo) + (chi-nhi);
	
			for (n=ilo; n<ihi; n++)
			{
				bin[n] += 1.0;
				cnt += 1.0;
			}
		}
	}
	long double measure = cnt / ((long double) nbins);
	printf ("# total count=%Lg measure=%Lg\n#\n", cnt, measure);
	for (i=0; i<nbins; i++)
	{
		long double bcnt = bin[i];
		bcnt /= (long double) cnt;
		bcnt *= (long double) nbins;
		long double x = ((long double) 2*i+1) / ((long double) 2*nbins);
		printf ("%5d	%8.6Lg	%Lg\n", i, x, bcnt);
	}
}

main(int argc, char *argv[])
{
	int i;

	if (argc <3)
	{
		fprintf (stderr, "Usage: %s <nbins> <pow> <z>\n", argv[0]);
		exit (1);
	}
	int nbins = atoi (argv[1]);
	int max = atoi (argv[2]);
	long double z = atof (argv[3]);

	bincount (nbins, max, z);
}

