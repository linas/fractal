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
bincount (int nbins, int pmax, double z)
{
	int i;

	printf ("#\n# bincount of projection of Cantor polynomial\n#\n");
	printf ("#\n# nbins=%d   pow=%d\n#\n",nbins,pmax);

#define BINSZ 45720
	double bin[BINSZ];
	for (i=0; i<BINSZ; i++)
	{
		bin[i] = 0.0;
	}

	int max = 1<<pmax;

	int n, d;
	double cnt = 0.0;
	for (d=0; d<max; d+=2)
	{
		int id;
		// work out the binary digit expansion of d / pmax
		int mask = 1<<(pmax-1);
		double clo = 0.0;
		double chi = 0.0;
		double zn = 1.0;
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

		clo *= 1.0-z;
		chi *= 1.0-z;

		// clo and chi are the enpoints of the cantor intervals
		// i.e. the stuff that's left after midpoint removal.
		// printf ("Cantor interval=%d/%d lo=%g hi=%g\n", d/2, max/2, clo, chi);

		// now bincount
		clo *= nbins;
		chi *= nbins;
		double nlo = ceil (clo);
		double nhi = floor (chi);

		int ilo = (int) nlo;
		int ihi = (int) nhi;

		// if the interval is spread across multiple bins ... 
		// or if the interval fits entirely within one bin
		if (ilo-1 == ihi)
		{
			// printf ("onebin (%g, %g) i=(%d, %d), \tdelt= %g\n", clo, chi, ilo, ihi, chi-clo);
			bin [ihi] += chi - clo;
			cnt += chi - clo;
		}
		else
		{ 
			if (ilo-1 >= 0) {
				bin [ilo-1] += nlo - clo;
			}
			bin [ihi] += chi - nhi;

			// printf ("duude (%g, %g) i=(%d, %d), \tdelt= %g %g\n", clo, chi, ilo, ihi, nlo-clo, chi-nhi);

			cnt += (nlo-clo) + (chi-nhi);
	
			for (n=ilo; n<ihi; n++)
			{
				bin[n] += 1.0;
				cnt += 1.0;
			}
		}
	}
	double measure = cnt / ((double) nbins);
	printf ("# total count=%g measure=%g\n", cnt, measure);
	for (i=0; i<nbins; i++)
	{
		double bcnt = bin[i];
		bcnt /= (double) cnt;
		bcnt *= (double) nbins;
		double x = ((double) 2*i+1) / ((double) 2*nbins);
		printf ("%5d	%8.6g	%g\n", i, x, bcnt);
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
	double z = atof (argv[3]);

	bincount (nbins, max, z);
}

