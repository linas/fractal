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
bincount (int nbins, int pmax, double c)
{
	int i;

	printf ("#\n# bincount of projection of Cantor polynomial\n#\n");
	printf ("#\n# nbins=%d   pow=%d\n#\n",nbins,pmax);

#define BINSZ 45720
	int bin[BINSZ];
	for (i=0; i<BINSZ; i++)
	{
		bin[i] = 0;
	}

	int max = 1<<pmax;

	int n, d;
	int cnt = 0;
	for (d=0; d<max; d+=2)
	{
		int id;
	printf ("duude d=%d\n", d);
		// work out the binary digits
		int mask = 1<<(pmax-1);
		for (id=0; id < pmax; id++)
		{
			if (mask & d)	{
				printf ("hit id=%d mask=0x%x\n", id, mask);
			}
			else
			{
				printf ("mis id=%d mask=0x%x\n", id, mask);
			}
			mask >>= 1;
		}
		for (n=0; n<=d; n++)
		{
		}
	}
	printf ("# total count=%d\n", cnt);
	for (i=0; i<nbins; i++)
	{
		double bcnt = bin[i];
		bcnt /= (double) cnt;
		bcnt *= (double) nbins;
		double x = ((double) i) / ((double) nbins);
		printf ("%5d	%8.6g	%g\n", i, x, bcnt);
	}
}

main(int argc, char *argv[])
{
	int i;

	if (argc <2)
	{
		fprintf (stderr, "Usage: %s <nbins> <maxiter>\n", argv[0]);
		exit (1);
	}
	int nbins = atoi (argv[1]);
	int max = atoi (argv[2]);

	bincount (nbins, max, 0.6);
}

