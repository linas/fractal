
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
		fprintf (stderr, "Usage: %s <nbins> <maxiter>\n", argv[0]);
		exit (1);
	}
	int nbins = atoi (argv[1]);
	int max = atoi (argv[2]);

#define BINSZ 45720
	int bin[BINSZ];
	for (i=0; i<BINSZ; i++)
	{
		bin[i] = 0;
	}

	int n, d;
	int cnt = 0;
	for (d=1; d<max; d++)
	{
		for (n=0; n<=d; n++)
		{
			int gcf = gcf32 (n,d);
			int nn = n/gcf;
			int dd = d/gcf;
#define DO_GCF_ELIM
#ifdef DO_GCF_ELIM
			if (gcf != 1) continue;
#endif

			double x = ((double) (nn*nbins))/ ((double) dd);
			int ib = (int) x;
			if (ib >= nbins) continue;
			bin [ib] ++;
			cnt ++;
// if (ib == nbins/2) { printf ("bin %d f=%d/%d\n", ib, n, d); }
// if (ib == nbins/2-1) { printf ("bin %d f=%d/%d\n", ib, n, d); }
		}
	}
#if 0
	cnt -= bin[0];
	cnt -= bin[nbins-1];
	cnt += 2*cnt/(nbins-2);
	bin[0] = cnt/(nbins-2);
	bin[nbins-1]= cnt/(nbins-2);
#endif

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
