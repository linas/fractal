/*
 * rdist.C
 *
 * Distribution of the rationals
 * converted to use gnu multiple precision libs
 *
 * Linas October 2004
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"

void
bincount (int nbins, int max)
{
	int i;

	printf ("#\n# bincount of rationals using plain math\n#\n");
	printf ("#\n# nbins=%d   maxiter=%d\n#\n",nbins,max);

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
// #define DO_GCF_ELIM
#ifdef DO_GCF_ELIM
			int gcf = gcf32 (n,d);
			int nn = n/gcf;
			int dd = d/gcf;
			if (gcf != 1) continue;
#else 
			int nn = n;
			int dd = d;
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

void
gmp_bincount (int nbins, int max)
{
	int i;

	printf ("#\n# bincount of rationals using gnu multi-precison\n#\n");
	printf ("#\n# nbins=%d   maxiter=%d\n#\n",nbins,max);

// #define DO_GCF_ELIM
#ifdef DO_GCF_ELIM
	printf ("# common factors were eliminated\n");
#else
	printf ("# NOOOOO common factors eliminatation\n");
#endif 

#define BINSZ 45720
	int bin[BINSZ];
	for (i=0; i<BINSZ; i++)
	{
		bin[i] = 0;
	}

	mpz_t gib, gn, gnbins;
	mpz_init (gib);
	mpz_init (gn);
	mpz_init (gnbins);
	mpz_set_ui (gnbins, nbins);

	unsigned int n, d;
	int cnt = 0;
	for (d=1; d<max; d++)
	{
		for (n=0; n<=d; n++)
		{
// #define DO_GCF_ELIM
#ifdef DO_GCF_ELIM
			int gcf = gcf32 (n,d);
			int nn = n/gcf;
			int dd = d/gcf;
			if (gcf != 1) continue;
#endif

			// implement the following bining in gmp:
			// double x = ((double) n)/ ((double) d);
			// x *= nbins;
			// int ib = (int) x;
	
			mpz_mul_ui (gn, gnbins, n);
			mpz_fdiv_q_ui (gib, gn, d);
	
			unsigned int ib = mpz_get_ui (gib);

			/* Note the following discard of the p/q= 1 bin affects 
			 * the normalization!  */
			if (ib >= nbins) continue;
			bin [ib] ++;
			cnt ++;
		}
	}
#if 0
	cnt -= bin[0];
	cnt -= bin[nbins-1];
	cnt += 2*cnt/(nbins-2);
	bin[0] = cnt/(nbins-2);
	bin[nbins-1]= cnt/(nbins-2);
#endif

	double renorm = ((double) nbins) / ((double) cnt);
	printf ("# total count=%d    renorm=%g\n", cnt, renorm);
	for (i=0; i<nbins; i++)
	{
		double bcnt = bin[i];
		bcnt *= renorm;
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

	// bincount (nbins, max);
	gmp_bincount (nbins, max);
}

