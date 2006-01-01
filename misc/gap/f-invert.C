
/*
 * fdist.C
 *
 * Distribution of the Farey Numbers on the unit interval
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


void bincount(int nbins, int max)
{
	int i;

	printf ("#\n# nbins=%d   maxiter=%d\n#\n",nbins,max);

	FareyIterator fi;

#define BINSZ 45400
	int bin[BINSZ];
	for (i=0; i<BINSZ; i++)
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
}

void 
gmp_bincount(int nbins, int max)
{
	int i;

	printf ("#\n# nbins=%d   maxiter=%d\n#\n",nbins,max);

#define BINSZ 45400
	int bin[BINSZ];
	for (i=0; i<BINSZ; i++)
	{
		bin[i] = 0;
	}
	bin[0] = 1;
	bin[nbins-1] = 1;

	mpz_t gib, gn, gnbins;
	mpz_init (gib);
	mpz_init (gn);
	mpz_init (gnbins);
	mpz_set_ui (gnbins, nbins);

	/* Compute the distribution by bining */
	unsigned int cnt =2;
	for (i=0; i<max; i++)
	{
		unsigned int n,d;
		GetNextDyadic (&n, &d);

		// implement the following bining in gmp:
		// double x = ((double) n)/ ((double) d);
		// x *= nbins;
		// int ib = (int) x;

		mpz_mul_ui (gn, gnbins, n);
		mpz_fdiv_q_ui (gib, gn, d);

		unsigned int ib = mpz_get_ui (gib);

		bin [ib] ++;
		cnt ++;
	}

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

