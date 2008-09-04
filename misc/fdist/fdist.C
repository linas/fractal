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


void bincount(int nbins, int depth)
{
	int i;

	int max = 1 << depth;
	printf ("#\n# nbins=%d   tree depth=%d\n#\n",nbins,depth);
	printf ("# Legend:\n");
	printf ("# i, x, bin_cnt, bin_cnt_sum, exact_farey\n");
	fflush (stdout);

	FareyIterator fi;

	int *bin = (int *) malloc (nbins * sizeof (int));
	for (i=0; i<nbins; i++)
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
	double egral = 0.0;
	double ejgral = 0.0;
	double dgral = 0.0;
	double fprev = 0.0;
	for (i=0; i<nbins; i++)
	{
		/* gral is the ordinary integral of the bin count */
		double rect = bin[i] / ((double) cnt);
		gral += rect;

		/* Be careful to bin-count the entropy 
		 * (use discrete not continuous formula) 
		 * This is the entropy of just ?'(x) and not of the jacobian */
		double entropy = - rect * log(rect);
		if (bin[i] == 0) entropy = 0;
		egral += entropy;

		double bcnt = rect * nbins;

		double x = ((double) i) / ((double) nbins);

   	f.SetRatio (i+1, nbins);
   	double far = f.ToFarey (); 

		/* Integral of the jacobian */
		double delt = (far - fprev) * nbins;
		// if (0.0 != delt)
		if (1.0e-8 < delt)
		{
			dgral += rect / delt;
			ejgral += entropy / delt;
		}

#if 1
		printf ("%6d	%8.6g	%8.6g	%8.6g	%8.6g	%8.6g	%8.6g	%8.6g %8.6g\n", 
			i, x, bcnt, gral, far, entropy, dgral, egral, ejgral);
		fflush (stdout);
#endif
		fprev = far;
	}

	printf ("#Total entropy for ?' =  %18.16g (log2 = %18.16g)\n", egral, egral * log(2.0));
	printf ("#Total entropy for ?'(?^{-1}) =  %18.16g (log2=%18.16g)\n", ejgral, ejgral*log(2.0));
}

void 
gmp_bincount(int nbins, int depth)
{
	int i;

	int max = 1 << depth;
	printf ("#\n# nbins=%d   tree depth=%d\n#\n",nbins,depth);

	int *bin = (int *) malloc (nbins * sizeof (int));
	for (i=0; i<nbins; i++)
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

	if (argc <3)
	{
		fprintf (stderr, "Usage: %s <nbins> <tree-depth> [<arg>]\n", argv[0]);
		exit (1);
	}
	int nbins = atoi (argv[1]);
	int depth = atoi (argv[2]);
	if (4 == argc)
	{
		double misc_arg = atof (argv[3]);
	}

	bincount (nbins, depth);
	// gmp_bincount (nbins, depth);
}

