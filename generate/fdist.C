/*
 * fdist.C
 *
 * Distribution of the Farey Numbers on the unit interval
 * AKA the Minkowski measure or multi-fractal measure.
 *
 * Linas October 2004
 */
#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"
#include "FareyTree.h"
#include "brat.h"

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

int nbins = 0;
double * bin = NULL;

void make_bins(int _nbins, int depth)
{
	int i;

	nbins = _nbins;

	int max = 1 << depth;
	printf ("#\n# nbins=%d   tree depth=%d\n#\n",nbins,depth);
	printf ("# Legend:\n");
	printf ("# i, x, bin_cnt, bin_cnt_sum, exact_farey\n");
	fflush (stdout);

	FareyIterator fi;

	double *bin = (double *) malloc (nbins * sizeof (double));
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

	/* now, renormalize */
	for (i=0; i<nbins; i++)
	{
		bin[i] /= ((double) cnt);
	}
}


static void bincount_c (double re_q, double im_q, double *prep, double *pimp)
{
	double tmp;

	double rcz = cos (re_q) * cosh (im_q);
	double icz = sin (re_q) * sinh (im_q);

	/* Compute the integral of the distribution */
	double gral = 0.0;
	for (i=0; i<nbins; i++)
	{
		/* gral is the ordinary integral of the bin count */
		gral += bin[i];

		double x = ((double) i) / ((double) nbins);

	}
}

	*prep = rcz;
	*pimp = icz;
}

static double bincount_series (double re_q, double im_q, int itermax, double param)
{
	double rep, imp;
	bincount_c (re_q, im_q, &rep, &imp);
	// return sqrt (rep*rep+imp*imp);
	// return rep;
	return (atan2 (imp,rep)+M_PI)/(2.0*M_PI);
}

DECL_MAKE_HEIGHT(bincount_series);

/* --------------------------- END OF LIFE ------------------------- */
