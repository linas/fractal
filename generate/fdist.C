/*
 * fdist.C
 *
 * Distribution of the Farey Numbers on the unit interval
 * AKA the Minkowski measure or multi-fractal measure.
 *
 * Show the Mellin transform.
 * See also hardy.C for the singular Poisson kernel.
 *
 * Linas October 2004
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"
#include "FareyTree.h"
#include "brat.h"

static int nbins = 0;
static double * bin = NULL;
static double * logs = NULL;

static void make_bins(int _nbins, int depth)
{
	int i;

	nbins = _nbins;

	int max = 1 << depth;

	FareyIterator fi;

	bin = (double *) malloc (nbins * sizeof (double));
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

	ContinuedFraction f;

	/* Cache of precomputed values */
	logs = (double *) malloc (nbins * sizeof (double));
	for (int i=0; i<nbins; i++)
	{
		/* x is the midpoint of the bin */
		// double x = ((double) 2*i+1) / ((double) 2*nbins);

		/* Likewise, the midpoint */
   	f.SetRatio (2*i+1, 2*nbins);
   	double far = f.ToFarey (); 

		logs[i] = log(far);
	}
}

// Compute the Mellin transform
static void mellin_c (double re_q, double im_q, double *prep, double *pimp)
{
	/* Compute the integral of the distribution */
	double re_gral = 0.0;
	double im_gral = 0.0;
	for (int i=1; i<nbins; i++)
	{
		// double x = ((double) i) / ((double) nbins);
		// double lgx = log(x);
		double lgx = logs[i];
		double r = exp (re_q * lgx);
		double re = r * cos(im_q * lgx);
		double im = r * sin(im_q * lgx);

		re_gral += re * bin[i];
		im_gral += im * bin[i];
	}

	*prep = re_gral;
	*pimp = im_gral;
}

static double bincount_series (double re_q, double im_q, int itermax, double param)
{
	if (0 == nbins) make_bins (9111, 24);

	double rep, imp;
	mellin_c (re_q, im_q, &rep, &imp);
	// return sqrt (rep*rep+imp*imp);
	// return rep;
	return (atan2 (imp,rep)+M_PI)/(2.0*M_PI);
}

DECL_MAKE_HEIGHT(bincount_series);

/* --------------------------- END OF LIFE ------------------------- */
