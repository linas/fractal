/*
 * hardy.C
 *
 * Poison integral of the Minkowski measure (i.e. of the
 * distribution of the Farey Numbers on the unit interval)
 * The Poison ingtegral takes the form of a so-called 
 * "singular inner function" in the theory of Hardy spaces.
 *
 * See directory fractal/misc/fdist for more details.
 *
 * Linas October 2004
 * Linas October 2008
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"
#include "FareyTree.h"

double *bin, *si, *co;

void bincount(int nbins, int depth)
{
	int i;
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

	/* renormalize */
	for (i=0; i<nbins; i++)
	{
		bin[i] *= ((double) nbins) / ((double) cnt);
	}
}

void init(int nbins)
{
	si = (double *) malloc (nbins * sizeof (double));
	co = (double *) malloc (nbins * sizeof (double));

	for (i=0; i<nbins; i++)
	{
		/* x is the midpoint of the bin */
		double x = ((double) 2*i+1) / ((double) 2*nbins);

		/* Likewise, the midpoint */
   	f.SetRatio (2*i+1, 2*nbins);
   	double far = f.ToFarey (); 

		si[i] = sin(2.0*M_PI*far);
		co[i] = cos(2.0*M_PI*far);
	}
}

void hardy(double re, double im, double *reh, double *imh)
{
	int i, n;

	/* Compute the integral of the distribution */

	double resum = 0.0;
	double imsum = 0.0;
	for (i=0; i<nbins; i++)
	{
		double red =  co[i] - re;
		double imd =  si[i] - im;

		double deno = red*red + imd*imd;

		double ren = co[i] + re;
		double imn = si[i] + im;

		double regr = ren*red + imn*imd;
		double imgr = imn*red - imd*ren;

		regr /= deno;
		imgr /= deno;

		resum += regr * bin[i];
		imsum += imgr * bin[i];
	}

	/* renormalize */
	resum /= (double) nbins;
	imsum /= (double) nbins;

	*reh = resum;
	*imh = imsum;
}

static double hardy_series (double re_q, double im_q, int itermax, double param)
{
	double rep, imp;
	rep = re_q;
	imp = im_q;
	hardy (re_q, im_q, &rep, &imp);
	
	return sqrt (rep*rep+imp*imp);
	// return rep;
	// return (atan2 (imp,rep)+M_PI)/(2.0*M_PI);
}

DECL_MAKE_HEIGHT(hardy_series);

/* --------------------------- END OF LIFE ------------------------- */
