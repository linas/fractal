/*
 * brat-gapper.C
 *
 * FUNCTION:
 * Explore dstricution of gaps in the continued fraction.
 * Specifically, plot q_(n-1)/q_n ratio of denominators.
 *
 * HISTORY:
 * quick hack -- Linas Vepstas October 1989
 * modernize -- Linas Vepstas March 1996
 * more stuff -- January 2000
 * more stuff -- October 2004
 */

#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "brat.h"
#include "gcf.h"
#include "Farey.h"
#include "FareyTree.h"

/*-------------------------------------------------------------------*/
/* The gapper tries to do the distribution of the gaps in the 
 * continued fraction
 */

void 
MakeHisto (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   double	height,
   int		itermax,
   double 	renorm)
{
   int		i,j;

   int globlen = sizex*sizey;
   for (i=0; i<globlen; i++) {
      glob [i] = 0.0;
   }

	int *tot_cnt = (int *) malloc (sizex*sizeof (int));
	int *bin_cnt = (int *) malloc (sizex*sizeof (int));
	for (i=0; i<sizex; i++)
	{
		tot_cnt[i] = 0;
		bin_cnt[i] = 0;
	}
   
	int d,n;  // denom, numerator
	ContinuedFraction f;

	FareyIterator fi;

	for (d=2; d<itermax; d++)
	{
		int dd = d;
		for (n=1; n<d; n++)
		{
			int nn = n;

// #define DO_RAND
#ifdef DO_RAND
			nn = rand() >> 10;
			dd = rand() >> 10;
			nn %= dd;
			if (0 == nn) continue;
			if (0 == dd) continue;
#endif
			// fi.GetNextFarey (&nn, &dd);

			double t = (double) nn/ (double) dd;

			i = (int) ((double) sizex * t);
			tot_cnt[i] ++;
if (i>=sizex) printf ("xxxxxxxxxxxxxxxxx\n");
			if (i>=sizex) continue;
			

			int gcf = gcf32 (nn,dd);
#ifndef DO_RAND
			if (1 != gcf) continue;
#endif
			nn /= gcf;
			dd /= gcf;
			bin_cnt [i] ++;
			
			f.SetRatio (nn,dd);

#ifdef DO_THE_GAPS
			double gap = f.ToGapEven() - f.ToGapOdd() - 1.0;
#endif

#define DO_CONVERGENTS 1
#ifdef DO_CONVERGENTS
			int nt = f.GetNumTerms ();
			double qn = f.GetConvDenom (nt);
			double qnm1 = f.GetConvDenom (nt-1);

			double gap = qnm1 / qn;
			gap *= 2.0;
printf ("duude gap=%g\n", gap);
#endif

			gap = 1.0-gap;
			j = (int) (gap * (double) sizey);
if ((j>=sizey) || (0>j)) printf ("badddddd j=%d gap=%g p/q=%d/%d\n", j, gap, nn,dd);
			if (0>j) j=0;
			if (j>=sizey) j=sizey-1;

			glob [j*sizex +i] ++;
		}
	}

   /* renormalize */
	int bin_tot = 0;
   for (i=0; i<sizex; i++) 
	{
		bin_tot += bin_cnt[i];
	}
   for (i=0; i<sizex; i++) 
	{
		double r = ((double) bin_cnt[i]) / ((double) bin_tot);
		double x = ((double) i +0.5)/((double) sizex);
		// printf ("%g	%g\n", x, r);
	}
   for (i=0; i<sizex; i++) 
	{
		if (bin_cnt[i])
		{
			double r = 1.0 / ((double) bin_cnt[i]);
			for (j=0; j<sizey; j++)
			{
				glob [j*sizex+i] *= r;
			}
		}
		// printf ("duude i=%d b=%d t=%d\n", i, bin_cnt[i], tot_cnt[i]);
   }

#if 0
   /* draw parabola */
   for (i=0; i<sizex; i++) 
	{
		double x = (double) i;
		x /= sizex;
		double y = x-0.5;
		y = 4.0*y*y;
		y *= sizey;
		j = y;
		j = sizey - j;
		glob [j*sizex+i] =15400;
	}
   for (i=0; i<sizex; i++) 
	{
		double x = (double) i;
		x += 0.5;
		x /= sizex;
		double y = x-0.5;
		y = 4.0*y*y;
		y *= sizey;
		j = y;
		j = sizey - j;
		glob [j*sizex+i] =15400;
	}
#endif

	/* Profile */
	double norm = 0.0;
	for (j=0; j<sizey; j++)
	{
		double tmp =0.0;
		for (i=0; i<sizex; i++)
		{
			tmp += glob [j*sizex+i];
		}
		tmp /= sizex;
		double y = ((double) j +0.5)/((double) sizey);
		printf ("%g	%g\n", y, tmp);
		norm += tmp;
	}
	printf ("# duude norm=%g\n", norm);

   free (bin_cnt);
   free (tot_cnt);
}

/* --------------------------- END OF LIFE ------------------------- */
