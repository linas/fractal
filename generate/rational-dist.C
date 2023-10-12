/*
 * rational-dist.C
 *
 * FUNCTION:
 * Distribution of rational numbers in bins. See also fdist.C
 * for the Minowski measure 
 *
 * HISTORY:
 * first version of rational stuff Linas October 2004
 * more stuff -- September 2005
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

int
bincount (int bin[], int nbins, int max)
{
	int i;

	// printf ("#\n# bincount of rationals using plain math\n#\n");
	// printf ("#\n# nbins=%d   maxiter=%d\n#\n",nbins,max);

	for (i=0; i<nbins; i++)
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
			int gcf = gcf64 (n,d);
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

	return cnt;
}

/*-------------------------------------------------------------------*/
/* This routine fills in the interior of the a square with rational
 * number distribution counts.
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
   int		i,j, globlen;
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 0.0;

#define BINSZ 4570
	int bin[BINSZ];
	
   for (i=0; i<sizey; i++) 
	{
      if (i%10==0) printf(" start row %d\n", i);
		
		int nbins = i+2;
		int cnt = bincount (bin, nbins, 10*nbins);

		double norm = ((double) nbins) / ((double) cnt);
		
		double delta = ((double) nbins) / ((double) sizex);
		double xbin = 0.0;
      for (j=0; j<sizex; j++) 
		{
			xbin += delta;
			int ibin = (int) xbin;
         glob [i*sizex +j] = norm * bin[ibin];
      }
   }
}

/* --------------------------- END OF LIFE ------------------------- */
