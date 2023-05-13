/*
 * exacto.C
 *
 * explore exact relations
 *
 * Linas Vepstas October 2004
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"

main(int argc, char * argv[]) 
{
	ContinuedFraction f;

	int nn = 23;
	int dd = 53;

	f.SetRatio (nn,dd);

	int nterms = f.GetNumTerms();

	int q = f.GetConvDenom (nterms);
	int qm = f.GetConvDenom (nterms-1);

	printf ("# have nterms=%d q=%d qq=%d\n", nterms, q, qm);

	double w;
	for (w=0.02; w<1.01; w+=0.02)
	{
		double xm = f.ToXMinus(w);
		xm -= ((double) nn) / ((double) dd);
		xm *= ((double) dd*dd);
		
		double gm = ((double) qm)/((double) q);
		gm = 1.0-gm;
		gm *= w;
		gm += 1.0;
		gm = w/gm;
		gm = -gm;

		double xp = f.ToXPlus(w);
		xp -= ((double) nn) / ((double) dd);
		xp *= ((double) dd*dd);
		
		double gp = ((double) qm)/((double) q);
		gp *= w;
		gp += 1.0;
		gp = w/gp;

		printf ("%8.6g	%8.6g	%8.6g	%8.6g\n", w, xp, gp, xp-gp);
	}
}
