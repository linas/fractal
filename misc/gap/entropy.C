
/* entropy.C
 *
 * computation of the entropy, per shanon-type question 
 *  about info content.
 *
 * Linas Vepstas October 2004
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"

main (int argc, char *argv[])
{
	double z;
	int i;
	ContinuedFraction f;
	// f.SetEvenize();

	int qmax = 3531;

	qmax = atoi (argv[1]);

	for (i=1; i<qmax; i++)
	{

		int p = i;
		int q = qmax;

#if 0
		p = rand();
		q = rand();
		p %= q;
#endif

		int gcf = gcf32 (p,q);
		p /= gcf;
		q /= gcf;

		double x = ((double) p)/ ((double) q);
		f.SetRatio (p,q);
		int n = f.GetNumTerms();
		
		double y= log (q) / ((double) n);;
		
		printf("%5d	%8.6g	%8.6g\n", i,x,y);
	}
}

