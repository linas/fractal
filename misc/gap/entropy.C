
/* entropy.C
 *
 * computation of the entropy, per shanon-type question 
 *  about info content.
 *
 * Linas Vepstas October 2004
 */

#include <math.h>
#include <stdio.h>

#include "Farey.h"

main (int argc, char *argv[])
{
	double z;
	int i;
	ContinuedFraction f;
	// f.SetEvenize();

	int qmax = 531;
	for (i=0; i<qmax; i++)
	{

		int p = i;
		int q = qmax;
		double x = ((double) p)/ ((double) q);
		f.SetRatio (p,q);
		int n = f.GetNumTerms();
		
		// double y= log (q) / ((double) n);;
		double y= 1.0 / ((double) n);;
		
		printf("%5d	%8.6g	%8.6g\n", i,x,y);
	}
}

