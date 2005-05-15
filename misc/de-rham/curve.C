
/* curve.C
 *
 * draw derham curves by iteration of functions
 *
 * Linas Vepstas may 2005
 */

#include <math.h>
#include <stdio.h>

#include "Farey.h"

void gee (double *x, double *y)
{
	*x = 0.5* *x;
	*y = 0.3* *y;
}
void are (double *x, double *y) 
{
	*x = 1.0- *x;
	*y = 1.0- *y;
}

void fixpt (double val, double *x, double *y)
{
	ContinuedFraction f;
	f.SetReal(val);

	int i = 0;
	int nt = f.GetNumTerms();
	for (i=0; i<nt; i++)
	{
		int k = f.GetTerm (i);
		int j;

		for (j=0; j<k; j++) 
		{
			gee (x, y);
		}
		are (x, y);
	}
} 

main (int argc, char *argv[])
{
	int i;
	int p,q;
	q  = 43;
	for (p=0; p<q; p++) 
	{
		double val = (double) p / (double) q;
		double x = 0.0;
		double y = 0.0;
		fixpt (val, &x, &y);

		printf ("%d\t%d\t%g	%g	%g\n", p,q,val,x,y);
	}
}
