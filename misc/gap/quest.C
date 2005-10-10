
/* 
 * quest.C
 *
 * draw the Question mark
 *
 * Linas October 2005
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"

ContinuedFraction far;

double mink(double x)
{
	far.SetReal (x);
	x = far.ToFarey();
	return x;
}

main (int argc, char *argv[])
{
	int i;

	int nmax = 512;

	for (i=0; i<nmax; i++)
	{
		double x = i/((double)nmax);

		double y = mink(x);

		printf ("%d	%8.6g	%8.6g\n", i, x, y);
		fflush (stdout);
	}
}
