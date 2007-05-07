
/* feigen.C
 *
 * eigen-functions:
 * bernoulli-map eigenfunctions concated with Minkowski map.
 * Basically borken right now.
 *
 * Linas Vepstas October 2004
 */

#include <math.h>
#include <stdio.h>

#include "Farey.h"
#include "question.h"

main (int argc, char *argv[])
{
	int i;
	ContinuedFraction f;
	f.SetEvenize();


#define NP 6
	double arr[NP];

	int nmax = 400;
	for (i=0; i<nmax; i++)
	{
		double x = ((double) 2*i+1)/ ((double) (2*nmax));

		f.SetRatio (2*i+1, 2*nmax);
		double xb = f.ToFarey();

		int j;
		double yb = xb;
		for (j=0; j<NP; j++)
		{
			yb *= 0.5;
			double y = question_inverse (yb);
			arr[j] = y;
		}
		printf ("%5d	%8.6g", i, x);
		for (j=0; j<NP; j++)
		{
			printf ("	%8.6g", arr[j]);
		}
		printf ("\n");

		fflush (stdout);
	}
}

