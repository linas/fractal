
/* 
 * takagi.C
 *
 * draw the Takagi curve
 *
 * Linas October 2004
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double triangle (double x)
{
	double t = x - floor(x);
	if (0.5 > t) return t;
	return 1.0-t;
}

double takagi (double w, double x)
{
	int k;
	double acc = 0.0;
	double tw = w;
	double tp = 1.0;
	for (k=0; k<1000; k++)
	{
		acc += tw* triangle (tp*x);
		tp *= 2.0;
		tw *= w;
		if (1.0e-14 > tw) break;
	}

	return acc;
}

main (int argc, char *argv[])
{
	int i;

	int nmax = 531;

	if (argc <2)
	{
		printf ("Usage: %s <w-value>\n", argv[0]);
		exit (1);
	}
	double w = atof(argv[1]);

	for (i=0; i<nmax; i++)
	{
		double x = i/((double)nmax);
		double tw = takagi (w, x);

		printf ("%d	%8.6g	%8.6g\n", i, x, tw);
	}
}
