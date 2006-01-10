
/*
 * rotations.c
 *
 * Try to build the question mark function out of 
 * hyperbolic roations.
 *
 * Linbas Vepstas January 2006
 */

#include <stdio.h>

/* hyperbolic rotatins to the left and to the right */
double rot_left (double x)
{
	if (x<0.25) return 2.0*x;
	if (x<0.5) return x+0.25;
	return 0.5*x+0.5;
}

double rot_right (double x)
{
	if (x<0.5) return 0.5*x;
	if (x<0.75) return x-0.25;
	return 2.0*(x-0.5);
}

/* right joined to left */
double scurve (double x)
{
	if (x<0.5) return 0.5 * rot_right (2.0*x);
	return 0.5 + 0.5 * rot_left (2.0*x-1.0);
}

/* recursive s-curve */
double rcurve (double x, int cnt)
{
	int i;
	if (1 == cnt) return scurve (x);

	for (i=1; i<cnt; i++)
	{
		if (x<0.5) {
			x = 0.5 * scurve (2.0*x);
		} else {
			x = 0.5 + 0.5 * scurve (2.0*x-1.0);
		}
	}

	return x;
}


main ()
{
	int i;
	int imax = 100;
	for (i=0; i<imax; i++)
	{
		double x = i / ((double) imax);
		double y = rcurve(x, 3);
		printf ("%d	%g	%g\n", i, x, y);
	}
}
