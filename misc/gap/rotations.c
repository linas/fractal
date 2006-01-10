
/*
 * rotations.c
 *
 * Try to build the question mark function out of 
 * hyperbolic roations.
 *
 * Linbas Vepstas January 2006
 */

#include <stdio.h>
#include <stdlib.h>

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
	if (0 == cnt) return scurve (x);

	if (x<0.5) {
		x = 0.5 * rcurve (2.0*x, cnt-1);
	} else {
		x = 0.5 + 0.5 * rcurve (2.0*x-1.0, cnt-1);
	}

	return x;
}

double qcurve (double x, int imax)
{
	int i;
	for (i=0; i<imax; i++)
	{
		x = rcurve (x, imax - i -1);
	}
	return x;
}


main (int argc, char *argv[])
{
	int i;

	int ir = atoi (argv[1]);

	int imax = 400;
	for (i=0; i<imax; i++)
	{
		double x = i / ((double) imax);
		double y = qcurve(x, ir);
		printf ("%d	%g	%g\n", i, x, y);
	}
}
