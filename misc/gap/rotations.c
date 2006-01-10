
/*
 * rotations.c
 *
 * Try to build the question mark function out of 
 * hyperbolic roations.
 *
 * Linbas Vepstas January 2006
 */

#include <stdio.h>

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

main ()
{
	int i;
	int imax = 100;
	for (i=0; i<imax; i++)
	{
		double x = i / ((double) imax);
		double y = rot_right(x);
		printf ("%d	%g	%g\n", i, x, y);
	}
}
