
/*
 * tent-map.C
 *
 * dyadic symbolic values of the tent map.
 * iterate the unit-height tent map. The result is a set of
 * symbolic values, left or right. use these symbolic values 
 * to construct a binary real number. Graphs this number.
 *
 * Linas June 2005
 */

#include <math.h>
#include <stdio.h>

long double tent_dyadic (long double x, long double left, long double right)
{
	int i;

	long double b = 0.0L;
	long double norm = 0.5L / (right-left);
	left *= norm;
	right *= norm;

	x -= floorl (x);
	for (i=0; i<30; i++)
	{
		/* convert dyadic symbol to binary number */
		if (x < 0.5L) 
		{
			b += left;
		}
		else
		{
			b += right;
		}
		left *= norm;
		right *= norm;

		// printf ("duude i=%d x=%Lg b=%Lg \n", i, x, b);
		/* iterate on the tent map */
		x *= 2.0L;
		if (x >= 1.0L) x = 2.0L-x;
	}

	return b;
}

main () 
{
	int i;

	int nmax = 723;
	for (i=0; i<nmax; i++)
	{
		long double x = ((long double) i) / ((long double) nmax);

		long double y = tent_dyadic (x, 0.0, 1.0);
		// long double y = tent_dyadic (x, -1.0, 1.0);

		printf ("%d	%Lg	%Lg\n", i, x,y);
	}
}
