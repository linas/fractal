/* 
 * perimeter.c
 *
 * Map tree to perimeter of the circle.  Used to sanity check
 * various results.
 *
 * Linas Vepstas May 2007
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "flt.h"

/**
 * real_to_mobius -- convertt real number to mobius xform
 *
 * Innput is real number x, assume 0<x<1.  Perform binary 
 * expansion of x. Use L,R for binary digits 0,1.  Return 
 * the result.
 */
mobius_t real_to_mobius (double x)
{
	if (0.0>x) x -= (int) x - 1;
	if (1.0<x) x -= (int) x;

	mobius_t ell = mobius_set (1,0,1,1);
	mobius_t are = mobius_set (1,1,0,1);

	mobius_t acc = mobius_set (1,0,0,1);
	int i;
	for (i=0; i<45; i++)
	{
		if (0.5 <= x)
		{
			acc = mobius_mul (acc, are);
			x -= 0.5;
		}
		else
		{
			acc = mobius_mul (acc, ell);
		}
		x *= 2.0;
	}

	return acc;
}

main()
{
	real_to_mobius (0.7333);
}
