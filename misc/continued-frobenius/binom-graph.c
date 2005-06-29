
/*
 * graph of binomial coeffcients on real number line 
 *
 * Linas June 2005
 */

#include <math.h>
#include <stdio.h>
#include "zetafn.h"

int
main () 
{
	int p,q;

	for (q=1; q<25; q++)
	{
		for (p=1; p<q; p++)  
		{
			double x = (double) p / (double) q;

			double y = binomial (p,q);
			printf ("%g	%g\n", x, y);
		}
	}

	return 0;
}
