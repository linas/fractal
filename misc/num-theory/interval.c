
/*
 * interval.c
 *
 * replot ont an interval
 *
 * Linas Vepstas August 2009
 */

#include <math.h>
#include <stdio.h>

#include "modular.h"

main()
{
	int pm = 8;
	int n = 1;
	for (int p=1; p<pm; p++)
	{
		m = 1<<p;
		for (int k=1; k<m; k+=2)
		{
			double x = ((double) k) / ((double) m);
			double y = divisor(n);

			printf("%d	%f	%f\n", n, x, y);
			n++;
		}
	} 
}
