
/*
 * interval.c
 *
 * replot ont an interval
 *
 * Linas Vepstas August 2009
 */

#include <math.h>
#include <stdio.h>

#include "moebius.h"

int main()
{
	int p, k;
	int pm = 8;
	int n = 1;
	for (p=1; p<pm; p++)
	{
		int m = 1<<p;
		for (k=1; k<m; k+=2)
		{
			double x = ((double) k) / ((double) m);
			double y = divisor(n);

			printf("%d	%f	%f\n", n, x, y);
			n++;
		}
	} 

	return 0;
}
