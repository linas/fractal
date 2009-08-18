
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
	int pm = 11;
	int n = 1;
	for (p=1; p<pm; p++)
	{
		int m = 1<<p;
		for (k=1; k<m; k+=2)
		{
			double x = ((double) k) / ((double) m);
			double y = divisor(n);
			// double y = moebius_mu(n);
			// double y = liouville_omega(n);
			// double y = liouville_lambda(n);
			// double y = mangoldt_lambda(n);

			printf("%f	%f\n", x, y);
			n++;
		}
	} 

	return 0;
}
