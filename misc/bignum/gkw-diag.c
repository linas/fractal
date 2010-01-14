
/* 
 * explore and fit gkw along diagonal
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mp-gkw.h"

int main () 
{
	int n;
	mpf_t acc;
	double g;

	int prec = 400;

	/* Set the precision (number of binary bits) = prec*log(10)/log(2) */
	mpf_set_default_prec (3.3*prec);

	mpf_init (acc);

	for (n=0; n<100; n++)
	{
		gkw(acc, n, n, prec);

		g = mpf_get_d (acc);
		printf("%d	%g\n", n, g);
	}

	return 0;
}
