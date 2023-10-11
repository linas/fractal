
/*
 * beta.c
 *
 * play funny binomial game on zeta terms 
 * Linas November 2004
 */

#include <math.h>
#include <stdio.h>
#include "ache.h"
#include "zetafn.h"

int
main ()
{
	int n;

	for (n=0; n<20; n++)
	{
		double an = a_sub_n (n);
		double tn = an + 0.5/((double) (n+1));
		// double rad = tn - 1.0 + M_GAMMA;
		double rad = tn + M_GAMMA;

		double bn = 1.0 + pow (rad, 1.0/((double) n));
		printf ("its %d \tan=%12.10g \trad=%12.10g \tbn=%12.10g\n", n, an, rad, bn);
	}

	return 0;
}
