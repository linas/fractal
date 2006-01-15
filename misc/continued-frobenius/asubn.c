/* 
 * asubn.c
 *
 *  Graph of a_sub_n -- for manuscript
 * but also explore the gneral case better
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ache.h"

int
main (int argc, char * argv[])
{
	int i;

	if (2> argc)
	{
		fprintf (stderr, "Usage: %s <a-val>\n", argv[0]);
		exit (1);
	}
	long double a = atof (argv[1]);

	for (i=1; i<40; i++)
	{
#if TRADITIONAL_A_SUB_N
		double x = a_sub_n (i);
		double r = sqrt (i+1);
		r = exp (-4.0*r);
		x /= r;
		printf ("%d	%8.6g\n", i, x);
#endif
#if GEN
		double x = t_sub_n (i, a);
		x -= -1.0+3.0/(2.0*a*((double)(i+1)));
		x = -x;
		printf ("%d	%8.6g\n", i, x);
#endif
#if 1
		double x = a_sub_n (i);
		double y = a_sub_n_poor_convergence(i);
		printf ("%d	%8.6g	%8.6g\n", i, x, y);
#endif
	}
	return 0;
}
