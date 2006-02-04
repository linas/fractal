/* 
 * asubn.c
 *
 * Graph of a_sub_n -- for manuscript
 * but also explore the gneral case better
 *
 * Linas Vepstas 2005
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

	double prev = 0.0;
	for (i=1; i<40; i++)
	{
#if TRADITIONAL_A_SUB_N
		double x = a_sub_n (i);
		double r = sqrt (i+1);
		r = exp (-4.0*r);
		x /= r;
		printf ("%d	%8.6g\n", i, x);
#endif
#if 0
		double x = t_sub_n (i, a);
		double y = a_sub_n (i);
		x -= 1.0/(2.0*a*((double)(i+1)));
		printf ("%d	%8.6g	%8.6g\n", i, x, y);
#endif
#define HURWITZ_A_SUB_N
#ifdef HURWITZ_A_SUB_N
		double x = a_sub_n (i);
		double y = hz_a_sub_n (i, a);
		double z = 1.0/(y-prev);
		// z /= i*(i+1);
		printf ("%d	%8.6g	%8.6g	%8.6g\n", i, x, y, z);
		prev = y;
#endif
	}
	return 0;
}
