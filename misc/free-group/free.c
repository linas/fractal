
/* 
 * free.c
 *
 * Freely generated Group games
 *
 * Linas Vepstas December 2004
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>


long double complex
recur (long double complex rstart, long double complex r,
       long double complex theta, int depth)
{
	int i;
	int imax = depth;
	long double rn = rstart;
	for (i=0; i<imax; i++)
	{
		long double nt = i;
		long double complex e = cexpl (I*nt*theta);
		e *= rn;
		rn *= r;

		printf ("%d	%Lg	%Lg\n", i, creall (e), cimagl (e));
	}
}

int
main (int argc, char * argv[])
{
	int i;

	long double r = 0.9L;

	if (2 > argc)  {
		fprintf (stderr, "Usage: %s <r>\n", argv[0]);
		exit (1);
	} 

	r = atof (argv[1]);
	printf ("#\n#  r = %Lg\n", r);

	long double theta = acosl (1.0L/6.0L);

	int imax = 500;
	long double rn = 1.0L;

	recur (rn, r, theta, imax);

	return 0;
}

