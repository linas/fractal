
/* 
 * free.c
 *
 * Freely generated Group games
 *
 * Linas Vepstas December 2004
 */

#include <math.h>
#include <stdlib.h>
#include <complex.h>


int
main (int argc, char * argv[])
{
	int i;

	long double r = 0.9L;

	if (1 >argc)  {
		fprintf (stderr, "Usage: %s <r>\n", argv[0]);
		exit (1);
	} 

	long double theta = acosl (1.0L,3.0L);

	int imax = 500;
	long double rn = 1.0L;
	for (i=0; i<imax; i++)
	{
		long double nt = i;
		long double complex e = cexpl (I*nt*theta);
		e *= rn;
		rn *= r;

		printf ("%d	%Lg	%Lg\n", i, creall (e), cimagl (e));
	}

	return 0;
}

