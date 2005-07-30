/*
 * bound.c
 * bound games
 * 
 * Linas Vepstas July 2005
 */

#include <math.h>
#include <stdio.h>


main ()
{
	int n;
	for (n=0; n<500; n+=20)
	{
		double dbn = 1.0/exp (-4.0*sqrt (n+1));
		
		printf ("%d\t%g\n",n, dbn);
		fflush (stdout);
	}

}

