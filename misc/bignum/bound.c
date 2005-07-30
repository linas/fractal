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
	for (n=1; n<500; n+=20)
	{
		double ref = exp (-4.0*sqrt (n+1));
		double dbn = exp (-4.0*sqrt (n));
		
		double x = n;
		double y = x/log(x);

		printf ("%d\t%g	%g	%g\n",n, ref, dbn, y );
		fflush (stdout);
	}

}

