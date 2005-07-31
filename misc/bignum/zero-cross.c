
/*
 * zero-cross.c
 *
 * get zero crossing fits for the data
 *
 * Linas July 2005
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

main ()
{
	int c;
	char buf[4000];
	double var[300];

	int disc = 1;	
	int i =0;
	int n=0;
	while( (c=getchar()) != EOF)
	{
		if (disc && c != '\t') continue;
		disc = 0;
		if (c == '\t') continue;

		buf[i] = c;

		if ( c == '\n')
		{
			buf[i] = 0;
			disc = 1;
		   i=-1;	

			var[n] = atof(buf);
			n++;
			printf ("%g\n", var[n-1]);
		}
		i++;
	}
}
