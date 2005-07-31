/* 
 * trunc.c
 *
 * truncate long string values to something that
 * gnuplot can read, e.g.about 18 decimal places
 *
 * Linas Vepstas July 2005
 */

#include <stdio.h>
#include <stdlib.h>

main ()
{
	int c;
	char buf[4000];

	/* read in floating point values in the first column */
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

			double var = atof(buf);
			printf ("%d\t%23.18g\n", n, var);
			n++;
		}
		i++;
	}
}

