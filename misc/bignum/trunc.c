/* 
 * trunc.c
 *
 * truncate long string values to something that
 * gnuplot can read, e.g.about 18 decimal places
 *
 * Linas Vepstas July 2005
 */

#include <math.h>
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
		if (c == '#')
		{
			printf ("%c", c);
			while( (c=getchar()) != '\n')
			{
				printf ("%c", c);
			}
			printf("\n");
			continue;
		}
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

			// double corr = exp (4.0* (sqrt(n+1)-sqrt(n)));
			double corr = exp (-4.0* (sqrt(n+1)));
			corr *= exp (3.602* (sqrt(n+1)));
			var *= corr;
			var /= sin(M_PI*(-2.125+sqrt(2.125*2.125+4.0*(n-2)/M_PI)));
		 	printf ("%d\t%23.18g\n", n, var);

			// double x = log (n);
			// var = log (var);
			// printf ("%d\t%23.18g\t%23.18g\n", n, x, var);
			n++;
		}
		i++;
	}
}

