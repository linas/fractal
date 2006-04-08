
/*
 * zero-cross.c
 *
 * Get zero crossing fits for the data.  Read in a file of oscillatory
 * data in the form (n, value).
 *
 * The result of this is that the zero-crossing second derive appears to
 * be exactly pi/2.  So we're going now.
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
	double var[4000];

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

			var[n] = atof(buf);
			// printf ("%d	%g\n", n, var[n]);
			n++;
		}
		i++;
	}

	int nmax = n;

	/* compute zero corssing by linear interpolation */
	int ncross = 0;
	double vlast = var[0];
	for (i=1; i<nmax; i++)
	{
		if (vlast*var[i] < 0.0) 
		{
			double cross = (i-1) + vlast / (vlast-var[i]);
			printf ("%d\t%20.10g\n", ncross, cross);
			ncross ++;
		}
		vlast = var[i];
	}

#if WTF
	/* numerical best fit */
	double cee_0 = 2.0;
	double cee_1 = 2.553;
	cee_1 = 13.0*M_PI / 16.0;

	for (i=0; i<nmax; i++)
	{
		double b2a = 2.0*cee_1 / M_PI + 0.5;
		double guess = -b2a + sqrt (b2a*b2a + 4.0*(((double)i)-cee_0)/M_PI);
		double x = guess;
		// guess = cee_0 + (cee_1+0.25*M_PI)*x + 0.25*M_PI*x*x;
#if 0
		double a = 0.25*M_PI;
		double b = cee_1 + 0.25*M_PI;
		double c = cee_0 - ((double)i);
		double guess = (-b + sqrt (b*b-4.0*a*c))/(2.0*a);
		double x = guess;
		// guess = a*x*x + b*x +c;
#endif

	
		printf ("%d	%g\t guess=%g\n", i, var[i], guess);;
	}

	/* now interpolate */
	double cross[100];
	int icross = 0;

	double prev = var[0];
	n=1;
	while (n < nmax)
	{
		double hur = var[n];
		if (prev*hur <= 0.0)
		{
			// printf ("%d %d %g  %g\n", icross, n, prev, hur);

			double x = hur - prev;
			x = prev / x;
			x = fabs (x);
			x += n-1;
			double guess = 0.25*M_PI*(icross)*(icross);
			guess += cee_0 + (cee_1+0.25*M_PI) * icross;
			printf ("crossings -- %d %g\t guess=%g\t err=%g\n", icross, x, guess, x-guess);
			cross[icross] = x;
			icross ++;
		}
		prev = hur;
		n++;
	} 

	/* now take diffs */
	int icmax = icross;
	for (i=1; i<icmax; i++)
	{
		double d = cross[i] - cross[i-1];
		double guess = 0.5*M_PI*i;
		guess += cee_1;
		printf ("first-order differences -- %d	%g\t\t guess=%g\t err=%g\n", i, d, guess, d-guess);
	}


	double avg = 0.0;
	int acnt = 0;
	for (i=1; i<icmax-1; i++)
	{
		double d = 2.0 * cross[i] - cross[i-1] - cross[i+1];
		printf ("second-order differences -- %d	%g\t guess=%g \terr=%g\n", i, d, 0.5*M_PI, d+0.5*M_PI);
		avg += d;
		acnt ++;
	}

	avg /= acnt;
	printf ("average value = %g\n", avg);
#endif
}
