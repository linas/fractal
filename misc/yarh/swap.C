
/* swap.C
 *
 * Swap terms in continued fraction
 * the earliest  yarh experiments ... 
 *
 * Linas Vepstas October 2004
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"

main (int argc, char *argv[])
{
	int i;
	ContinuedFraction f;
	// f.SetEvenize();
	
	if (argc<3)
	{
		printf ("Usage: %s <p> <q>\n", argv[0]);
		exit (1);
	}
	int p = atoi (argv[1]);
	int q = atoi (argv[2]);
	printf ("# swapping terms %d and %d\n", p,q);

	int doit = 1;
	int nmax = 2521;
	for (i=nmax; 0<=i; i--)
	{
		double x = ((double) i)/ ((double) (nmax));

#define BASIC_GRAPH
#ifdef BASIC_GRAPH
		f.SetRatio (i,nmax);
		f.SwapTerms (p,q);
		double y = f.ToReal();

#ifndef DETAIL_ATTEMPT
		if ((y<0.15) && (1==doit))
		{
			doit = 0;
			int j;
#define NDETAIL 4
			int ndet = NDETAIL*nmax + 13;
			int ioff = NDETAIL*i;
			for (j=-NDETAIL+1; j<NDETAIL; j++)
			{
				x = ((double) ioff-j)/ ((double) (ndet));
				f.SetRatio (ioff-j,ndet);
				f.SwapTerms (p,q);
				double y = f.ToReal();

				/* Cheat to make graph prettier */
				if ((x<0.1) && (y<0.1)) y=0.0;
				if ((x<0.02) && (y<0.3)) y=0.0;

				printf ("%5d	%8.6g	%8.6g\n", i, x, y);
			}
		}
		else
		{
			printf ("%5d	%8.6g	%8.6g\n", i, x, y);
			doit = 1;
		}
#else // DETAIL_ATTEMPT
		printf ("%5d	%8.6g	%8.6g\n", i, x, y);
#endif // DETAIL_ATTEMPT
#endif

// #define SHOW_SELF_SIM 1
#ifdef SHOW_SELF_SIM
		int n = 5;
		int j = 2*nmax - i;
		int k = 2*nmax;
		f.SetRatio (j,k);
		f.SwapTerms (p,q);
		double y = f.ToReal();

		f.SetRatio (-(n-1)*j+n*k, -n*j+(n+1)*k);
		f.SwapTerms (p,q);
		double z = f.ToReal();

		double d = z - y/(n*y+1.0);

		printf ("%5d	%8.6g	%8.6g	%8.6g	%8.6g\n", i, x, y,z, d);
#endif

// #define OTHER_SYM
#ifdef OTHER_SYM
		f.SetRatio (i,nmax);
		f.SwapTerms (p,q);
		double y = f.ToReal();

		f.SetRatio (i,i+nmax);
		f.SwapTerms (p,q);
		double z = f.ToReal();
		double zo = z;

		// z *= (1.0+x)*(1.0-x);
		// z *= (1.0+x)*y;
		// z *= (2.0-y)*y;

		// z = z / (1.0+z);
		// z *= (2.0-x);
		// z *= (1.0+y);

#if 1
		y *= (2.0-x)/ (2.0-2.0*x);
		double a = 2.0;
		double c;
		c = a - 1.0/a;
		double d = (a-c);
		z = (a*y) / (c*y+d);

		z *= (2.0-2.0*x)/ (2.0-x);
#endif


		printf ("%5d	%8.6g	%8.6g	%8.6g	%8.6g\n", i, x, y, z, zo);
#endif

	}
}

