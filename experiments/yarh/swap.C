
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

	// int nmax = 2521;  //  works very well for s12
	int nmax = 30031;
	int a1 = 1;
	int a2 = 1;
	for (i=nmax; 0<=i; i--)
	{
		double x = ((double) i)/ ((double) (nmax));

// #define BASIC_GRAPH
#ifdef BASIC_GRAPH
		f.SetRatio (i,nmax);
		f.SwapTerms (p,q);
		double y = f.ToReal();

// #define CLEANUP_S12
#ifdef CLEANUP_S12
		double clamp = 1.0 / ((double) a1);
		if (x<clamp)
		{
			printf ("%5d	%8.6g	%8.6g\n", i, clamp, 0.0);
			y = ((double) a1-1) / ((double) a1);
			printf ("%5d	%8.6g	%8.6g\n", i, clamp, y);
			a1 = (int) floor(1.0/x) + 1;
		}
#endif /* CLEANUP_S12 */

// #define CLEANUP_S13
#ifdef CLEANUP_S13
		double clamp = 1.0 / ((double) a1);
		if (x<clamp)
		{
			printf ("%5d	%8.6g	%8.6g\n", i, clamp, 0.0);
			printf ("%5d	%8.6g	%8.6g\n", i, clamp, 1.0);
			a1 = (int) floor(1.0/x) + 1;
			a2 = 1;
		}
#endif /* CLEANUP_S13 */

		printf ("%5d	%8.6g	%8.6g\n", i, x, y);
#endif /* BASIC_SWAP */

#define MOBIUS_GRAPH
#ifdef MOBIUS_GRAPH
		f.SetRatio (i,nmax);
		f.MobTerms (p,q, 1, 0, 1, 1);
		double y = f.ToReal();

		printf ("%5d	%8.6g	%8.6g\n", i, x, y);
#endif /* MOBIUS_SWAP */

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

