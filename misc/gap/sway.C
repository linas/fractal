
/*
 * sway.C
 *
 * simple graph of swaying
 *
 * Linas Vepstas October 2004
 */

#include "Farey.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


main (int argc, char *argv[])
{
	ContinuedFraction f;
	double w, t;
	int n,d;

   if (argc <2) {
      printf ("Usage: %s <denom> \n", argv[0]);
      exit (1);
   }

	d = atoi (argv[1]);
	printf ("#\n# denom=%d\n#\n", d);
  
#ifdef SHOW_SYMMETRY
	for (n=1; n<d; n++)
	{
		t = ((double) n) / ((double) d);
		int gcf = gcf32(n,d);
		int dd = d/gcf;
		int nn = n/gcf;
   	f.SetRatio (nn,dd);

		w = 1.0;
		double xp = f.ToXEven(w);

		xp -= t;
		xp *= dd*dd;
		xp -= w - w*w + 0.5*w*w*w;
		xp += t;
		
   	f.SetRatio ((dd-nn),dd);
		double xm = f.ToXOdd(w);
		xm -= (1.0-t);
		xm *= dd*dd;
		xm -= w - w*w + 0.5*w*w*w;
		xm += (1.0-t);
		
		printf ("%8.6g	%4d	%8.6g	%8.6g	%g\n", t, dd, xm, xp, xm+xp);
	}
#endif

	int i;
	for (w=0.01; w<1.0; w+= 0.01)
	{
		n = 1; d=200;
		t = ((double) n) / ((double) d);
		int gcf = gcf32(n,d);
		int dd = d/gcf;
		int nn = n/gcf;
   	f.SetRatio (nn,dd);

		double xp = f.ToXEven(w);

		xp -= t;
		xp *= dd*dd;
		xp -= w - w*w + 0.5*w*w*w;
		xp += t;

		double a = xp-t;
		
		n = 1; d=100;
		t = ((double) n) / ((double) d);
		gcf = gcf32(n,d);
		dd = d/gcf;
		nn = n/gcf;
   	f.SetRatio (nn,dd);

		xp = f.ToXEven(w);

		xp -= t;
		xp *= dd*dd;
		xp -= w - w*w + 0.5*w*w*w;
		xp += t;

		double b = xp-t;
		printf ("%8.6g	%8.6g	%8.6g	%g\n", w, a, b, b-a);
	}

}
