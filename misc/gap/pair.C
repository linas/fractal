
/*
 * pair.C
 *
 * Games with final pair.
 *
 * Linas Vepstas November 2004
 */

#include "Farey.h"
#include "gcf.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


void
prt_pair (int nn, int dd)
{
	ContinuedFraction f;

	int gcf, rn, rd;

	gcf = gcf32 (nn,dd);
	rn = nn/gcf;
	rd = dd/gcf;

	double t = ((double) nn) / ((double) dd);
  	f.SetRatio (rn,rd);
	double xp = f.ToLastPair();
	double at = 2.0 *atan (xp)/M_PI;

	printf ("%5d	%5d	%8.6g	%8.6g	%8.6g\n", rn, rd, t, xp, at);
}

main (int argc, char *argv[])
{
	int n,d;

   if (argc <2) {
      printf ("Usage: %s <denom> \n", argv[0]);
      exit (1);
   }

	d = atoi (argv[1]);
	printf ("#\n# denom=%d\n#\n", d);
  
	int nn = 1;
	int dd = d;
	for (n=1; n<d; n++)
	{
// #define SHOW_RAND
#ifdef SHOW_RAND
		nn = rand() >> 18;
		dd = rand() >> 18;
		nn = nn%dd;
		if (0 == nn) continue;
		prt_pair (nn,dd);
#endif

		prt_pair (n,d);
	}
}
