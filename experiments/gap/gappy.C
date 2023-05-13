
/*
 * gappy.C
 *
 * Games with Gaps.
 *
 * In other words,  
 * gap(p/q)(w) = (w/q^2) [2 - w + w^2/2 + (w^2/2) T(p/q)]
 * where abs T(p/q)(w) <= 1
 * for w<=1
 * and where 0<= T(p/q)(0) <= 1  when w=0
 *
 * Linas Vepstas October 2004
 */

#include "Farey.h"
#include "gcf.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NV 1


void
prt_gap (int nn, int dd)
{
	ContinuedFraction f;
	double w, t, gap[NV];
	int iw;

	int gcf, rn, rd;

	gcf = gcf32 (nn,dd);
	rn = nn/gcf;
	rd = dd/gcf;

	t = ((double) nn) / ((double) dd);
  	f.SetRatio (rn,rd);

	for (iw=0; iw<NV; iw++)
	{
		w = (double) (iw+1);
		w /= (double) NV;
		w *= 0.1;
		// w *= 100.0;
		double xm = f.ToXOdd(w);
		double xp = f.ToXEven(w);
		gap[iw] = xp - xm;
		gap[iw] *= (double) rd*rd;
		gap[iw] /= w;
		if (1.1 > w) {
			gap[iw] -= 2.0;
			gap[iw] += w;
			gap[iw] -= 0.5*w*w;
			// biggest remaining variation is 0<~ gap < 0.5 *w*w
			// although for larger w, it should be -0.5*w*w <= gap < 0.5 *w*w
			gap[iw] /= 0.5*w*w;

		}

		/* This should never happen */
		if ((gap[0] > 1.0) || (gap[0] < -1.0)) {
			fprintf (stderr, "xxxxxxxxxxxxxxxxxxxxxxxxxxx found one: p/q = %d/%d\n", rn, rd);
			exit(1);
		}
	}
	double egap = f.ToGapEven() - f.ToGapOdd() - 1.0;
	double err = 100.0 * (gap[NV-1] - egap)/gap[NV-1];

	printf ("%5d	%5d	%8.6g	%8.6g	%g %g\n", rn, rd, t, egap, gap[NV-1], err);
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
		prt_gap (nn,dd);
#endif

// #define PARABOLA
#ifdef PARABOLA
		prt_gap (1, n+1);
		prt_gap (n, n+1);
		prt_gap (n,d);
#endif

		prt_gap (2,2*n+1);
		prt_gap (3,2*n+3);
		prt_gap (5,2*n+5);
		
	}
}
