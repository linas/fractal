
/*
 * gappy.C
 *
 * Games with Gaps.
 *
 * Linas Vepstas October 2004
 */

#include "Farey.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NV 1

main (int argc, char *argv[])
{
	ContinuedFraction f;
	double w, t, gap[NV];
	int n,d,iw;

   if (argc <2) {
      printf ("Usage: %s <denom> \n", argv[0]);
      exit (1);
   }

	d = atoi (argv[1]);
	printf ("#\n# denom=%d\n#\n", d);
  
	for (n=1; n<d; n++)
	{
		int gcf, rn, rd;

		int nn = n;
		int dd = d;
#if SHOW_RAND
		nn = rand() >> 18;
		dd = rand() >> 18;
		nn = nn%dd;
		if (0 == nn) nn=1;
#endif

		// nn = 3;
		// dd = 3*n+1;

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
/*
 * In other words,  
 * gap(p/q) = (w/q^2) [2 - w + w^2/2 + (w^2/2) T(p/q)]
 * where abs T(p/q) < 1
 * for w<=1
 */
			}
if ((gap[0] > 1.0) || (gap[0] < -1.0)) {
printf ("found one: p/q = %d/%d\n", rn, rd);
}
		}
		double egap = f.ToGapEven() - f.ToGapOdd() - 1.0;
		double err = 100.0 * (gap[NV-1] - egap)/gap[NV-1];

		printf ("%g	%g	%g %g	%d	%d\n", 
			t, gap[NV-1], egap, err, rn , rd);
	}
}
