
/* 
 * FUNCTION:
 * print out a bunch of points for a fairy graph
 *
 * HISTORY:
 * Linas Vepstas January 16 1994
 */

#include "Farey.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ------------------------------------------------------------ */

main (int argc, char *argv[])
{
   ContinuedFraction f;
   double x, y, t;
   int n, d;

   if (argc <2) {
      printf ("Usage: %s <denom> \n", argv[0]);
      exit (1);
   }

   d = atoi (argv[1]);

	for (n=0; n<d; n++)
	{
   	t = ((double) n) / ((double) d);
   	/* f.SetReal (t); */
   	f.SetRatio (n, d);
   	x = f.ToFarey (); 
		y = f.ToPAdicFarey(3);
		printf ("%g	%g	%g\n", t, x, y);	
	}

#ifdef SHOW_SYMMETRY

   f.SetRatio (n, n+d);

   y = 2.0 * f.ToFarey (); 
#endif

   // printf ("---------- f(%d/%d) = %f %f\n", n,d,x, y); 

}

/* ---------------------- END OF FILE ------------------------- */

