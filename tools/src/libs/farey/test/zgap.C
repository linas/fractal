
/* 
 * FUNCTION:
 * Test operation of farey number converter
 *
 * HISTORY:
 * Linas Vepstas Januery 16 1994
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"

/* ------------------------------------------------------------ */

main (int argc, char *argv[])
{
   ContinuedFraction f;
   double r, x, y, z, t;
   double v, vprev, gap;
   int i, n;
   int nume, deno;

   if (argc <2) {
      printf ("Usage: %s <denom> \n", argv[0]);
      exit (1);
   }

   n = atoi (argv[1]);
   deno = 2*n;

   for (i = 0; i <= n; i ++) 
   {

      nume = 1;
      deno = i+2;
      x = ((double) nume) / ((double) deno);
      f.SetRatio (nume, deno);
      v = f.ToZRealGap ();
      // z = (f.ToZReal(1.001) - f.ToZReal(1.0)) / 0.001;
      // f.Print ();

      printf ("%g	%g	%g\n", x, 4.0*v, 2.0 /((double) deno*deno));
      fflush (stdout);
   }

}

/* ---------------------- END OF FILE ------------------------- */
