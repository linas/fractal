
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
   int i, n, m, d;
   int nume, deno;

   if (argc <4) {
      printf ("Usage: %s <n> <nume> <m> <a>\n", argv[0]);
      exit (1);
   }

   n = atoi (argv[1]);

   nume = atoi (argv[2]);
   m = atoi (argv[3]);
   d = atoi (argv[4]);

   for (i = 1; i <= n; i ++) 
   {
      deno = m*i+d;
      x = ((double) nume) / ((double) deno);
      f.SetRatio (nume, deno);
      v = f.ToZRealGap ();
      // z = (f.ToZReal(1.001) - f.ToZReal(1.0)) / 0.001;
      // f.Print ();
      t = v - 1.0 / ((double) deno*deno);

      if ((1.0e-16 < t) || (-1.0e-16 > t)) {
         printf ("problem: i=%d deno=%d x=%g gap=%g diff=%g\n", i, deno, x, v, t);
      }
      // printf ("%d	%g	%g	%g\n", deno, x, v, t);
      fflush (stdout);
   }

   printf ("excellent!\n");
}

/* ---------------------- END OF FILE ------------------------- */

