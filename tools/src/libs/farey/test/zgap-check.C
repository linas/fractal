
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
   int i, j, n, m, d;
   int nume, deno;

   if (argc <4) {
      printf ("Usage: %s <n> <nume> <m>\n", argv[0]);
      exit (1);
   }

   n = atoi (argv[1]);
   nume = atoi (argv[2]);

   m = atoi (argv[3]);

   for (i = 1; i <= n; i ++) 
   {
      for (j=1; j<m; j++) {
         deno = m*i+j;
         x = ((double) nume) / ((double) deno);
         f.SetRatio (nume, deno);
         v = f.ToZRealGap ();
         t = v - 1.0 / ((double) deno*deno);
         if ((1.0e-16 < t) || (-1.0e-16 > t)) {
            printf ("problem: i=%d j=%d deno=%d x=%g gap=%g diff=%g\n", i,j, deno, x, v, t);
         }
      }
      for (j=1; j<m; j++) {
         deno = m*i-j;
         x = ((double) nume) / ((double) deno);
         f.SetRatio (nume, deno);
         v = f.ToZRealGap ();
         t = v - 1.0 / ((double) deno*deno);
         if ((1.0e-16 < t) || (-1.0e-16 > t)) {
            printf ("problem: i=%d j=%d deno=%d x=%g gap=%g diff=%g\n", i,j, deno, x, v, t);
         }
      }
      fflush (stdout);
   }

}

/* ---------------------- END OF FILE ------------------------- */
