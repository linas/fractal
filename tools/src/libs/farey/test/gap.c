
/* 
 * FUNCTION:
 * Test operation of farey number converter
 *
 * HISTORY:
 * Linas Vepstas Januery 16 1994
 */

#include "Farey.h"
#include "Prime.h"
#include <stdio.h>
#include <math.h>

/* ------------------------------------------------------------ */

main (argc, argv)
int argc;
char *argv[];
{
   struct Farey *f;
   struct Prime *p;
   double r, x, y, z, t;
   double v, vprev, gap;
   int i, j, n;
   int nume, deno;

   if (argc <3) {
      printf ("Usage: %s <number of terms> <base> \n", argv[0]);
      exit (1);
   }

   f = CreateFarey();
   p = CreatePrime ();

   n = atoi (argv[1]);
   z = atof (argv[2]);

   r = 0.0;
   v = 0.0;
   gap = 0.0;
   j = 0;
   for (i = 0; i < n; i ++){
      deno = GetPrime (p, i);
      for (nume = 1; nume < deno; nume++){


   x = ((double) nume) / ((double) deno);
   RatioToContinuedFraction (f, nume, deno);
   v = ContinuedFractionToEFareyGap (f, z);
   r += v;
   j ++;
/*   x = (double) j; */

         printf ("i %g f %g \n", x, r);
         fflush (stdout);
      }
   }

}

/* ---------------------- END OF FILE ------------------------- */
