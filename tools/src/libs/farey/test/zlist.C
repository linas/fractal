
/* 
 * FUNCTION:
 * Test operation of farey number converter
 *
 * HISTORY:
 * Linas Vepstas Januery 16 1994
 */

#include "Farey.h"
#include <stdio.h>
#include <math.h>

/* ------------------------------------------------------------ */

main (argc, argv)
int argc;
char *argv[];
{
   struct Farey *f;
   double x, y, z, t;
   int i, n;
   int nume, deno;

   if (argc <3) {
      printf ("Usage: %s <number of terms> <base> \n", argv[0]);
      exit (1);
   }

   f = CreateFarey();

   n = atoi (argv[1]);
   t = atof (argv[2]);

   /* deno = 6899 * n +1; */
   deno = 65536 *n;
   for (i=0; i<n; i++){
      /* nume = 6889*i; */
      nume = 65536*i + rand();

      x = ((double) (nume))/ ((double) deno);
      /* SetReal (f, x); */

      RatioToContinuedFraction (f, nume, deno);
      y = ContinuedFractionToZReal (f, t);

/*
      RatioToContinuedFraction (f, i, n);
      z = ContinuedFractionToZReal (f, t);

      y -= z;
*/

      printf ("i %g f %g \n", x, y);
      fflush (stdout);
   }

}

/* ---------------------- END OF FILE ------------------------- */
