
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
   double x, y;
   int i, n;

   if (argc <2) {
      printf ("Usage: %s <number> \n", argv[0]);
      exit (1);
   }

   f = CreateFarey ();

   x = atof (argv[1]);

   n = 20;

   for (i=0; i<n; i++){
      SetReal (f, x);
      x = GetFarey (f);

      x *= 2.0;
      x -= (double) ((int) x);

      printf ("i %d f %g \n", i, x);
      fflush (stdout);
   }

}

/* ---------------------- END OF FILE ------------------------- */
