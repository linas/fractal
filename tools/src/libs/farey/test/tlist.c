
/* 
 * FUNCTION:
 * print t**k style farey numbers
 *
 * HISTORY:
 * Linas Vepstas January 16 1994
 * update July 1995 -- linas
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
   double x, y, t;
   int i, n;

   if (argc <3) {
      printf ("Usage: %s <number of terms> <t-value> \n", argv[0]);
      exit (1);
   }

   f = CreateFarey();

   n = atoi (argv[1]);
   t = atof (argv[2]);

   for (i=0; i<n; i++){
      x = ((double) (i+1))/ ((double) n);
      SetReal (f, x);
      y = ContinuedFractionToTFarey (f, t);

      printf ("i %g f %g \n", x, y);
      fflush (stdout);
   }

}

/* ---------------------- END OF FILE ------------------------- */
