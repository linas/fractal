
/* 
 * FUNCTION:
 * Test operation of farey number converter
 *
 * HISTORY:
 * Linas Vepstas Januery 16 1994
 */

#include <stdio.h>
#include <math.h>

#include "Farey.h"

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

   f = CreateFarey();

   n = atoi (argv[1]);

   for (i=0; i<n; i++){
/*
   for (i=2*n/3-2; i<2*n/3+2; i++){
*/
      x = ((double) (i+1))/ ((double) n);
      SetReal (f, x);
      y = ContinuedFractionToXFarey (f, -1.0);

      printf ("i %g f %g \n", x, y);
      fflush (stdout);
   }

}

/* ---------------------- END OF FILE ------------------------- */
