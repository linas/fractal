
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
   int num, den;

   if (argc <2) {
      printf ("Usage: %s <number of terms> \n", argv[0]);
      exit (1);
   }

   f = CreateFarey ();

   n = atoi (argv[1]);

   /* den = n*6899+1; */
   den = n* 65536;
   for (i=0; i<n; i++){
      /* num = i*6899; */
      num = 65536*i + rand();
      
      x = ((double) (num))/ ((double) den);

      RatioToContinuedFraction (f, num, den);
      y = ContinuedFractionTerms (f);

      printf ("i %g f %g \n", x, y);
      fflush (stdout);
   }

}

/* ---------------------- END OF FILE ------------------------- */
