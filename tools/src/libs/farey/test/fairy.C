
/* 
 * FUNCTION:
 * Test operation of farey number converter
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

   if (argc <3) {
      printf ("Usage: %s <num> <denom> \n", argv[0]);
      exit (1);
   }

   n = atoi (argv[1]);
   d = atoi (argv[2]);

   t = ((double) n) / ((double) d);
   /* f.SetReal (t); */
   f.SetRatio (n, d);

   x = f.ToFarey (); 

   f.SetRatio (n, n+d);

   y = 2.0 * f.ToFarey (); 

   printf ("---------- f(%d/%d) = %f %f\n", n,d,x, y); 

}

/* ---------------------- END OF FILE ------------------------- */

