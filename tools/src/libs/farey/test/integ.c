
/* 
 * FUNCTION:
 * Integrate farey numbers
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
   double sum = 0.0;
   double sumsum = 0.0;

   if (argc <2) {
      printf ("Usage: %s <number> \n", argv[0]);
      exit (1);
   }

   f = CreateFarey();

   n = atoi (argv[1]);

   for (i=0; i<n; i++){
      x = ((double) (i+1))/ ((double) n);
      SetReal (f, x);
      y = ContinuedFractionToFarey (f);
      sum += y / ((double) n);
      sumsum += sum / ((double) n);

      printf ("x %g y %g I %g In %g II %g IIn %g \n", 
           x, y, sum, 2.0 *sum/(x*x), sumsum, 6.0*sumsum/(x*x*x) );
      fflush (stdout);
   }

}

/* ---------------------- END OF FILE ------------------------- */
