
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
   double x, y, t;
   int i, n;

   if (argc <3) {
      printf ("Usage: %s <number of terms> <base> \n", argv[0]);
      exit (1);
   }

   f = (struct Farey *) malloc (sizeof (struct Farey));

   n = atoi (argv[1]);
   t = atof (argv[2]);

   for (i=0; i<n; i++){
      x = ((double) (i+1))/ ((double) n);
      SetReal (f, x);
      /* y = ContinuedFractionToEFarey (f, t); */
      y = ContinuedFractionToZReal (f, exp(-t));
      SetReal (f, y);
      /* y = ContinuedFractionToEFarey (f, (2.0*log(2.0)-t)); */
      y = ContinuedFractionToEFarey (f, t);

      printf ("i %g f %g \n", x, y);
      fflush (stdout);
   }

}

/* ---------------------- END OF FILE ------------------------- */
