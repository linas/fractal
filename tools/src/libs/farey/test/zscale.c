
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

   deno = n* 6899 +1;
   for (i=0; i<n; i++){
      nume = 6899*i;

      x = ((double) (nume))/ ((double) deno);

#ifdef FLOAT_PT_COMPS
      x = ((double) (i))/ ((double) n);
#ifdef WORKS
      y = 1.0 / (2.0 -x);
      SetReal (f, y);
      y = ContinuedFractionToZReal (f, t);
      y = 1.0 / y;
      y = 2.0 - y;
/*
      y -= 0.8;
      y *= 5.0;
*/
      y -= (1.0-t);
      y /= t;
      x = 1.0 - x;
      y = 1.0 - y;

      SetReal (f, x);
      z = ContinuedFractionToZReal (f, t);
      y -= z;
#endif

      y = 1.0 / (1+x);
      SetReal (f, y);
      z = ContinuedFractionToZReal (f, t);

      SetReal (f, x);
      y = ContinuedFractionToZReal (f, t);
      y = 1.0 / (1.0 + t*y);

      y -= z;
#endif

      RatioToContinuedFraction (f, deno, deno+nume);
      z = ContinuedFractionToZReal (f, t);

      RatioToContinuedFraction (f, nume, deno);
      y = ContinuedFractionToZReal (f, t);
      y = 1.0 / (1.0 + t*y); 

      y -= z;
      

      printf ("i %g f %g \n", x, y);
      fflush (stdout);
   }

}

/* ---------------------- END OF FILE ------------------------- */
