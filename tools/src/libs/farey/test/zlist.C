
/* 
 * FUNCTION:
 * Test operation of farey number converter
 *
 * HISTORY:
 * Linas Vepstas Januery 16 1994
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"

/* ------------------------------------------------------------ */

main (int argc, char *argv[])
{
   ContinuedFraction f;
   double x, y, z, t;
   int i, n;
   int nume, deno;
   Complex ct, cy;

   if (argc <3) {
      printf ("Usage: %s <number of terms> <base> \n", argv[0]);
      exit (1);
   }

   n = atoi (argv[1]);
   t = atof (argv[2]);

   deno = 6899 * n +1;
   // deno = 65536 *n;
   for (i=0; i<n; i++){
      nume = 6889*i;
      // nume = 65536*i + rand();

      x = ((double) (nume))/ ((double) deno);
      /* f.SetReal (x); */

      f.SetRatio (nume, deno);
      y = f.ToZReal (t);
      ct = t;
      cy = f.cToZReal (ct);
      z = cy.real();

/*
      f.SetRatio (i, n);
      z = f.ToZReal (t);

      y -= z;
*/

      printf ("i %g f %g %g %g \n", x, y, z, y-z);
      fflush (stdout);
   }

}

/* ---------------------- END OF FILE ------------------------- */
