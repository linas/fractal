
/* 
 * FUNCTION:
 * Test operation of farey number converter
 *
 * HISTORY:
 * Linas Vepstas Januery 16 1994
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"

/* ------------------------------------------------------------ */

main (int argc, char *argv[])
{
   ContinuedFraction f;
   double x, y, z, t;
   double v, vprev, gap;
   int i, j, n, m, d;
   int nume, deno;
   int bn, bd;
   double vhi;
   int p,q,r,s;
   int dn, dd;

   if (argc <2) {
      printf ("Usage: %s <n> \n", argv[0]);
      exit (1);
   }

   n = atoi (argv[1]);

   // the goasl here is to find the largest gap in the middle.
   // p/q < nume/deno < r/s
   p = 1;
   q = 5;
   r = 3;
   s = 5;
   dn = r*q - s*p;
   dd = s*q;

   deno = q*n*dd;
   

   vhi = -1.0;
   for (i = 1; i < n; i ++) 
   {
      nume = p*n*dd + i*q*dn;
      x = ((double) nume) / ((double) deno);
      f.SetRatio (nume, deno);
      v = f.ToZRealGap ();
      if (v > vhi) {
         vhi = v;
         bn = nume;
         bd = deno;
         printf ("best: x=%g n=%d d=%d v=%g\n", x, bn,bd, v);
      }
   }

}

/* ---------------------- END OF FILE ------------------------- */
