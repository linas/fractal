
#include <math.h>

main ()
{
   double s, c, sn , cn, sm, cm, r;
   int i;

   s = sin (0.00123);
   c = cos (0.00123);
   sm = 0.0;
   cm = 1.0;

   for (i=0; i<10000; i++) {
      sn = c*sm + s*cm;
      cn = c*cm - s*sm;
      sm = sn;
      cm = cn;
      r = 1.0 - (cn*cn + sn*sn);
      printf ("i = %d c = %f s = %f r = %e \n", i, cn, sn, r);
   }
}

