
/* 
 * FUNCTION:
 * Test operation of farey number converter
 *
 * HISTORY:
 * Linas Vepstas January 16 1994
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
   double x, t;
   int n, d;
   int i;
   double delta_t;

   if (argc <3) {
      printf ("Usage: %s <num> <denom> \n", argv[0]);
      exit (1);
   }

   f = CreateFarey ();

   n = atoi (argv[1]);
   d = atoi (argv[2]);
   /* SetReal (f, x); */
   /* x = GetFarey (f); */
d=61*3200+1;
t = -0.28812499999997809;
t = 1.0 - 2.0 * (4060.0 / 5100.0);

t = 1.0;
delta_t = -2.0 / 5100.0;
t -= delta_t;

for(i=0; i<2061; i++) { t += delta_t;}

for (i=0; i<3200; i++) {
n = 61*i;
   RatioToContinuedFraction (f, n, d);
   /* PrintContinuedFraction (f); */
   /* x = GetReal (f); */
   /* x = ContinuedFractionToEFraction (f, 1.0); */
   /* x = ContinuedFractionToEFarey (f, log (2.0)); */
   x = ContinuedFractionToZReal (f, t);
}

   /* printf ("---------- its %f \n", x); */

}

/* ---------------------- END OF FILE ------------------------- */
