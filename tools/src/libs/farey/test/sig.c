
#include <signal.h>
#include <math.h>

main () {
   double p, x, y;

   sigsetmask (0xffff);

   p = pow (1.1, 26000.0);

   x = 26000.0 * log (1.1);

   y = exp (x);
   
   x = log (1.0e300);

  printf (" %g %g %g \n", p, x, y);
}
