
/* 
 * HISTORY:
 * Linas Vepstas February 1994
 */

#include <math.h>
#include <gl/gl.h>
#include "Farey.h"

#define move2d(a,b)  move2 ((float) (a), (float) (b));
#define draw2d(a,b)  draw2 ((float) (a), (float) (b));
#define pnt2d(a,b)  pnt2 ((float) (a), (float) (b));

main () {
   int i, length;
   int nume, deno;
   struct Farey *f;
   double z, zlast;
   double r, theta;
   double x, y;
   double ratio;

   prefsize (400, 400);
   winopen ("Farey Numbers");
   color (BLACK);
   clear ();
   color (WHITE);
   ortho2 (-3.0, 3.0, -3.0, 3.0);

   f = (struct Farey *) malloc (sizeof (struct Farey));

   length = 30000;

   deno = length * 367 +1;
   zlast = 0.0;
   move2d (0.0, 0.0);

   for (i=0; i< length; i++) {
      nume = i * 367 +1;
      ratio = ((double) nume) / ((double) deno);
      RatioToContinuedFraction (f, nume, deno);

      z = GetFarey (f);
      theta = z;
      r = z - zlast;
      r *= ((double) length);
      r += 1.0;
      r = log (r);
      zlast = z;

      theta *= 2.0 * M_PI;

      x = r * cos (theta);
      y = r * sin (theta);

/*
      draw2d (x, y);
*/
      pnt2d (x, y);
  }

   sleep (1000);
}
