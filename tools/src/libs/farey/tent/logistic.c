
/* 
 * FUNCTION:
 * Test operation of farey number converter
 *
 * HISTORY:
 * Linas Vepstas Januery 16 1994
 */

#include "gl.h"
#include <stdio.h>
#include <math.h>

double omega = 0.0;
double scale = 1.0;
double offset = 0.0;

/* ------------------------------------------------------------ */
double logistic_map (lambda, x)
double lambda, x;
{
   double y;

   x = scale * (x-0.5) + 0.5 - omega;

/*
   x -= (double) ((int) x);
*/

   y = lambda * x * (1.0-x);
   y += offset;
   y -= (double) ((int) y);
   if (y<0.0) y+=1.0;

   return (y);
}
/* ------------------------------------------------------------ */
double bernoulli_map (omega, x)
double omega, x;
{
   double y;

   if (x<0.5) {
      y = 2.0*x;
   } else {
      y = 2.0 * x - 1.0;
   }  
   y += omega;
   y -= (double) ((int) y);
   return (y);
}

/* ------------------------------------------------------------ */
double tent_map (lambda, x)
double lambda, x;
{
   if (x<0.5) {
      return (2.0*lambda*x);
   } else {
      return (2.0*lambda*(1.0-x));
   }  
}
/* ------------------------------------------------------------ */


/*
#define LOW (0.495)
#define HIGH (0.55)
*/

#define LOW (0.0)
#define HIGH (6.0)

redraw ()
{
   int i, j;
   double lambda;
   double x;

/*
   clear ();
*/

   ortho2 (LOW, HIGH, 0.0, 1.0);

   for (lambda=LOW; lambda<HIGH; lambda += (HIGH-LOW)*0.002) {
/*
offset = 0.5 - 0.25 * lambda;
*/

      for (j=0; j<120; j++) {
         x = ((float) j) / 120.0;
         for (i=0; i<60; i++) {
            x = logistic_map (lambda, x);        
         }
         v2d (lambda, x);
      }
      gsync ();
   }
}

/* ------------------------------------------------------------ */

main (argc, argv)
int argc;
char *argv[];
{
   double r, x, y, z, t;
   int i, n;
   int nume, deno;
   int ix, iy;

   if (argc <2) {
      printf ("Usage: %s <omega> \n", argv[0]);
      exit (1);
   }

   prefsize (250, 250);
   winopen ();

   offset = atof (argv[1]);

   while (TRUE) {
      XEvent ev;
      XNextEvent (dpy, &ev);

      switch (ev.type) {
         case ConfigureNotify:
         case Expose:
            printf (" Expose \n");
            redraw ();
            break;
         default:
            break;
      }
   }

}

/* ---------------------- END OF FILE ------------------------- */
