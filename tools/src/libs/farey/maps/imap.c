
/* 
 * FUNCTION:
 * generate pixmap of continued fraction data.
 *
 * HISTORY:
 * Linas Vepstas Januery 16 1994
 */

#include "Farey.h"
#include <stdio.h>
#include <math.h>

/* ------------------------------------------------------------ */
/* very crude inversion by search & previous guess */

float straight_zmap (struct Farey *f, double where, double t, float guess) 
{
   int n;
   float y;
   RealToContinuedFraction (f, guess);
   y = ContinuedFractionToZReal (f, t);

   if (where > y) {
      while (where > y) {
         guess += 0.001;
         RealToContinuedFraction (f, guess);
         y = ContinuedFractionToZReal (f, t);
      }
   } else {
      while (where < y) {
         guess -= 0.001;
         RealToContinuedFraction (f, guess);
         y = ContinuedFractionToZReal (f, t);
      }
   }

   n = (int) y;
   if (y<0.0) n--;
   y -= (double) n;
   return y;
}

/* ------------------------------------------------------------ */

main (argc, argv)
int argc;
char *argv[];
{
   struct Farey *f;
   double x, y, z, t;
   double delta_t;
   int i, j, n;
   FILE *fil;
   int width, height;
   float *pixel_row;

   if (argc <4) {
      printf ("Usage: %s <filename> <width> <height> \n", argv[0]);
      exit (1);
   }

   f = CreateFarey();

   fil = fopen (argv[1], "w");
   fprintf (fil, "%s %s", argv[2], argv[3]);

   /* null terminated string */
   { char zip=0; fwrite (zip, sizeof(char), 1, fil); }
   fflush (fil);
   
   width = atoi (argv[2]);
   height = atoi (argv[3]);

   pixel_row = (float *) malloc (width*sizeof (float));
   for (j=0; j<width; j++) pixel_row[j] = ((double) j) / ((double) width);

   t = 1.0;
   delta_t = - 2.0 / ((double) height);
   t -= delta_t;

   for (i=0; i<height; i++) {
      t += delta_t;
   
      for (j=0; j<width; j++){
         x = ((double) j) / ((double) width);
         pixel_row[j] = straight_zmap (f, x, t, pixel_row[j]);
      }
      fwrite (pixel_row, sizeof(float), width, fil);
      fflush (fil);
   }
}

/* ---------------------- END OF FILE ------------------------- */
