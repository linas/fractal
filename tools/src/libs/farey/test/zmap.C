
/* 
 * FUNCTION:
 * generate pixmap of continued fraction data.
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
   double delta_t;
   int i, j, n;
   int nume, deno;
   FILE *fil;
   int width, height;
   float *pixel_row;

   if (argc <4) {
      printf ("Usage: %s <filename> <width> <height> \n", argv[0]);
      exit (1);
   }

   fil = fopen (argv[1], "w");
   fprintf (fil, "%s %s", argv[2], argv[3]);

   /* null terminated string */
   { char zip=0; fwrite (&zip, sizeof(char), 1, fil); }
   fflush (fil);
   
   width = atoi (argv[2]);
   height = atoi (argv[3]);

   pixel_row = (float *) malloc (width*sizeof (float));

   t = 1.0;
   delta_t = - 2.0 / ((double) height);
   t -= delta_t;

   for (i=0; i<height; i++) {
      t += delta_t;
   
      deno = 6899 * width +1;
      for (j=0; j<width; j++){
         nume = 6889*i; 
   
         f.SetRatio (nume, deno);
         y = f.ToZReal (t);
         n = (int) y;
         y -= (double) n;
         pixel_row[i] = y;

      }
      fwrite (pixel_row, sizeof(float), width, fil);
      fflush (fil);
   }
}

/* ---------------------- END OF FILE ------------------------- */
