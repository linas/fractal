
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

float straight_zmap (struct Farey *f, int nume, int deno, double t) 
{
   int n;
   float y;
   RatioToContinuedFraction (f, nume, deno);
   y = ContinuedFractionToZReal (f, t);
   n = (int) y;
   if (y<0.0) n--;
   y -= (double) n;
   return y;
}

/* ------------------------------------------------------------ */

float symmetric_zmap (struct Farey *f, int nume, int deno, double t) 
{
   int n;
   float y, z;
   RatioToContinuedFraction (f, nume, deno);
   y = ContinuedFractionToZReal (f, t);

   RatioToContinuedFraction (f, deno-nume, deno);
   z = ContinuedFractionToZReal (f, t);

   /* y = 0.5 * (y+z-1.0); */
   y = 0.5 * (y-z);

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
   int nume, deno;
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

   t = 1.0;
   delta_t = - 2.0 / ((double) height);
   t -= delta_t;

   for (i=0; i<height; i++) {
      t += delta_t;
   
      deno = 6899 * width +1;
      for (j=0; j<width; j++){
         nume = 6889*j; 
         /* pixel_row[j] = straight_zmap (f, nume, deno, t); */
         pixel_row[j] = symmetric_zmap (f, nume, deno, t);
      }
      fwrite (pixel_row, sizeof(float), width, fil);
      fflush (fil);
   }
}

/* ---------------------- END OF FILE ------------------------- */
