
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
   { char zip=0; fwrite (&zip, sizeof(char), 1, fil); }
   fflush (fil);
   
   width = atoi (argv[2]);
   height = atoi (argv[3]);

   pixel_row = (float *) malloc (width*sizeof (float));

   t = -1.0;
   delta_t = 2.0;
   delta_t /= ((double) height);
   t -= delta_t;

   for (i=0; i<height; i++) {
      t += delta_t;
      printf ("horking %d of %d \n", i, height);
   
      for (j=0; j<width; j++)pixel_row[j] = 0.0;

#define NUMITER (3)
#define IRR 689 
      deno = IRR * NUMITER * width +1;
      for (j=0; j<NUMITER* width; j++){
         nume = IRR*j + IRR/2; 
   

/*
if (0.0 > t) {
*/
         RatioToContinuedFraction (f, nume, deno);
         y = ContinuedFractionToZReal (f, 1.5 - t);
         RatioToContinuedFraction (f, deno - nume, deno);
         y += 1.0 - ContinuedFractionToZReal (f, 1.5 - t);
         y *= 0.5;
/*
} else {
   y = nume / ((double) deno);
}

         RatioToContinuedFraction (f, nume, deno);
         y = ContinuedFractionToZReal (f, t) -0.5;
         RatioToContinuedFraction (f, deno - nume, deno);
         y *= 0.5 - ContinuedFractionToZReal (f, t);
         y *= 4.0;
*/

         n = (int) y;
         if (0.0 > y) n--;
         y -= (double) n;

         x = y;
         RealToContinuedFraction (f, x);
         y = ContinuedFractionToZReal (f, t);

         RealToContinuedFraction (f, 1.0 - x);
         y += 1.0 - ContinuedFractionToZReal (f, t);
         y *= 0.5;

         n = (int) y;
         if (0.0 > y) n--;
         y -= (double) n;

         y *= width;
         n = (int) y;
         if (0>n) n=0;
         if (width <=n) n = width-1;
         pixel_row[n] += 1.0/((double) NUMITER);

      }
      fwrite (pixel_row, sizeof(float), width, fil);
      fflush (fil);
   }

   fclose (fil);

   free (pixel_row);
}

/* ---------------------- END OF FILE ------------------------- */
