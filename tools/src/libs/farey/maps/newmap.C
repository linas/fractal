
/* 
 * FUNCTION:
 * generate pixmap of continued fraction data.
 *
 * HISTORY:
 * Linas Vepstas January 16 1994
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"

/* ------------------------------------------------------------ */

double straight_zmap (ContinuedFraction &f, int nume, int deno, double t) 
{
   int n;
   double x,y;
   Complex z;

   f.SetRatio (nume, deno);
   /* y = f.ToExpReal (t);  -- emap */
   /* y = f.ToEFarey (t);  -- fmap */
   /* y = f.ToInvZReal (t);  -- invfrac */

   // palms is just plain zreal
   // y = f.ToZReal (t); 

   // explore the unit circle. set up coords so that bottom 
   // of image is 0, top is +pi
   t = 0.5 * (1.0-t);
   z = polar (1.0, M_PI*t);
   z = f.cToZReal (z); 
   // y = z.real();
   // y = z.imag();
   // y = abs(z);
   y = arg (z) / M_PI;
   if (0.0 != t) y /= -t;
   y = y +1.0;


   // y *= ((double) nume) / ((double) deno);
   // y *= ((double) deno) / ((double) nume);
   // y *= t + (1.0-t)*(((double) nume) / ((double) deno));
   // y *= t + (1.0-t)*(((double) deno) / ((double) nume));
   // x = ((double) nume) / ((double) deno);
   // y = t*y + (1.0-t)*(x*(1.0+y) -y)/(y*y);
   /* y = f.ToEReal (t);  -- ereal */
   /* y = f.ToEFraction (t); -- efrac */

   n = (int) y;
   if (y<0.0) n--;
   y -= (float) n;
   return y;
}

/* ------------------------------------------------------------ */
/* compute measure density in the range */

void bin_sort (ContinuedFraction &f, float *arr, int narr, double t)
{
   int j, n, width;
   int nume, deno;
   float y;
   int bin;

   for (bin=0; bin<narr; bin++) arr[bin] = 0.0;

#define NUM_SAMPLES (76)
   width = NUM_SAMPLES * narr;

   deno = 61 * width +1;

   for (j=0; j<width; j++){
      nume = 61*j; 
      y = (float) straight_zmap (f, nume, deno, t);
      /* y = straight_zmap (f, deno, nume, t); */ 
      bin = (int) (y * ((double) narr));
      if (bin>narr-1) bin = narr-1;
      if (bin<0) bin = 0;
      arr[bin] += 1.0 / ((double) NUM_SAMPLES);
   }

   for (bin=0; bin<narr; bin++) {
      y = arr[bin];
      n = (int) y;
      if (0.0>y) n--;
      y -= (float) n;
      arr[bin] = y;
   }
}

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
   fprintf (fil, "%s %s\n", argv[2], argv[3]);

   /* null terminated string */
   /*  { char zip=0; fwrite (&zip, sizeof(char), 1, fil); } */
   fflush (fil);
   
   width = atoi (argv[2]);
   height = atoi (argv[3]);

   pixel_row = (float *) malloc (width*sizeof (float));

   t = -1.0;
   delta_t =  2.0;
   delta_t /= ((double) height);
   t -= delta_t;

   for (i=0; i<height; i++) {
      t += delta_t;
      if (0 == i%10) printf ("horking %d of %d \n", i, height);
   
      bin_sort (f, pixel_row, width, t);

      fwrite (pixel_row, sizeof(float), width, fil);
      fflush (fil);
   }
}

/* ---------------------- END OF FILE ------------------------- */
