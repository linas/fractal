
/* 
 * FUNCTION:
 * generate pixmap of continued fraction data.
 *
 * HISTORY:
 * Linas Vepstas January 16 1994
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
#ifdef WIDE
   n /= 2;
   n *= 2;
#endif 
   y -= (float) n;
   return y;
}

/* ------------------------------------------------------------ */
/* compute measure density in the range */

void bin_sort (struct Farey *f, float *arr, int narr, int nume, int deno)
{
   int j, n, width;
   float y;
   int bin;
   double t;

   for (bin=0; bin<narr; bin++) arr[bin] = 0.0;

#define NUM_SAMPLES (10)
   width = NUM_SAMPLES * narr;

   for (j=0; j<width; j++){
      t = ((double) j) / ((double) width);
      y = straight_zmap (f, nume, deno, t); 
      bin = (int) (y * ((double) narr));
      if (bin>narr-1) bin = narr-1;
      if (bin<0) bin = 0;
      arr[bin] += 1.0 / ((double) NUM_SAMPLES);
   }

   /* wrap */
   for (bin=0; bin<narr; bin++) {
      y = arr[bin];
      n = (int) y;
      if (y<0.0) n--;
      y -= (float) n;
      arr[bin] = y;
   }
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
   fprintf (fil, "%s %s\n", argv[2], argv[3]);

   /* null terminated string */
   /*  { char zip=0; fwrite (zip, sizeof(char), 1, fil); } */
   fflush (fil);
   
   width = atoi (argv[2]);
   height = atoi (argv[3]);

   pixel_row = (float *) malloc (width*sizeof (float));

   deno = 61*height +1;
   for (i=0; i<height; i++) {
printf ("horking %d of %d\n", i, height);
   
      bin_sort (f, pixel_row, width, 61*i+27, deno);

      fwrite (pixel_row, sizeof(float), width, fil);
      fflush (fil);
   }
}

/* ---------------------- END OF FILE ------------------------- */
