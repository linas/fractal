
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

float straight_cmap (struct Farey *f, int nume, int deno, double nu) 
{
   int n;
   float y;

   RatioToContinuedFraction (f, nume, deno);
   /* y = ContinuedFractionToSincReal (f, nu);  -- jmap */
   /* y = ContinuedFractionToCosReal (f, nu); -- cmap */

   y = ContinuedFractionToCnReal (f, nu);

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

void bin_sort (struct Farey *f, float *arr, int narr, double nu)
{
   int j, n, width;
   int nume, deno;
   float y;
   int bin;

   for (bin=0; bin<narr; bin++) arr[bin] = 0.0;

#define NUM_SAMPLES (10)
   width = NUM_SAMPLES * narr;

   deno = 61 * width +1;
#ifdef WIDE
   deno /= 2;
#endif
   for (j=0; j<width; j++){
      nume = 61*j; 
      y = straight_cmap (f, nume, deno, nu); 
      bin = (int) (y * ((double) narr));
      if (bin>narr-1) bin = narr-1;
      if (bin<0) bin = 0;
      arr[bin] += 1.0 / ((double) NUM_SAMPLES);
   }

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
   double x, y, z, nu;
   double delta_nu;
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

   nu = 0.0;
   delta_nu = 1.0;
   delta_nu /= ((double) height);
   nu -= delta_nu;

   for (i=0; i<height; i++) {
      nu += delta_nu;
      printf ("Horking %i of %i \n", i, height);
   
      bin_sort (f, pixel_row, width, nu);

      fwrite (pixel_row, sizeof(float), width, fil);
      fflush (fil);
   }
}

/* ---------------------- END OF FILE ------------------------- */
