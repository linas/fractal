
/* 
 * FUNCTION:
 * generate pixmap of continued fraction data.
 *
 * HISTORY:
 * Linas Vepstas January 16 1994
 * Updates Linas March, April 1996
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Farey.h"

/* ------------------------------------------------------------ */
/* compute measure density in the range */

void splatter (float *arr, int width, int height, double omega)
{
   int i, j, n;
   float x,y;
   int xbin, ybin ;
   double r, zee;

#define PTS 860

   zee = 0.1;
   for (j=0; j<PTS; j++){ 
      r = ((double) j) / ((double) PTS);

      x = r;

      xbin = (int) (x * ((double) width));
      if (xbin>width-1) xbin = width-1;
      if (xbin<0) xbin = 0;

      y = InvZReal (r, zee);

      ybin = (int) (y * ((double) height));
      if (ybin>height-1) ybin = height-1;
      if (ybin<0) ybin = 0;


      arr[ybin*width+xbin] += 1.0;
   }
}

/* ------------------------------------------------------------ */

main (argc, argv)
int argc;
char *argv[];
{
   double x, y, z, omega;
   double delta_omega;
   int i, j, n;
   int nume, deno;
   FILE *fil;
   int width, height;
   float *pixel_map;

   if (argc <4) {
      printf ("Usage: %s <filename> <width> <height> [<t>]\n", argv[0]);
      exit (1);
   }

   width = atoi (argv[2]);
   height = atoi (argv[3]);
   omega = atof (argv[4]);

   printf (" width=%d height=%d %s t=%f \n", width, height, argv[4], omega);

   pixel_map = (float *) malloc (width*height*sizeof (float));
   for (i=0; i<width*height; i++) pixel_map[i] = 0.00000000001;


   splatter (pixel_map, width, height, omega);

/*
   if (!strcmp ("cprod", argv[0])) splatter (f, pixel_map, width, height, omega);
   if (!strcmp ("cden", argv[0])) {
      for (omega=0.0; omega <6.28; omega +=0.003) {
         printf ("horking %f \n", omega);
         splatter (f, pixel_map, width, height, omega);
      }
   }
*/

   fil = fopen (argv[1], "w");
   fprintf (fil, "%s %s\n", argv[2], argv[3]);

   /* null terminated string */
   /*  { char zip=0; fwrite (zip, sizeof(char), 1, fil); } */
   fflush (fil);
   
   fwrite (pixel_map, sizeof(float), width*height, fil);
   fflush (fil);
   fclose (fil);

   free (pixel_map);
}

/* ---------------------- END OF FILE ------------------------- */
