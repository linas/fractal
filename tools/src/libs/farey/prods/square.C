
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

main (int argc, char *argv[])
{
   ContinuedFraction f;
   double x, y, z, omega;
   int i, j, n;
   int nume, deno;
   FILE *fil;
   int sizex, sizey;
   double width, height, delta;
   double re_start, im_start;
   double re_center, im_center;
   double re_position, im_position;
   double re, im;
   float *pixel_map;

   if (argc <4) {
      printf ("Usage: %s <filename> <width> <height> [<t>]\n", argv[0]);
      exit (1);
   }

   sizex = atoi (argv[2]);
   sizey = atoi (argv[3]);
   // omega = atof (argv[4]);

   // printf (" width=%d height=%d %s t=%f \n", sizex, sizey, argv[4], omega);
   printf (" width=%d height=%d \n", sizex, sizey);


   width = 1.0;
   re_center = 0.5;
   im_center = 0.5;
   height = width * ((double) sizey) / ((double) sizex);
   delta = width / (double) sizex;
   re_start = re_center -  width / 2.0;
   im_start = im_center + height / 2.0;
   
   pixel_map = (float *) malloc (sizex*sizey*sizeof (float));
   for (i=0; i<sizex*sizey; i++) pixel_map[i] = 0.0;

   im_position = im_start;
   for (i=0; i<sizey; i++) {
      if (i%10==0) printf(" start row %d\n", i);
      re_position = re_start;
      for (j=0; j<sizex; j++) {
         re = re_position;
         im = im_position;

         f.SetReal (re);
         re = f.ToFarey ();

         f.SetReal (im);
         im = f.ToFarey ();

         x = re*re - im*im;
         y = 2 * re * im;

         re = InvZReal (x, 1.0);
         im = InvZReal (y, 1.0);
      

         pixel_map [i*sizex +j] = sqrt (re*re + im*im);
         re_position += delta;
      }
      im_position -= delta;  /*top to bottom, not bottom to top */
   }

   fil = fopen (argv[1], "w");
   fprintf (fil, "%s %s\n", argv[2], argv[3]);

   /* null terminated string */
   /*  { char zip=0; fwrite (zip, sizeof(char), 1, fil); } */
   fflush (fil);
   
   fwrite (pixel_map, sizeof(float), sizex*sizey, fil);
   fflush (fil);
   fclose (fil);

   free (pixel_map);

   exit (0);
}

/* ---------------------- END OF FILE ------------------------- */
