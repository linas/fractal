
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

#define PTS 590

/* ------------------------------------------------------------ */

double straight_zmap (struct Farey *f, int nume, int deno, double t) 
{
   int n;
   double y;

   RatioToContinuedFraction (f, nume, deno);
   y = ContinuedFractionToZReal (f, t);
   n = (int) y;
   if (y<0.0) n--;
   y -= (float) n;
   return y;
}

/* ------------------------------------------------------------ */
/* compute measure density in the range */

void circle_splatter (struct Farey *f, float *arr, 
                     int width, int height, float param)
{
   int i, j, n;
   int nume, deno;
   double x,y,fz;
   double s,c,co;
   int xbin, ybin;
   double omega;

   deno = 61 * PTS +1;
   for (j=0; j<PTS; j++){ 
      nume = 61*j+29; 


         x = ((double)nume) / ((double)deno);
         fz = straight_zmap (f, nume, deno, 0.000001); 

         y = (x*(1+fz)-fz)/(fz*fz);

printf ("horking %f f0=%f i=%f\n", x,fz,y);
   }
}

/* ------------------------------------------------------------ */

main (argc, argv)
int argc;
char *argv[];
{
   struct Farey *f;
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

   f = CreateFarey();

   width = atoi (argv[2]);
   height = atoi (argv[3]);
   // omega = atof (argv[4]);

   // printf (" width=%d height=%d %s t=%f \n", width, height, argv[4], omega);
   printf (" width=%d height=%d \n", width, height);

   pixel_map = (float *) malloc (width*height*sizeof (float));
   for (i=0; i<width*height; i++) pixel_map[i] = 0.00000000001;


   circle_splatter (f, pixel_map, width, height, omega);

   fil = fopen (argv[1], "w");
   fprintf (fil, "%s %s\n", argv[2], argv[3]);

   /* null terminated string */
   /*  { char zip=0; fwrite (zip, sizeof(char), 1, fil); } */
   fflush (fil);
   
   fwrite (pixel_map, sizeof(float), width*height, fil);
   fflush (fil);
   fclose (fil);

   free (pixel_map);

   exit (0);
}

/* ---------------------- END OF FILE ------------------------- */
