
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

#define ZZZ 0.8
#define PTS 4060
#define DOME 0.0001

/* ------------------------------------------------------------ */

float cn_map (struct Farey *f, int nume, int deno, double omega) 
{
   int n;
   float y;

   RatioToContinuedFraction (f, nume, deno);
   y = 0.5 * ContinuedFractionToCnReal (f, omega); 
   /* y = 0.5 * ContinuedFractionToFCnReal (f, omega);  */
   /* y = 0.5 * ContinuedFractionToZCnReal (f, omega, ZZZ); */

   n = (int) y;
   if (y<0.0) n--;
   y -= (float) n;
   return y;
}

/* ------------------------------------------------------------ */

float sn_map (struct Farey *f, int nume, int deno, double omega) 
{
   int n;
   float y;

   RatioToContinuedFraction (f, nume, deno);
   y = 0.5 * ContinuedFractionToSnReal (f, omega); 
   /* y = 0.5 * ContinuedFractionToFSnReal (f, omega); */
   /* y = ContinuedFractionToZSnReal (f, omega, ZZZ); */

   n = (int) y;
   if (y<0.0) n--;
   y -= (float) n;
   return y;
}

/* ------------------------------------------------------------ */

float straight_cmap (struct Farey *f, int nume, int deno, double omega) 
{
   int n;
   float y;

   RatioToContinuedFraction (f, nume, deno);
   y = ContinuedFractionToCosReal (f, omega);
   n = (int) y;
   if (y<0.0) n--;
   y -= (float) n;
   return y;
}

/* ------------------------------------------------------------ */

float straight_zmap (struct Farey *f, int nume, int deno, double omega) 
{
   int n;
   float y;

   RatioToContinuedFraction (f, nume, deno);
   y = ContinuedFractionToZReal (f, omega);
   n = (int) y;
   if (y<0.0) n--;
   y -= (float) n;
   return y;
}

/* ------------------------------------------------------------ */
/* compute measure density in the range */

void splatter (struct Farey *f, float *arr, int width, int height, double omega)
{
   int i, j, n;
   int nume, deno;
   float x,y;
   int xbin, ybin ;

   deno = 61 * PTS +1;
   for (j=0; j<PTS; j++){ 
   /* for (j=PTS/2; j<2*PTS/3; j++){ */
   /* {j = PTS/2; */
      nume = 61*j+29; 
      if (0 == j%20) printf ("buzzing %d of %d \n", j, PTS);
/*
      x = straight_cmap (f, nume, deno, 0.3333*omega); 
      x = cn_map (f, nume, deno, omega); 
*/
      x = straight_zmap (f, nume, deno, omega); 

      x = 1.0 + cos (2.0 * M_PI * x);
      x *= 0.5;

      xbin = (int) (x * ((double) width));
      if (xbin>width-1) xbin = width-1;
      if (xbin<0) xbin = 0;

/*
      y = straight_cmap (f, nume, deno, omega); 
      y = sn_map (f, nume, deno, omega); 
*/
      y = straight_zmap (f, nume, deno, 1.0-omega); 

      y = 1.0 + sin (2.0 * M_PI * y);
      y *= 0.5;

      ybin = (int) (y * ((double) height));
      if (ybin>height-1) ybin = height-1;
      if (ybin<0) ybin = 0;


      /* printf ("got %f  %f %f \n", nume/(float)deno, x, y); */

      arr[ybin*width+xbin] += 1.0;
   }
}

/* ------------------------------------------------------------ */
/* compute measure density in the range */

void circle_splatter (struct Farey *f, float *arr, 
                     int width, int height, float param)
{
   int i, j, n;
   int nume, deno;
   float x,y;
   int xbin, ybin;
   double omega;

   deno = 61 * PTS +1;
   for (j=0; j<PTS; j++){ 
   /* for (j=PTS/2; j<2*PTS/3; j++){ */
   /* { j = 2*PTS/3-2; */
      nume = 61*j+29; 

printf ("horking %d of %d \n", j, PTS);

      for (omega = 0.0; omega < 2.0; omega +=DOME) {
/*
printf ("bzz %f \n", omega);
         x = straight_cmap (f, nume, deno, 0.3333*omega); 
         x = cn_map (f, nume, deno, omega); 
*/
         x = straight_zmap (f, nume, deno, omega); 

         xbin = (int) (x * ((double) width));
         if (xbin>width-1) xbin = width-1;
         if (xbin<0) xbin = 0;

/*
         y = straight_cmap (f, nume, deno, omega); 
         y = sn_map (f, nume, deno, omega); 
*/
         y = straight_zmap (f, nume, deno, param - omega); 

         ybin = (int) (y * ((double) height));
         if (ybin>height-1) ybin = height-1;
         if (ybin<0) ybin = 0;

         arr[ybin*width+xbin] += 1.0;
      }
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
   omega = atof (argv[4]);

   printf (" width=%d height=%d %s t=%f \n", width, height, argv[4], omega);

   pixel_map = (float *) malloc (width*height*sizeof (float));
   for (i=0; i<width*height; i++) pixel_map[i] = 0.00000000001;


   if (!strcmp ("cprod", argv[0])) splatter (f, pixel_map, width, height, omega);
   if (!strcmp ("cden", argv[0])) circle_splatter (f, pixel_map, width, height, omega);

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
