//
// FILE:
// wave.C
//
// FUNCTION:
// Sinai's Billards ray tracer.  Sinai's Billiards is a well-known
// example of a chaotic system.
//
// This has a twist: we do a feynmann path integral over rays.
//
// HISTORY:
// Linas Vepstas
// November 2001


#include <stdio.h>
#include <stdlib.h>

#include "ray.h"
#include "vvector.h"

/* ==================================== */
// The PathIntegral class generates the pretty pictures, converts 
// colors, etc.
//
// The Trace() method performs ray-tracing using the classic Sinai
// reflective boundary conditions.
//
// The TraceToroid() method performs ray-tracing using periodic
// toroidial boundary conditions.

class PathIntegral
   : public SinaiBox, 
     public Pixels
{
   public:
      PathIntegral (int, int);
      void Trace (void);
      void TraceToroid (void);
      void ToPixels (void);

   public:
      int nx;  // dimension of pixel grid
      int ny;

      double shooter[3]; // position of source

      double *phase;
      double omega;  // energy
      int oversample;   // samples per pixel
};

/* ==================================== */

PathIntegral::PathIntegral (int px, int py)
   : Pixels (px, py)
{
   int i;

   nx = px;
   ny = py;

   shooter[0] = 0.0;
   shooter[1] = 0.0;
   shooter[2] = 5.0;

   phase = new double [nx*ny];

   for (i=0; i<nx*ny; i++) phase[i] = 0;

   oversample = 2;
}

/* ==================================== */

void
PathIntegral::TraceToroid(void)
{
   SinaiRay sr;
   int i,j;

   int px = oversample *nx;
   int py = oversample *ny;

   // project rays from source
   for (i=0; i<px; i++) 
   {
      for (j=0; j<py; j++) 
      {

         // define initial ray direction
         double pixel[3];
         double dir[3];
         pixel[0] = 2.0 * (((double) i) + 0.51333)/ ((double) px) - 1.0;
         pixel[1] = 2.0 * (((double) j) + 0.51666)/ ((double) py) - 1.0;
         pixel[2] = 1.0;
         VEC_DIFF (dir, pixel, shooter);
         sr.Set (pixel, dir);
         sr.last_wall = 4;

         // ray trace
         SinaiBox::TraceToroid (sr);

         // map final ray direction to pixel
         // do this by projection
         // treat the four side walls as identical,
         // ignore the front and back walls for now
         xxxxxxxxxxxx
         use_ray = 0;
         if ((sr.direction[0] > 0.0)  && 
             (fabs (sr.direction[1]) < sr.direction[0]) &&
             (fabs (sr.direction[2]) < sr.direction[0]))
         {
            use_ray = 1;
         }

         if (use_ray)
         {
             double x = sr.direction[1] / sr.direction[0];
             double y = sr.direction[2] / sr.direction[0];
printf ("duude xy= %f %f \n", x, y);
         }
      }
   }

}

/* ==================================== */

void 
PathIntegral::ToPixels (void)
{

   for (int i=0; i<nx*ny; i++)
   {
      double red = 0.0;
      double green = 0.0;
      double blue = 0.0;

      abgr[i] = 0xff & ((unsigned int) red);
      abgr[i] |= (0xff & ((unsigned int) green)) << 8;
      abgr[i] |= (0xff & ((unsigned int) blue)) << 16;
   }
}

/* ==================================== */

main (int argc, char * argv[])
{
   PathIntegral v (400,400);

   v.radius = 0.6;
   v.max_distance = 640.0;

   if (3 > argc) {
      printf ("Usage: %s <radius> <maxdist> [<niterations> [<max manhattan>]]\n", argv[0]);
      exit (1);
   }

   double radius = atof (argv[1]);
   double maxdist= atof (argv[2]);

   int niter = 1000000;
   if (4 == argc) niter = atoi (argv[3]);

   int manhat = 1000000;
   if (5 == argc) manhat = atoi (argv[4]);

   v.radius = radius;
   v.max_distance = maxdist;
   v.niterations = niter;
   v.max_manhattan = manhat;

   // v.Trace();
   v.TraceToroid();
   v.ToPixels ();
   v.WriteMTV ("junk.mtv");
}

/* ===================== end of file ====================== */
