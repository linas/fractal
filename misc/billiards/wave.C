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


#include <complex.h>
#include <math.h>
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

      double omega;  // energy
      int oversample;   // samples per pixel

      double shooter[3]; // position of source

      complex<double> *amplitude;
      int *norm;
};

/* ==================================== */

PathIntegral::PathIntegral (int px, int py)
   : Pixels (px, py)
{
   int i;

   oversample = 6;
   omega = 1.0;

   nx = px;
   ny = py;

   shooter[0] = 0.0;
   shooter[1] = 0.0;
   shooter[2] = 5.0;

   amplitude = new complex<double> [nx*ny];
   norm = new int [nx*ny];

   for (i=0; i<nx*ny; i++)
   {
      amplitude[i] = 0.0;
      norm[i] = 0;
   }

}

/* ==================================== */

complex<double> myexp (double const &ph)
{
    complex<double> retval (cos(ph), sin (ph));
    return retval;
}

/* ==================================== */
void
PathIntegral::TraceToroid(void)
{
   SinaiRay sr;
   int i,j;

   int ox = oversample *nx;
   int oy = oversample *ny;

   // project rays from source
   for (i=0; i<ox; i++) 
   {
      for (j=0; j<oy; j++) 
      {

         // define initial ray direction
         double pixel[3];
         double dir[3];
         pixel[0] = 2.0 * (((double) i) + 0.51333)/ ((double) ox) - 1.0;
         pixel[1] = 2.0 * (((double) j) + 0.51666)/ ((double) oy) - 1.0;
         pixel[2] = 1.0;
         VEC_DIFF (dir, pixel, shooter);
         sr.Init();
         sr.Set (pixel, dir);
         sr.last_wall = 4;

         // ray trace
         SinaiBox::TraceToroid (sr);

         // map final ray direction to pixel
         // do this by projection
         // treat the four side walls as identical,
         // ignore the front and back walls for now
         int use_ray = 0;
         if ((sr.direction[0] > 0.0)  && 
             (fabs (sr.direction[1]) < sr.direction[0]) &&
             (fabs (sr.direction[2]) < sr.direction[0]))
         {
            use_ray = 1;
         }
#define ROTATED_FACES
#ifdef ROTATED_FACES
         else // 90 degrees
         if ((sr.direction[1] > 0.0)  && 
             (fabs (sr.direction[0]) < sr.direction[1]) &&
             (fabs (sr.direction[2]) < sr.direction[1]))
         {
            double tmp = sr.direction[0];
            sr.direction[0] = sr.direction[1];
            sr.direction[1] = -tmp;
            use_ray = 1;
         }
         else   // 180 degrees
         if ((sr.direction[0] < 0.0)  && 
             (fabs (sr.direction[1]) < -sr.direction[0]) &&
             (fabs (sr.direction[2]) < -sr.direction[0]))
         {
            sr.direction[0] = - sr.direction[0];
            sr.direction[1] = - sr.direction[1];
            use_ray = 1;
         }
         else // 270 degrees
         if ((sr.direction[1] < 0.0)  && 
             (fabs (sr.direction[0]) < -sr.direction[1]) &&
             (fabs (sr.direction[2]) < -sr.direction[1]))
         {
            double tmp = sr.direction[1];
            sr.direction[1] = sr.direction[0];
            sr.direction[0] = -tmp;
            use_ray = 1;
         }
#endif

         if (use_ray)
         {

             // first convert ray direction to grid coords
             double x = sr.direction[1] / sr.direction[0];
             double y = sr.direction[2] / sr.direction[0];

             int px = (int) (((double) nx) * 0.5 * (x+1.0));
             int py = (int) (((double) ny) * 0.5 * (y+1.0));

             if ((0 > px) || (px >=nx) || (0 > py) || (py >=ny))
             {
                printf ("duude out of bounds !! %d %d \n", px, py);
             }

             // next perform phase summation
             double phase = omega * sr.distance;
             amplitude [nx*py+px] += myexp (phase);
             norm [nx*py+px] ++;
          
#if 0
printf ("duude ph= %f ", phase);
printf ("ampd= %f %f ", real(amplitude [nx*py+px]), imag(amplitude [nx*py+px]));
printf ("duude amp= %f\n", abs<double>(amplitude [nx*py+px]));
#endif
         }
      }

      if (0 == i % oversample) { printf ("."); fflush (stdout); }
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

      if (0 < norm[i])
      {
         red = 255.0 * abs<double> (amplitude[i]) / ((double) norm[i]);
      }

      abgr[i] = 0xff & ((unsigned int) red);
      abgr[i] |= (0xff & ((unsigned int) green)) << 8;
      abgr[i] |= (0xff & ((unsigned int) blue)) << 16;
   }
}

/* ==================================== */

main (int argc, char * argv[])
{
   PathIntegral v (400,400);

   if (6 > argc) {
      printf ("Usage: %s <fileout> <radius> <omega> <samples> <maxdist> [<niterations> [<max manhattan>]]\n", argv[0]);
      exit (1);
   }

   char * outfile = argv[1];
   double radius = atof (argv[2]);
   double omega = atof (argv[3]);
   int samples = atoi (argv[4]);
   double maxdist= atof (argv[5]);

   int niter = 1000000;
   if (7 == argc) niter = atoi (argv[6]);

   int manhat = 1000000;
   if (8 == argc) manhat = atoi (argv[7]);

   v.radius = radius;
   v.omega = omega;
   v.oversample = samples;
   v.max_distance = maxdist;
   v.niterations = niter;
   v.max_manhattan = manhat;

   // v.Trace();
   v.TraceToroid();
   v.ToPixels ();
   v.WriteMTV (outfile);
}

/* ===================== end of file ====================== */
