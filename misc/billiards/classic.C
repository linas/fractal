//
// FILE:
// classic.C
//
// FUNCTION:
// Sinai's Billards ray tracer.  Sinai's Billiards is a well-known
// example of a chaotic system.
//
// HISTORY:
// Linas Vepstas
// November 2001


#include <stdio.h>
#include <stdlib.h>

#include "ray.h"
#include "vvector.h"

/* ==================================== */
// The SinaiView class generates the pretty pictures, converts 
// colors, etc.
//
// The Trace() method performs ray-tracing using the classic Sinai
// reflective boundary conditions.
//
// The TraceToroid() method performs ray-tracing using periodic
// toroidial boundary conditions.

class SinaiView
   : public SinaiBox, 
     public Pixels
{
   public:
      SinaiView (int, int);
      void Trace (void);
      void TraceToroid (void);
      void ColoredMirrors (void);
      void ToPixels (void);

   public:
      int nx;  // dimension of pixel grid
      int ny;

      double eye[3]; // position of eye

      double reflectivity;
      double density;

   private:
      SinaiRay *sr;
      
};

/* ==================================== */

SinaiView::SinaiView (int px, int py)
   : Pixels (px, py)
{
   int i, j;

   density = 0.0;
   reflectivity = 1.0;

   nx = px;
   ny = py;

   eye[0] = 0.0;
   eye[1] = 0.0;
   eye[2] = 5.0;

   sr = new SinaiRay [nx*ny];

   // project rays from eyepoint 
   for (i=0; i<nx; i++) {
      for (j=0; j<ny; j++) {
         double pixel[3];
         double dir[3];
         pixel[0] = 2.0 * (((double) i) + 0.51333)/ ((double) nx) - 1.0;
         pixel[1] = 2.0 * (((double) j) + 0.51666)/ ((double) ny) - 1.0;
         pixel[2] = 1.0;
         VEC_DIFF (dir, pixel, eye);
         sr[nx*j+i].Set (pixel, dir);
         sr[nx*j+i].last_wall = 4;

      }
   }
}

/* ==================================== */

void
SinaiView::Trace(void)
{
   int i;
   for (i=0; i<nx*ny; i++) 
   {
      SinaiBox::Trace (sr[i]);
   }
}

/* ==================================== */

void
SinaiView::TraceToroid(void)
{
   int i;
   for (i=0; i<nx*ny; i++) 
   {
      SinaiBox::TraceToroid (sr[i]);
   }
}

/* ==================================== */

void 
SinaiView::ToPixels (void)
{
   double absorbtivity = 1.0 - reflectivity;

   for (int i=0; i<nx*ny; i++)
   {
      double red = 0.0;
      double green = 0.0;
      double blue = 0.0;

      if (0 == sr[i].last_wall) red = 255.0;
      if (1 == sr[i].last_wall) { green = 255.0; blue = 255.0; }
      if (2 == sr[i].last_wall) green = 255.0;
      if (3 == sr[i].last_wall) { blue = 255.0; red = 255.0; }
      if (4 == sr[i].last_wall) blue = 255.0;
      if (5 == sr[i].last_wall) { red = 255.0; green = 255.0; }

#if 0
      int ib = 0;
      for (int iw=0; iw<6; iw++) ib += sr[i].bounces[iw];

      red *= exp (-absorbtivity * ib);
      green *= exp (-absorbtivity * ib);
      blue *= exp (-absorbtivity * ib);

      red *= exp (-density * sr[i].distance);
      green *= exp (-density * sr[i].distance);
      blue *= exp (-density * sr[i].distance);
#endif
      
      abgr[i] = 0xff & ((unsigned int) red);
      abgr[i] |= (0xff & ((unsigned int) green)) << 8;
      abgr[i] |= (0xff & ((unsigned int) blue)) << 16;
   }
}

void 
SinaiView::ColoredMirrors (void)
{
   double absorbtivity = 1.0 - reflectivity;

   for (int i=0; i<nx*ny; i++)
   {
      double red = 255.0;
      double green = 255.0;
      double blue = 255.0;

      red *= exp (-absorbtivity * sr[i].bounces[0]);
      red *= exp (-absorbtivity * sr[i].bounces[1]);
      red *= exp (-absorbtivity * sr[i].sphere_hits);
      // red *= exp (-density * sr[i].distance);

      green *= exp (-absorbtivity * sr[i].bounces[2]);
      green *= exp (-absorbtivity * sr[i].bounces[3]);
      green *= exp (-absorbtivity * sr[i].sphere_hits);
      // green *= exp (-density * sr[i].distance);

      blue *= exp (-absorbtivity * sr[i].bounces[4]);
      blue *= exp (-absorbtivity * sr[i].bounces[5]);
      blue *= exp (-absorbtivity * sr[i].sphere_hits);
      // blue *= exp (-density * sr[i].distance);
      
      abgr[i] = 0xff & ((unsigned int) red);
      abgr[i] |= (0xff & ((unsigned int) green)) << 8;
      abgr[i] |= (0xff & ((unsigned int) blue)) << 16;
   }
}

/* ==================================== */

main (int argc, char * argv[])
{
   SinaiView v (400,400);

   v.density = 0.005;
   v.density = 0.0;

   v.radius = 0.6;
   v.reflectivity = 0.9975;
   v.max_distance = 1280.0;

   v.reflectivity = 0.92;
   v.max_distance = 40.0;

   v.reflectivity = 0.995;
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
   v.reflectivity = 1.0 - 0.01 * 320.0 / maxdist;  // ad hoc lighting
   v.niterations = niter;
   v.max_manhattan = manhat;

   // v.Trace();
   v.TraceToroid();
   // v.ColoredMirrors();
   v.ToPixels ();
   v.WriteMTV ("junk.mtv");
}

/* ===================== end of file ====================== */
