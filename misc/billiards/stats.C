//
// FILE:
// stats.C
//
// FUNCTION:
// Sinai's Billards ray tracer.  Sinai's Billiards is a well-known
// example of a chaotic system.  This tracer performs a 'classical
// mechanics' newtonian trace of rays, and colors resulting pixels
// appropriately.  It collects statistics about the trace.
//
// HISTORY:
// Linas Vepstas
// March 2002


#include <stdio.h>
#include <stdlib.h>

#include "ray.h"
#include "vvector.h"

/* ==================================== */
//
// The Trace() method performs ray-tracing using the classic Sinai
// reflective boundary conditions.
//
// The TraceToroid() method performs ray-tracing using periodic
// toroidial boundary conditions.

class SinaiStats
   : public SinaiBox 
{
   public:
      SinaiStats (int, int);
      void Trace (void);
      void TraceToroid (void);

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

SinaiStats::SinaiStats (int px, int py)
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
SinaiStats::Trace(void)
{
   int i;
   for (i=0; i<nx*ny; i++) 
   {
      SinaiBox::Trace (sr[i]);
   }
}

/* ==================================== */

void
SinaiStats::TraceToroid(void)
{
   int i;
   for (i=0; i<nx*ny; i++) 
   {
      SinaiBox::TraceToroid (sr[i]);
   }
}


/* ==================================== */

main (int argc, char * argv[])
{
   SinaiStats v (400,400);


   v.radius = 0.6;
   v.reflectivity = 0.9975;
   v.max_distance = 1280.0;

   v.max_distance = 40.0;

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
   // v.ColoredMirrors();
}

/* ===================== end of file ====================== */
