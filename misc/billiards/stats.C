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
		double MeanFreePath (void);

   public:
      int nx;  // dimension of pixel grid
      int ny;

      double eye[3]; // position of eye

   private:
      SinaiRay *sr;
      
};

/* ==================================== */

SinaiStats::SinaiStats (int px, int py)
{
   int i, j;

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

double
SinaiStats::MeanFreePath(void)
{
   int i;
	double mean_path = 0.0;
	
   for (i=0; i<nx*ny; i++) 
   {
		mean_path += sr[i].distance / ((double) sr[i].sphere_hits);
   }

   mean_path /= (double) nx*ny;

   return mean_path;
}


/* ==================================== */

main (int argc, char * argv[])
{
   SinaiStats v (400,400);


   v.radius = 0.6;

   if (3 > argc) {
      printf ("Usage: %s <radius> <maxhits>\n", argv[0]);
      exit (1);
   }

   double radius = atof (argv[1]);
   int maxhits = atoi (argv[2]);

   v.radius = radius;
   v.max_distance = 1.0e100;
   v.niterations = 1123123123;
   v.max_manhattan = 123456789;
   v.max_sphere_hits = maxhits;

   // v.Trace();
   v.TraceToroid();

   double fp = v.MeanFreePath ();

   printf ("free path=%f\n", fp);
}

/* ===================== end of file ====================== */
