//
// FILE:
// ray.C
//
// FUNCTION:
// Sinai's Billards ray tracer.  Sinai's Billiards is a well-known
// example of a chaotic system.
//
// Design:
// put sphere in center of cube.
// cube is 2 units on a side.
// sphere has radius < 1.0
//
// HISTORY:
// Linas Vepstas
// November 2001


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "vvector.h"
#include "intersect.h"

// Intersect() determines the intersection between the ray and the plane.
// It returns the distance to the intersection if the intersection
// is forward along the ray.  Otherwise it returns minus the distance.
//
// Bounce() determines if a forward trace of the ray will intersect the
// plane.  If it does, then the ray will be traced forward and reflected
// off the plane.  
//
// Forward() determines if a forward trace of the ray will intersect the
// plane.  If it does, then the ray will be traced forward to the plane.
//
// Sphere() returns a positive number if the ray intersects with a sphere
// of radius R centered at the origin; otherwise negative.  The returned
// value is the distance the ray travelled to hit the sphere.

class Ray
{
   public:
      void Set (double x[3], double t[3]);  // assign value
      void Set (double, double, double, double, double, double);
      double Intersect (Ray& plane);
      void Bounce (Ray& plane);
      void Forward (Ray& plane);
      double Sphere (double r);
   public:
      double position[3];   
      double direction[3];
};

/* ==================================== */

void 
Ray::Set (double x[3], double t[3])
{
   VEC_COPY(position, x);
   VEC_COPY(direction,  t);
   VEC_NORMALIZE(direction);
}

void
Ray::Set (double x0, double x1, double x2, 
          double t0, double t1, double t2)
{
   position[0] = x0;
   position[1] = x1;
   position[2] = x2;
   direction[0] = t0;
   direction[1] = t1;
   direction[2] = t2;
   VEC_NORMALIZE(direction);
}

/* ==================================== */

double 
Ray::Intersect (Ray& plane)
{
   double p2[3];
   VEC_SUM (p2, position, direction);

   int valid;
   double result[3];
   INTERSECT (valid, result, 
             plane.position, plane.direction,
             position, p2);

   if (FALSE == valid) return -1.0/DEGENERATE_TOLERANCE;
   
   // The 'direction' should be parallel to the result, and thus,
   // the dot product should be the length.
   
   double dot;
   VEC_DIFF (p2, result, position);
   VEC_DOT_PRODUCT (dot, p2, direction);

   return dot;
}

/* ==================================== */

void
Ray::Bounce (Ray& plane)
{
   double p2[3];
   VEC_SUM (p2, position, direction);

   int valid;
   double result[3];
   INTERSECT (valid, result, 
             plane.position, plane.direction,
             position, p2);

   if (FALSE == valid) return;
   VEC_COPY (position, result);
   VEC_REFLECT (direction, direction, plane.direction);  
}

/* ==================================== */

void
Ray::Forward (Ray& plane)
{
   double p2[3];
   VEC_SUM (p2, position, direction);

   int valid;
   double result[3];
   INTERSECT (valid, result, 
             plane.position, plane.direction,
             position, p2);

   if (FALSE == valid) return;

   VEC_COPY (position, result);
}

/* ==================================== */

double
Ray::Sphere (double radius)
{
   double b[3], pact, r2, s;
   double radial[3], dist, oldpos[3];

   VEC_PERP (b, position, direction);  
   VEC_LENGTH (pact, b);  // impact paramter

   r2 = radius*radius;
   if (r2 < pact) return -1.0;
   s = - sqrt (r2 - pact);

   VEC_COPY (oldpos, position);
   VEC_BLEND (radial, 1.0, b, s, direction);
   VEC_COPY (position, radial);
   VEC_NORMALIZE (radial);
   VEC_REFLECT (direction, direction, radial);

   VEC_DIFF (oldpos, position, oldpos);
   VEC_LENGTH (dist, oldpos);

   return dist;
}

/* ==================================== */
// class inherits from ray; collects statistics about
// the path of the ray as its traced, including distance 
// traveled, number of bounces, etc.

class SinaiRay
  : public Ray
{
   public:
      SinaiRay (void);
      double distance;  // distance travelled by ray   
      int last_wall;    // last wall of intersection
      int bounces[6];   // how many bounces have occured off each wall
      int sphere_hits;  // how many times ray bounced off sphere
};


SinaiRay::SinaiRay (void)
{
   sphere_hits = 0;
   last_wall = -1;
   distance = 0.0;
   for (int i=0; i<6; i++) bounces[i] = 0;
}

/* ==================================== */
// The class SinaiBox defines the boundry conditions, the max
// iterations, the shape of the actual box.
//
// The Trace() method performs ray-tracing using the classic Sinai
// reflective boundary conditions.
//
// The TraceToroid() method performs ray-tracing using periodic
// toroidial boundary conditions.

class SinaiBox
{
   public:
      SinaiBox (void);
      void Trace (SinaiRay &);
      void TraceToroid (SinaiRay &);
   public:

      int niterations;
      double max_distance;
      int max_manhattan;
      double radius;

   private:
      Ray walls[6]; // left, right, top, bottom, front, back;
};

SinaiBox::SinaiBox (void)
{
   radius = 0.5;
   niterations = 10000;
   max_manhattan = 10000;
   max_distance = 10000.0;

   walls[0].Set(-1.0, 0.0, 0.0, 1.0, 0.0, 0.0);
   walls[1].Set(1.0, 0.0, 0.0, -1.0, 0.0, 0.0);
   walls[2].Set(0.0, 1.0, 0.0, 0.0, -1.0, 0.0);
   walls[3].Set(0.0, -1.0, 0.0, 0.0, 1.0, 0.0);
   walls[4].Set(0.0, 0.0, 1.0, 0.0, 0.0, -1.0);
   walls[5].Set(0.0, 0.0, -1.0, 0.0, 0.0, 1.0);
}

/* ==================================== */

#define NNN (-1.0 + 100.0*DEGENERATE_TOLERANCE)
#define PPP (1.0 - 100.0*DEGENERATE_TOLERANCE)

void 
SinaiBox::Trace (SinaiRay &sr)
{
   for (int n=0; n<niterations; n++)
   {
      double dist = 0.0;

      // first, see if ray bounces off sphere
      dist = sr.Sphere (radius);
      if (0.0 < dist) 
      {
         sr.sphere_hits ++;
         sr.distance += dist;
         sr.last_wall = -1;
      }

      // next, trace ray to wall.
      int next_wall = -1;
      double nearest = 1000000.0;
      for (int iwall=0; iwall<6; iwall ++)
      {
         if (iwall == sr.last_wall) continue;

         double dist = sr.Intersect (walls[iwall]);
         if (0.0 > dist) continue;
         if (nearest > dist) { nearest = dist; next_wall = iwall; }
      }
      sr.bounces[next_wall] ++;
      sr.last_wall = next_wall;
      sr.distance += nearest;

      sr.Bounce (walls[next_wall]);

#if BROKEN
      // check for corner conditions
      if (((NNN >= sr.position[0]) ||
           (PPP <= sr.position[0])) &&
           (0 != sr.last_wall) &&
           (1 != sr.last_wall))
      {
         sr.direction[0] = -  sr[i].direction[0];
      }

      if (((NNN >= sr.position[1]) ||
           (PPP <= sr.position[1])) &&
           (2 != sr.last_wall) &&
           (3 != sr.last_wall))
      {
         sr.direction[1] = -  sr[i].direction[1];
      }

      if (((NNN >= sr.position[2]) ||
           (PPP <= sr.position[2])) &&
           (4 != sr.last_wall) &&
           (5 != sr.last_wall))
      {
         sr.direction[2] = -  sr[i].direction[2];
      }
#endif

      if (sr.distance > max_distance) break;
      if (sr.bounces[next_wall] >= max_manhattan) break;
   }
}

/* ==================================== */

void 
SinaiBox::TraceToroid (SinaiRay &sr)
{
   for (int n=0; n<niterations; n++)
   {
      double dist = 0.0;

      // first, see if ray bounces off sphere
      dist = sr.Sphere (radius);
      if (0.0 < dist) 
      {
         sr.sphere_hits ++;
         sr.distance += dist;
         sr.last_wall = -1;
      }

      // next, trace ray to wall.
      int next_wall = -1;
      double nearest = 1000000.0;
      for (int iwall=0; iwall<6; iwall ++)
      {
         if (iwall == sr.last_wall) continue;

         double dist = sr.Intersect (walls[iwall]);
         if (0.0 > dist) continue;
         if (nearest > dist) { nearest = dist; next_wall = iwall; }
      }
      sr.bounces[next_wall] ++;
      sr.last_wall = next_wall;
      sr.distance += nearest;

      sr.Forward (walls[next_wall]);

      // enforce periodic boundary conditions
      sr.position[0] += 2.0 * walls[next_wall].direction[0];
      sr.position[1] += 2.0 * walls[next_wall].direction[1];
      sr.position[2] += 2.0 * walls[next_wall].direction[2];

      // make sure we set the flags correctly
      switch (sr.last_wall) {
         case 0: sr.last_wall = 1; break;
         case 1: sr.last_wall = 0; break;
         case 2: sr.last_wall = 3; break;
         case 3: sr.last_wall = 2; break;
         case 4: sr.last_wall = 5; break;
         case 5: sr.last_wall = 4; break;
         default: break;
      }
         
      if (sr.distance > max_distance) break;
      if (sr.bounces[next_wall] >= max_manhattan) break;
   }
}

/* ==================================== */
// The Pixels class handles colors
//

class Pixels
{
   public:
      Pixels (int, int);
      void TestPattern (void);
      void WriteMTV (const char *filename);

   public:
      int nx;  // dimension of pixel grid
      int ny;

      unsigned int *abgr;
   private:
      unsigned int *pack;
};

/* ==================================== */

Pixels::Pixels (int px, int py)
{
   int i;

   nx = px;
   ny = py;

   abgr = new unsigned int [nx*ny];
   pack = new unsigned int [3*nx*ny/4+1];
   for (i=0; i<nx*ny; i++) abgr[i] = 0;
   for (i=0; i<3*nx*ny/4; i++) pack[i] = 0;
}

void 
Pixels::TestPattern (void)
{
   int i,j;

   for (j=0; j<ny; j++) {
      for (i=0; i<nx; i++) {
         abgr[nx*j+i] = i%255 | (j%255 << 8);
      }
   }
}

/* ==================================== */

void 
Pixels::WriteMTV (const char * filename)
{
   // first, pack the pixels byte by byte
   unsigned int *p = pack;
   for (int i=0; i<nx*ny; i++) {
      *p = abgr[i] & 0xffffff;
      ((char *) p) += 3;
   }

   // write out MTV format pixmap.
   FILE *fh = fopen (filename, "w");
   fprintf (fh, "%d %d\n", nx, ny);
   fwrite (pack, 4, 3*nx*ny/4+1, fh);
   fclose (fh);
}

/* ===================== end of file ====================== */
