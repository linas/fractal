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
//
// The Trace() method performs ray-tracing using the classic Sinai
// reflective boundary conditions.
//
// The TraceToroid() method performs ray-tracing using periodic
// toroidial boundary conditions.

class SinaiView
   : public SinaiBox
{
   public:
      SinaiView (int, int);
      void Trace (void);
      void TraceToroid (void);
      void TestPattern (void);
      void ColoredMirrors (void);
      void ToPixels (void);
      void WriteMTV (const char *filename);

   public:
      int nx;  // dimension of pixel grid
      int ny;

      double eye[3]; // position of eye

      double reflectivity;
      double density;

   private:
      SinaiRay *sr;
      
      unsigned int *abgr;
      unsigned int *pack;

};

/* ==================================== */

SinaiView::SinaiView (int px, int py)
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
         sr[i].Set (pixel, dir);
         sr[i].last_wall = 4;

      }
   }

   abgr = new unsigned int [nx*ny];
   pack = new unsigned int [3*nx*ny/4+1];
   for (i=0; i<nx*ny; i++) abgr[i] = 0;
   for (i=0; i<3*nx*ny/4; i++) pack[i] = 0;
}

/* ==================================== */

#define NNN (-1.0 + 100.0*DEGENERATE_TOLERANCE)
#define PPP (1.0 - 100.0*DEGENERATE_TOLERANCE)

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
SinaiView::TestPattern (void)
{
   int i,j;

   for (j=0; j<ny; j++) {
      for (i=0; i<nx; i++) {
         abgr[nx*j+i] = i%255 | (j%255 << 8);
      }
   }
}

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

void 
SinaiView::WriteMTV (const char * filename)
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
